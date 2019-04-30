# code based on https://github.com/keras-team/keras-contrib/blob/master/examples/improved_wgan.py
# algorithm based on https://github.com/luslab/scRNAseq-WGAN-GP/blob/master/scripts/WGAN-GP_minimal.py

"""An implementation of the improved WGAN described in https://arxiv.org/abs/1704.00028

The improved WGAN has a term in the loss function which penalizes the network if its
gradient norm moves away from 1. This is included because the Earth Mover (EM) distance
used in WGANs is only easy to calculate for 1-Lipschitz functions (i.e. functions where
the gradient norm has a constant upper bound of 1).

The original WGAN paper enforced this by clipping weights to very small values
[-0.01, 0.01]. However, this drastically reduced network capacity. Penalizing the
gradient norm is more natural, but this requires second-order gradients. These are not
supported for some tensorflow ops (particularly MaxPool and AveragePool) in the current
release (1.0.x), but they are supported in the current nightly builds
(1.1.0-rc1 and higher).

To avoid this, this model uses strided convolutions instead of Average/Maxpooling for
downsampling. If you wish to use pooling operations in your discriminator, please ensure
you update Tensorflow to 1.1.0-rc1 or higher. I haven't tested this with Theano at all.

The model saves images using pillow. If you don't have pillow, either install it or
remove the calls to generate_images.
"""
import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
from PIL import Image
import argparse
import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from keras.models import Model, Sequential, load_model
from keras.layers import Input, Dense, Reshape, Flatten
from keras.layers.merge import _Merge
from keras.layers.convolutional import Convolution2D, Conv2DTranspose
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import LeakyReLU
from keras.optimizers import Adam
from keras.datasets import mnist
from keras import backend as K
from functools import partial
import sys
import tensorflow as tf
from keras.backend.tensorflow_backend import set_session

config = tf.ConfigProto()
config.gpu_options.allow_growth = True  # dynamically grow the memory used on the GPU
sess = tf.Session(config=config)
set_session(sess)  # set this TensorFlow session as the default session for Keras

# The training ratio is the number of discriminator updates
# per generator update. The paper uses 5.
GRADIENT_PENALTY_WEIGHT = 10  # As per the paper

# Deterministic output.
# Tired of seeing the same results every time? Remove the line below.
np.random.seed(1000)

inflate_to_size = 600

disc_internal_size = 200


def load_data(input_file):
    df = pd.read_csv(input_file)
    print(df.head())
    X_train = df.values
    # normalize it to the range [0, 1]
    X_train = X_train.astype(np.float32)
    X_case_control = X_train[:,-1:]
    min_max_scaler = MinMaxScaler()
    X_expression_scaled = min_max_scaler.fit_transform(X_train[:,0:-1])
    X_train_preprocessed = np.concatenate((X_expression_scaled, X_case_control), axis=1)
    return X_train_preprocessed


def wasserstein_loss(y_true, y_pred):
    """Calculates the Wasserstein loss for a sample batch.

    The Wasserstein loss function is very simple to calculate. In a standard GAN, the
    discriminator has a sigmoid output, representing the probability that samples are
    real or generated. In Wasserstein GANs, however, the output is linear with no
    activation function! Instead of being constrained to [0, 1], the discriminator wants
    to make the distance between its output for real and generated samples as
    large as possible.

    The most natural way to achieve this is to label generated samples -1 and real
    samples 1, instead of the 0 and 1 used in normal GANs, so that multiplying the
    outputs by the labels will give you the loss immediately.

    Note that the nature of this loss means that it can be (and frequently will be)
    less than 0."""
    return K.mean(y_true * y_pred)


def gradient_penalty_loss(y_true, y_pred, averaged_samples,
                          gradient_penalty_weight):
    """Calculates the gradient penalty loss for a batch of "averaged" samples.

    In Improved WGANs, the 1-Lipschitz constraint is enforced by adding a term to the
    loss function that penalizes the network if the gradient norm moves away from 1.
    However, it is impossible to evaluate this function at all points in the input
    space. The compromise used in the paper is to choose random points on the lines
    between real and generated samples, and check the gradients at these points. Note
    that it is the gradient w.r.t. the input averaged samples, not the weights of the
    discriminator, that we're penalizing!

    In order to evaluate the gradients, we must first run samples through the generator
    and evaluate the loss. Then we get the gradients of the discriminator w.r.t. the
    input averaged samples. The l2 norm and penalty can then be calculated for this
    gradient.

    Note that this loss function requires the original averaged samples as input, but
    Keras only supports passing y_true and y_pred to loss functions. To get around this,
    we make a partial() of the function with the averaged_samples argument, and use that
    for model training."""
    # first get the gradients:
    #   assuming: - that y_pred has dimensions (batch_size, 1)
    #             - averaged_samples has dimensions (batch_size, nbr_features)
    # gradients afterwards has dimension (batch_size, nbr_features), basically
    # a list of nbr_features-dimensional gradient vectors
    gradients = K.gradients(y_pred, averaged_samples)[0]
    # compute the euclidean norm by squaring ...
    gradients_sqr = K.square(gradients)
    #   ... summing over the rows ...
    gradients_sqr_sum = K.sum(gradients_sqr,
                              axis=np.arange(1, len(gradients_sqr.shape)))
    #   ... and sqrt
    gradient_l2_norm = K.sqrt(gradients_sqr_sum)
    # compute lambda * (1 - ||grad||)^2 still for each single sample
    gradient_penalty = gradient_penalty_weight * K.square(1 - gradient_l2_norm)
    # return the mean as loss over all the batch samples
    return K.mean(gradient_penalty)


def make_generator(input_dimension, random_dim):
    """Creates a generator model that takes a 100-dimensional noise vector as a "seed",
    and outputs images of size 28x28x1."""
    inputs = Input(shape=(random_dim,))
    x = Dense(inflate_to_size)(inputs)
    x = LeakyReLU()(x)
    x = Dense(inflate_to_size)(x)
    x = LeakyReLU()(x)
    x = Dense(input_dimension)(x)
    outputs = LeakyReLU()(x)
    model = Model(inputs, outputs)
    return model


def make_discriminator(input_dimension):
    """Creates a discriminator model that takes an image as input and outputs a single
    value, representing whether the input is real or generated. Unlike normal GANs, the
    output is not sigmoid and does not represent a probability! Instead, the output
    should be as large and negative as possible for generated inputs and as large and
    positive as possible for real inputs.

    Note that the improved WGAN paper suggests that BatchNormalization should not be
    used in the discriminator."""
    inputs = Input(shape=(input_dimension,))
    x = Dense(disc_internal_size)(inputs)
    x = LeakyReLU()(x)
    x = Dense(disc_internal_size)(x)
    x = LeakyReLU()(x)
    outputs = Dense(1)(x)
    model = Model(inputs, outputs)
    return model


class RandomWeightedAverage(_Merge):
    """Takes a randomly-weighted average of two tensors. In geometric terms, this
    outputs a random point on the line between each pair of input points.
    Inheriting from _Merge is a little messy but it was the quickest solution I could
    think of. Improvements appreciated."""

    def __init__(self, BATCH_SIZE):
        super(_Merge, self).__init__()
        self.BATCH_SIZE = BATCH_SIZE
        
    def _merge_function(self, inputs):
        weights = K.random_uniform((self.BATCH_SIZE, 1))
        return (weights * inputs[0]) + ((1 - weights) * inputs[1])


def generate_noise(lam, size):
    return np.random.normal(0, 1, size).astype(np.float32) + np.random.poisson(lam, size).astype(np.float32)

def generate_samples(generator_model, outfile, random_dim, lam, n_samples=10):
    """Feeds random seeds into the generator and tiles and saves the output to a PNG
    file."""
    test_sample_stack = generator_model.predict(generate_noise(lam, (n_samples, random_dim)))
    np.savetxt(outfile, test_sample_stack, delimiter=",")


# Save the generator and discriminator networks (and weights) for later use
def save_models(gan, discriminator, generator, output_dir, epoch):
    gan.save(output_dir + "/models/gan_epoch_" + str(epoch) + ".h5")
    generator.save(output_dir + "/models/generator_epoch_" + str(epoch) + ".h5")
    discriminator.save(output_dir + "/models/discriminator_epoch_" + str(epoch) + ".h5")


# Plot the loss from each batch
def plotLosses(output_dir, dLosses, gLosses, epoch):
    plt.figure(figsize=(10, 8))
    plt.plot(dLosses, label='Discriminitive loss')
    plt.plot(gLosses, label='Generative loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(output_dir + "/images/gan_loss_epoch_" + str(epoch) + ".png")

    
def train(n_epochs, bookkeeping_interval, TRAINING_RATIO, BATCH_SIZE, random_dim, lam, input_file, output_dir):
    # First we load the expression data
    X_train = load_data(input_file)
    print(X_train.shape)

    # input dimension
    input_dimension = X_train.shape[1]
    n_genes = input_dimension - 1

    # Now we initialize the generator and discriminator.
    generator = make_generator(input_dimension, random_dim)
    discriminator = make_discriminator(input_dimension)

    # The generator_model is used when we want to train the generator layers.
    # As such, we ensure that the discriminator layers are not trainable.
    # Note that once we compile this model, updating .trainable will have no effect within
    # it. As such, it won't cause problems if we later set discriminator.trainable = True
    # for the discriminator_model, as long as we compile the generator_model first.
    discriminator.trainable = False
    generator_input = Input(shape=(random_dim,))
    generator_layers = generator(generator_input)
    discriminator_layers_for_generator = discriminator(generator_layers)
    generator_model = Model(inputs=[generator_input],
                            outputs=[discriminator_layers_for_generator])
    # We use the Adam paramaters from Gulrajani et al.
    generator_model.compile(optimizer=Adam(0.0001, beta_1=0.5, beta_2=0.9),
                            loss=wasserstein_loss, metrics=["accuracy"])
    
    # Now that the generator_model is compiled, we can make the discriminator
    # layers trainable.
    discriminator.trainable = True
    generator.trainable = False
    
    # The discriminator_model is more complex. It takes both real image samples and random
    # noise seeds as input. The noise seed is run through the generator model to get
    # generated images. Both real and generated images are then run through the
    # discriminator. Although we could concatenate the real and generated images into a
    # single tensor, we don't (see model compilation for why).
    real_samples = Input(shape=X_train.shape[1:])
    generator_input_for_discriminator = Input(shape=(random_dim,))
    generated_samples_for_discriminator = generator(generator_input_for_discriminator)
    discriminator_output_from_generator = discriminator(generated_samples_for_discriminator)
    discriminator_output_from_real_samples = discriminator(real_samples)
    
    # We also need to generate weighted-averages of real and generated samples,
    # to use for the gradient norm penalty.
    averaged_samples = RandomWeightedAverage(BATCH_SIZE)([real_samples,
                                                generated_samples_for_discriminator])
    # We then run these samples through the discriminator as well. Note that we never
    # really use the discriminator output for these samples - we're only running them to
    # get the gradient norm for the gradient penalty loss.
    averaged_samples_out = discriminator(averaged_samples)
    
    # The gradient penalty loss function requires the input averaged samples to get
    # gradients. However, Keras loss functions can only have two arguments, y_true and
    # y_pred. We get around this by making a partial() of the function with the averaged
    # samples here.
    partial_gp_loss = partial(gradient_penalty_loss,
                              averaged_samples=averaged_samples,
                              gradient_penalty_weight=GRADIENT_PENALTY_WEIGHT)
    # Functions need names or Keras will throw an error
    partial_gp_loss.__name__ = 'gradient_penalty'
    
    # Keras requires that inputs and outputs have the same number of samples. This is why
    # we didn't concatenate the real samples and generated samples before passing them to
    # the discriminator: If we had, it would create an output with 2 * BATCH_SIZE samples,
    # while the output of the "averaged" samples for gradient penalty
    # would have only BATCH_SIZE samples.
    
    # If we don't concatenate the real and generated samples, however, we get three
    # outputs: One of the generated samples, one of the real samples, and one of the
    # averaged samples, all of size BATCH_SIZE. This works neatly!
    discriminator_model = Model(inputs=[real_samples,
                                        generator_input_for_discriminator],
                                outputs=[discriminator_output_from_real_samples,
                                         discriminator_output_from_generator,
                                         averaged_samples_out])
    # We use the Adam paramaters from Gulrajani et al. We use the Wasserstein loss for both
    # the real and generated samples, and the gradient penalty loss for the averaged samples
    discriminator_model.compile(optimizer=Adam(0.0001, beta_1=0.5, beta_2=0.9),
                                loss=[wasserstein_loss,
                                      wasserstein_loss,
                                      partial_gp_loss], metrics=["accuracy"])
    # We make three label vectors for training. positive_y is the label vector for real
    # samples, with value 1. negative_y is the label vector for generated samples, with
    # value -1. The dummy_y vector is passed to the gradient_penalty loss function and
    # is not used.
    positive_y = np.ones((BATCH_SIZE, 1), dtype=np.float32)
    negative_y = -positive_y
    dummy_y = np.zeros((BATCH_SIZE, 1), dtype=np.float32)

    for epoch in range(n_epochs):
        np.random.shuffle(X_train)
        print("Epoch: ", epoch)
        print("batch size:", BATCH_SIZE)
        print("training ratio:", TRAINING_RATIO)
        minibatches_size = BATCH_SIZE * TRAINING_RATIO
        n_batches = max(1, int(X_train.shape[0] // minibatches_size))
        print("Number of batches: ", n_batches)
        discriminator_loss = []
        generator_loss = []
        for i in range(n_batches):
            discriminator_minibatches = X_train[np.random.randint(0, X_train.shape[0], size=minibatches_size)]
            for j in range(TRAINING_RATIO):
                image_batch = discriminator_minibatches[j * BATCH_SIZE:
                                                        (j + 1) * BATCH_SIZE]
                noise = generate_noise(lam, (BATCH_SIZE, random_dim))
                d_loss = discriminator_model.train_on_batch(
                    [image_batch, noise],
                    [positive_y, negative_y, dummy_y])
                discriminator_loss.append(d_loss[0])
            g_loss = generator_model.train_on_batch(np.random.rand(BATCH_SIZE,
                                                        random_dim),
                                                        positive_y)
            generator_loss.append(g_loss[0])
            # Still needs some code to display losses from the generator and discriminator,
            # progress bars, etc.
        print("g_loss =", g_loss, "d_loss =", d_loss)
        if epoch % bookkeeping_interval == 0:
            generate_samples(generator, os.path.join(output_dir, "samples", 'epoch_{}.csv'.format(epoch)), random_dim, lam, X_train.shape[0])
            save_models(generator_model, discriminator_model, generator, output_dir, epoch)
            plotLosses(output_dir, discriminator_loss, generator_loss, epoch)
    generate_samples(generator, os.path.join(output_dir, "samples", 'epoch_{}.csv'.format(n_epochs)), random_dim, lam, X_train.shape[0])
    save_models(generator_model, discriminator_model, generator, output_dir, n_epochs)
    plotLosses(output_dir, discriminator_loss, generator_loss, n_epochs)

        
def generate(model_file, output_file, random_dim, lam, n_samples):
    generator_model = load_model(model_file)
    generate_samples(generator_model, output_file, random_dim, lam, n_samples)


def visualize(model_file, output_file, latent_vector):
    generator_model = load_model(model_file)
    test_sample_stack = generator_model.predict(np.array(latent_vector.split(",")).astype(np.float32).reshape(1,200))
    np.savetxt(output_file, test_sample_stack, delimiter=",")
    
    
def mkdirs(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def train2(args):
    mkdirs(args.output_dir + "/samples")
    mkdirs(args.output_dir + "/images")
    mkdirs(args.output_dir + "/models")
    train(args.n_epochs, args.bookkeeping_interval, args.training_ratio, args.batch_size, args.latent_space_dimension, args.poisson_lambda, args.input_file, args.output_dir)

                         
def generate2(args):
    generate(args.model_file, args.output_file, args.latent_space_dimension, args.poisson_lambda, args.n_samples)

    
def visualize2(args):
    visualize(args.model_file, args.output_file, args.latent_vector)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')

    subparsers = parser.add_subparsers(help='sub-command help')
    parser_train = subparsers.add_parser('train', help='train help')
    parser_train.add_argument('--input_file', required=True, help='input file')
    parser_train.add_argument('--output_dir', required=True, help='output dir')
    parser_train.add_argument('--n_epochs', required=False, type=int, default=100, help='number of epochs')
    parser_train.add_argument('--batch_size', required=False, type=int, default=64, help='batch size')
    parser_train.add_argument('--latent_space_dimension', required=False, type=int, default=200, help='latent space dimension')
    # The results are a little better when the dimensionality of the random vector is only 10.
    # The dimensionality has been left at 200 for consistency with other GAN implementations.
    parser_train.add_argument('--poisson_lambda', required=False, type=float, default=1, help='the lambda parameter of the poisson distribution for generating noise for latent space')
    parser_train.add_argument('--training_ratio', required=False, type=int, default=5, help='train ratio')
    parser_train.add_argument('--bookkeeping_interval', required=False, type=int, default=100, help='train ratio')

    parser_train.set_defaults(func=train2)

    parser_generate = subparsers.add_parser("generate", help="generate help")
    parser_generate.add_argument("--model_file", required=True, help="model file")
    parser_generate.add_argument("--n_samples", required=True, help="number of synthetics samples to be generated", type=int)
    parser_generate.add_argument('--output_file', required=True, help='output file')
    parser_generate.set_defaults(func=generate2)

    parser_generate = subparsers.add_parser("generate_with_latent_vector", help="visualize help")
    parser_generate.add_argument("--model_file", required=True, help="model file")
    parser_generate.add_argument("--latent_vector", required=True, help="csv latent vector 200 dim")
    parser_generate.add_argument('--output_file', required=True, help='output file')
    parser_generate.set_defaults(func=visualize2)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

