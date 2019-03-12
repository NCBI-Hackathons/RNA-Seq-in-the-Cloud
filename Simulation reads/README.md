
python3

```
pip install numpy pandas sklearn keras tensorflow-gpu matplotlib pillow argparse
```

## training
```
python src/exp_gan.py train --input_file data/test.csv --output_dir=train --n_epoch=500 --batch_size=10 --training_ratio=1
```

## generating samples
```
python src/exp_gan.py generate --model_file train/models/generator_epoch_499.h5 --output_file output/test.csv --n_samples 10
```
