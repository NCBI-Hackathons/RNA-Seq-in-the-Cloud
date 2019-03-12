FROM bioconductor/release_core2:R3.5.2_Bioc3.8
ENV DEBIAN_FRONTEND=noninteractive 
# RUN apt-get -y update ; 

RUN apt-get -y  install \
    python3-pip 

RUN pip3 install \
    gcsfs \
    matplotlib \
    numpy \
    pandas \
    seaborn \ 
    sklearn


RUN R -q -e 'BiocManager::install("DESeq2", version = "3.8")' && \
    rm -rf /tmp/Rtmp*

# import sys
# import sklearn
# import pandas as pd
# import gcsfs
# import seaborn as sns
# import numpy as np
# import matplotlib.pyplot as plt
# import math
# import random
# from sklearn.preprocessing import StandardScaler
# from sklearn.decomposition import PCA
# from sklearn.pipeline import Pipeline
