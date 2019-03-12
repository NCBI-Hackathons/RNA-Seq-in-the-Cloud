FROM bioconductor/release_core2:R3.5.2_Bioc3.8
ENV DEBIAN_FRONTEND=noninteractive 
# RUN apt-get -y update ; 

RUN apt-get -y  install \
    libxml2-utils \
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

RUN wget -O -  ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz | tar xz && \
    /edirect/setup.sh 

ENV PATH="/edirect:${PATH}"
