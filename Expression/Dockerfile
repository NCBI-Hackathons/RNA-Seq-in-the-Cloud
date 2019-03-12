FROM bioconductor/release_core2:R3.5.2_Bioc3.8
ENV DEBIAN_FRONTEND=noninteractive 
# RUN apt-get -y update ; 

RUN apt-get -y  install \
    libxml2-utils \
    python3-pip 

RUN wget -O -  ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz | tar xz && \
    /edirect/setup.sh 

ENV PATH="/edirect:${PATH}"

RUN wget -O /tmp/mdt.zip https://zenodo.org/record/1158314/files/NCBI-Hackathons/MetadataTable-v1.0.0.zip?download=1 && \
    unzip /tmp/mdt.zip

RUN pip3 install -r /NCBI-Hackathons-MetadataTable-fe9bfec/requirements.txt && \
    ln -s /NCBI-Hackathons-MetadataTable-fe9bfec /MetadataTable

RUN pip3 install \
    gcsfs \
    matplotlib \
    numpy \
    pandas \
    seaborn \ 
    sklearn


RUN R -q -e 'BiocManager::install("DESeq2", version = "3.8")' && \
    R -q -e 'install.packages("optparse")' && \
    rm -rf /tmp/Rtmp*

ADD deseq2/deseq2_contrast_groups.R /
