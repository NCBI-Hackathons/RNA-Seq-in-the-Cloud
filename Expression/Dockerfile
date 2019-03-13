FROM bioconductor/release_core2:R3.5.2_Bioc3.8

# builds 

ENV DEBIAN_FRONTEND=noninteractive 

# Install general prerequisites.  
# The Bioconductor core already includes a specific version of BioC and R

RUN apt-get -y  install \
    xsltproc \ 
    libxml2-utils \
    python3-pip 

# Install Entrez Direct: E-utilities on the UNIX Command Line
# adapted from https://www.ncbi.nlm.nih.gov/books/NBK179288/
# For use with efetch to get xml metadata
# For future expansion/use


RUN wget -O -  ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz | tar xz && \
    /edirect/setup.sh 

ENV PATH="/edirect:${PATH}"

# Install MetadataTable from previous hackathon for parsing xml
# For future expansion/use

RUN wget -O /tmp/mdt.zip https://zenodo.org/record/1158314/files/NCBI-Hackathons/MetadataTable-v1.0.0.zip?download=1 && \
    unzip /tmp/mdt.zip

RUN pip3 install -r /NCBI-Hackathons-MetadataTable-fe9bfec/requirements.txt && \
    ln -s /NCBI-Hackathons-MetadataTable-fe9bfec /MetadataTable

# Install python prerequisites for our tools

RUN pip3 install \
    gcsfs \
    matplotlib \
    numpy \
    pandas \
    seaborn \ 
    sklearn

# Install R prerequisites for our tools


RUN R -q -e 'BiocManager::install("DESeq2", version = "3.8")' && \
    R -q -e 'install.packages("optparse")' && \
    rm -rf /tmp/Rtmp*

ADD deseq2/deseq2_contrast_groups.R pooling_pipeline.py dan/counts_pca.py /
