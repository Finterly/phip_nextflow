FROM rocker/r-ubuntu:22.04

# Install  utilities
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    python3 python3-dev python3-pip python3-venv git git-lfs default-jdk ant \
    libbz2-dev libsdl1.2-dev liblzma-dev libcurl4-openssl-dev zlib1g-dev libxml2-dev \
	r-cran-tidyverse bwa samtools multiqc datamash && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Install conda packages
RUN conda install -c conda-forge -c bioconda cutadapt
RUN conda install -c conda-forge -c bioconda fastqc

# Set up R
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'

# Update apt-get
RUN Rscript -e 'install.packages("remotes", version = "2.4.2")'
RUN Rscript -e 'remotes::install_cran("rmarkdown",upgrade="never", version = "2.19")'
