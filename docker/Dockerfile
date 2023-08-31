FROM ubuntu:20.04
LABEL maintainer = "Jose Gabriel Nino Barreat <jose.ninobarreat@biology.ox.ac.uk>"

# Set paths to the minoconda bin directory
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# Get updates before installing software
RUN apt-get update

# Install wget and remove /var/lib/apt/lists/* (ubuntu package lists, not needed)
RUN apt-get -y install wget && rm -rf /var/lib/apt/lists/*

# Download and install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh -O ~/miniconda.sh \
    && mkdir /root/.conda \
    && bash ~/miniconda.sh -b \
    && rm -f ~/miniconda.sh

# Copy environment.yml into image (must be in build context directory)
COPY environment.yml .

# Update conda and create the hi-fever environment
RUN conda install conda=23.7.3 \
    && conda env create -f environment.yml \
    && conda clean --all --yes

# Load environment inside image and set path variables
RUN echo "source activate hi-fever" > ~/.bashrc
ENV PATH /root/miniconda3/envs/hi-fever/bin:$PATH
ENV CONDA_PREFIX=/root/miniconda3/envs/hi-fever
