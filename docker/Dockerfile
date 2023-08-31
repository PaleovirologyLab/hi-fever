FROM ubuntu:20.04
MAINTAINER Jose Gabriel Nino Barreat <jose.ninobarreat@biology.ox.ac.uk>

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update

RUN apt-get -y install wget && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh -O ~/miniconda.sh \
    && mkdir /root/.conda \
    && bash ~/miniconda.sh -b \
    && rm -f ~/miniconda.sh

COPY environment.yml .

RUN conda install conda=23.7.3 \
    && conda env create -f environment.yml \
    && conda clean --all --yes
RUN echo "source activate hi-fever" > ~/.bashrc
ENV PATH /root/miniconda3/envs/hi-fever/bin:$PATH
ENV CONDA_PREFIX=/root/miniconda3/envs/hi-fever
