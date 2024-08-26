FROM continuumio/miniconda3

ENV PYTHON_VERSION=3.11
ENV RDKIT_VERSION=2023.03.1

RUN conda install -c conda-forge python>=$PYTHON_VERSION
RUN conda install -c conda-forge rdkit>=RDKIT_VERSION

COPY . /lw-reg
WORKDIR /lw-reg

RUN pip install .