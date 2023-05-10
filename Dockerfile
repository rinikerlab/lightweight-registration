FROM continuumio/miniconda3

ENV PYTHON_VERSION=3.11
ENV RDKIT_VERSION=2023.03.1

# RUN apt update && apt install libtiff5 -y
RUN conda install -c conda-forge python=$PYTHON_VERSION
RUN conda install -c conda-forge rdkit=2023.03.1

COPY . /lw-reg
WORKDIR /lw-reg

RUN pip install --editable .