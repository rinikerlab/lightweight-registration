FROM mambaorg/micromamba:latest


# prevent buffering of output which can cause missed logs
ENV PYTHONUNBUFFERED=1 

COPY --chown=$MAMBA_USER:$MAMBA_USER . /home/mambauser/lw-reg

WORKDIR /home/mambauser/lw-reg

COPY --chown=$MAMBA_USER:$MAMBA_USER  environment.yml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN pip install .
RUN lwreg greet
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"  # Not a preferred method!
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]