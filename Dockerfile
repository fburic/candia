FROM continuumio/miniconda3:4.8.2

ADD paradias_env.yaml /tmp/paradias_env.yaml

# Update conda to at least 4.8.3 https://github.com/conda/conda/issues/9681
# then install packages from env file
RUN apt-get install -y libgomp1 \
    && conda update conda \
    && conda env update --name base --file /tmp/paradias_env.yaml  --prune \
    && conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && echo "source activate" > ~/.bashrc

ENV PATH /opt/conda/bin:$PATH

WORKDIR /app
COPY scripts /app/scripts
COPY test /app/test

CMD ["/bin/bash", "--login", "-c"]

