FROM continuumio/miniconda3:4.11.0 as production
ENV ARGS=''

WORKDIR /app
RUN apt update && apt upgrade -y && apt install -y build-essential mpich
# Install and compile FPocket
ARG FPOCKET_VERSION=4.0.1
RUN git clone -b ${FPOCKET_VERSION} https://github.com/Discngine/fpocket.git
RUN cd fpocket && make && make install

COPY rmsd_v2 /app/rmsd_v2

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

CMD python run.py ${ARGS}

##############################
FROM production as cookbook

# RUN pip install jupyterlab pip
COPY cookbook/cookbook_environ.yaml .
RUN conda env update --file cookbook_environ.yaml

COPY cookbook /app/cookbook
RUN ln -s /app/rmsd_v2 /app/cookbook/rmsd_v2
RUN ln -s /app/rmsd_v2 /opt/conda/lib/python3.9/site-packages/rmsd_v2

RUN ln -s /app/cookbook /opt/conda/lib/python3.9/site-packages/cookbook

WORKDIR /app/cookbook/notebooks
RUN jupyter notebook --generate-config
CMD jupyter lab --ip 0.0.0.0 --allow-root
