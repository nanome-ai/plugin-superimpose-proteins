FROM continuumio/miniconda3:4.11.0 as production
ENV ARGS=''

WORKDIR /app
RUN apt update && apt upgrade -y && apt install -y build-essential mpich
# Install and compile FPocket
ARG FPOCKET_VERSION=4.0.1
RUN git clone -b ${FPOCKET_VERSION} https://github.com/Discngine/fpocket.git
RUN cd fpocket && make && make install

# Install Python2 environment and SiteMotif
RUN git clone https://github.com/santhoshgits/MAPP-3D /app/site_motif
RUN conda create --name py27 python=2.7
RUN conda install -n py27 numpy mpi4py

COPY rmsd_v2 /app/rmsd_v2

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt


CMD python run.py ${ARGS}

##############################
FROM production as cookbook

# Install Node and NPM.
# RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.0/install.sh | bash
# ENV NVM_DIR=/root/.nvm
# ENV NODE_VERSION=16.13.1
# RUN . "$NVM_DIR/nvm.sh" && nvm install ${NODE_VERSION}
# RUN . "$NVM_DIR/nvm.sh" && nvm use v${NODE_VERSION}
# RUN . "$NVM_DIR/nvm.sh" && nvm alias default v${NODE_VERSION}
# ENV PATH="/root/.nvm/versions/node/v${NODE_VERSION}/bin/:${PATH}"
# COPY package.json .
# RUN npm install

# RUN pip install jupyterlab pip
COPY cookbook/cookbook_environ.yaml .
RUN conda env update --file cookbook_environ.yaml

COPY cookbook /app/cookbook
RUN ln -s /app/rmsd_v2 /app/cookbook/rmsd_v2
RUN ln -s /app/rmsd_v2 /opt/conda/lib/python3.9/site-packages/rmsd_v2

RUN ln -s /app/cookbook /opt/conda/lib/python3.9/site-packages/cookbook

RUN ln -s /app/site_motif /opt/conda/lib/python3.9/site-packages/site_motif
RUN ln -s /app/site_motif /app/cookbook/site_motif

WORKDIR /app/cookbook/notebooks
RUN jupyter notebook --generate-config
CMD jupyter lab --ip 0.0.0.0 --allow-root
