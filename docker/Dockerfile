FROM python:3 as production
ENV ARGS=''

WORKDIR /app

# Install and compile FPocket
ARG FPOCKET_VERSION=4.0.1
RUN git clone -b ${FPOCKET_VERSION} https://github.com/Discngine/fpocket.git
RUN cd fpocket && make && make install

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY rmsd_v2 /app/rmsd_v2

CMD python run.py ${ARGS}


##############################


FROM production as cookbook
RUN apt-get --allow-releaseinfo-change update && apt-get -y upgrade
RUN apt-get install -y curl

# Install Node and NPM.
RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.0/install.sh | bash
ENV NVM_DIR=/root/.nvm
ENV NODE_VERSION=16.13.1
RUN . "$NVM_DIR/nvm.sh" && nvm install ${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm use v${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm alias default v${NODE_VERSION}
ENV PATH="/root/.nvm/versions/node/v${NODE_VERSION}/bin/:${PATH}"
COPY package.json .
RUN npm install

RUN pip install jupyterlab pip
COPY cookbook/requirements-cookbook.txt .
RUN pip install -r requirements-cookbook.txt

COPY cookbook /app/cookbook
RUN ln -s /app/rmsd_v2 /app/cookbook/rmsd_v2
RUN ln -s /app/cookbook /usr/local/lib/python3.10/site-packages/cookbook
RUN ln -s /app/rmsd_v2 /usr/local/lib/python3.10/site-packages/rmsd_v2

WORKDIR /app/cookbook/notebooks
RUN jupyter notebook --generate-config
CMD jupyter lab --ip 0.0.0.0 --allow-root
