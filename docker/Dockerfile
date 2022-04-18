FROM python:3
ENV ARGS=''

WORKDIR /app

# Install and compile FPocket
ARG FPOCKET_VERSION=4.0.1
RUN git clone -b ${FPOCKET_VERSION} https://github.com/Discngine/fpocket.git
RUN cd fpocket && make && make install

# TODO: install non-pip dependencies here
ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .



CMD python run.py ${ARGS}
