FROM python:3.8
ENV ARGS=''

WORKDIR /app
RUN apt update && apt upgrade -y && apt install -y build-essential mpich

# Install and compile FPocket
ARG FPOCKET_VERSION=4.0.1
RUN git clone -b ${FPOCKET_VERSION} https://github.com/Discngine/fpocket.git
RUN cd fpocket && make && make install

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY plugin /app/plugin

COPY . .
CMD python run.py ${ARGS}
