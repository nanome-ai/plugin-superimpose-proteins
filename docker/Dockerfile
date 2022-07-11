FROM continuumio/miniconda3:4.11.0 as production
ENV ARGS=''

WORKDIR /app
# RUN apt update && apt upgrade -y && apt install -y build-essential mpich
# Install and compile FPocket
# ARG FPOCKET_VERSION=4.0.1
# RUN git clone -b ${FPOCKET_VERSION} https://github.com/Discngine/fpocket.git
# RUN cd fpocket && make && make install

# Install Site-Motif
# COPY site-motif /app/site-motif
# RUN chmod +x /app/site-motif/bin/site-motif
# ENV PATH="/app/site-motif/bin/:${PATH}"

COPY plugin /app/plugin

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .
CMD python run.py ${ARGS}
