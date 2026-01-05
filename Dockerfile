FROM ubuntu:26.04

RUN apt update \
    && apt install -y \
        autoconf \
        build-essential \
        gengetopt \
        libbz2-dev \
        libz-dev \
    && apt clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . . 

RUN autoreconf -i \
    && ./configure CFLAGS="-O3 -ffast-math" \
    && make \
    && make install

CMD [ "debyer", "--help"]