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

# TODO: restrict what we copy into the image
WORKDIR /app
COPY . .

RUN autoreconf -i \
    && ./configure CFLAGS="-O3 -ffast-math" \
    && make \
    && make install

CMD [ "debyer", "--help"]