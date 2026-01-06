FROM ubuntu:26.04 AS builder

# --- Development / Build Stage ---
# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    build-essential \
    gengetopt \
    libbz2-dev \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# A .dockerignore file is recommended to exclude unnecessary files
COPY . .

RUN autoreconf -i \
    && ./configure CFLAGS="-O3 -ffast-math" \
    && make \
    && make install

# --- Production Stage ---
FROM ubuntu:26.04 as production

# Install only runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-1.0 \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

# Copy the built binaries from the builder stage
COPY --from=builder /usr/local/bin/debyer /usr/local/bin/
COPY --from=builder /usr/local/bin/dbr_* /usr/local/bin/

CMD ["debyer", "--help"]