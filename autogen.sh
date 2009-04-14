#!/bin/sh

# keep it simple
autoreconf -iv
./configure "$@"
