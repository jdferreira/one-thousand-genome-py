#!/usr/bin/env bash

HREFS="input/1kgp_hrefs.txt"
DATADIR="data"

# Prepare the directory structure
prepare() {
    mkdir -p "$DATADIR"
}

# Read the contents of the HREFS file and download the files
download() {
    awk 'length($0) > 0 && $0 !~ /^#/' "$HREFS" | while read line; do
        echo "$line"
        (cd "$DATADIR"; curl -O "$line")
    done;
}

# Run the script
prepare
download
