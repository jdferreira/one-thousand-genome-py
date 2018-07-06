#!/usr/bin/env bash

DATADIR="data"

# Get the command line arguments
INPUT=$1; shift

# Extract the headers of the VCF files
extract_headers() {
    # Read the vcf files (gzipped)
    for fname in "$DATADIR"/*.vcf.gz; do
        # Where to output the files
        ofname="${fname##$DATADIR/}"
        ofname="${ofname%%.vcf.gz}"
        zcat "$fname" \
            | awk '/^#/ { print } /^[^#]/ { exit }' \
            > "$DATADIR/$ofname.header.txt"
    done
}

# Make sure that the headers of the chromosomes are sane
validate_headers() {
    valid="T"
    pivot=$(ls "$DATADIR"/*.header.txt | head -n 1)
    for fname in "$DATADIR"/*.header.txt; do
        diff "$pivot" "$fname" > /dev/null
        if [ "$?" = 1 ]; then
            echo "Header for $fname differs" >&2
            valid="F"
        fi
    done
    
    if [ "$valid" = "F" ]; then
        echo -n "Some headers differ. Continue? (Y/N) " >&2
        read answer
        if [ "$answer" != "Y" ]; then
            exit
        fi
    fi
}

# Print the (possibly unique) header
print_headers() {
    cat "$(ls "$DATADIR"/*.header.txt | head -n 1)"
    rm "$DATADIR"/*.header.txt
}

# Extract the relevant polymorphism data from the VCF files
extract_relevant() {
    # Make a regular expression pattern that matches the relevant polymorphism
    # identifiers
    PATTERN="\b(?:"$(awk 'length($0) > 0 && $0 !~ /^#/' "$INPUT" | paste -s -d'|')")\b"
    
    for fname in "$DATADIR"/*.vcf.gz; do
        # Where to output the files
        ofname="${fname##$DATADIR/}"
        ofname="${ofname%%.vcf.gz}"
        
        # Pipe through pv to get a progress bar
        name="${fname##$DATADIR/}"
        name="${name%%vcf.gz}"
        if [ "${#name}" -gt 14 ]; then
            name="$(echo "$name" | cut -c 1-10) ..."
        elif [ "${#name}" -lt 14 ]; then
            name="$(printf "%-14s" "$name")"
        fi
        
        pv -N "[$name]" "$fname" | zcat | grep -P "$PATTERN"
    done
}

# Run the script
extract_headers
validate_headers
print_headers
extract_relevant
