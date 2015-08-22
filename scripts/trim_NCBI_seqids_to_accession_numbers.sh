#!/usr/bin/env sh

if [[ $# -ne 1 ]]; then
  echo "Trims NCBI seqids in the given FASTA down to accession numbers."
  echo "Usage: $0 <fasta file>"
  exit
fi

cat $1 | ruby -pe 'gsub(/^>gi\|\d+\|\w+\|(.+)\| /) {">#{$1} "}'
