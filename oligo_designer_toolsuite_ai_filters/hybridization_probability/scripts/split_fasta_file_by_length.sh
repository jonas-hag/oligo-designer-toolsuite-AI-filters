#!/bin/env bash
# split_fasta_by_length_simple.sh
# Split a FASTA file into two files by sequence length (single-line sequences)
# Usage: ./split_fasta_by_length_simple.sh input.fasta short_out.fasta long_out.fasta [threshold]
# Default threshold is 30

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 INPUT.fasta SHORT_OUT.fasta LONG_OUT.fasta [THRESHOLD]" >&2
  exit 1
fi

infile="$1"
short_out="$2"
long_out="$3"
threshold="${4:-30}"

# Empty outputs first
: > "$short_out"
: > "$long_out"

# Process file
awk -v thr="$threshold" -v outS="$short_out" -v outL="$long_out" '
{
  if ($0 ~ /^>/) {
    header = $0
    getline seq
    if (length(seq) < thr) {
      print header >> outS
      print seq >> outS
    } else {
      print header >> outL
      print seq >> outL
    }
  }
}' "$infile"