#!/bin/bash

# Usage: ./script.sh [optional_input_file]
# If no file is provided, read from stdin.

input="${1:-/dev/stdin}"

round_file() {
  awk '
  BEGIN { OFS = "\t" }
  {
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^[+-]?[0-9]+$/) {
        # Integer — keep as-is
        continue
      } else if ($i ~ /^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/) {
        # Float — round to 3 decimals
        $i = sprintf("%.3f", $i)
        sub(/\.000$/, "", $i)  # Strip ".000" from whole numbers
      }
    }
    print
  }
  ' "$1"
}

round_file "$input"


