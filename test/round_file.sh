#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 file1"
  exit 1
fi

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
        sub(/\.000$/, "", $i)  # Optional: strip ".000" from floats like "5.000"
      }
    }
    print
  }
  ' "$1"
}

round_file "$1"
