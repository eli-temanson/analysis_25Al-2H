#!/bin/bash

pandoc "$1" \
    -V geometry:papersize="{25cm,15cm}" \
    -V geometry:margin="2cm" \
    -V fontsize="12pt" \
    --output="$2"
