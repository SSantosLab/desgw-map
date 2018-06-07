#!/bin/sh

# for use without installing desdw in the official pythonic way
# I prefer this method for research code, anyway
#

#export DESGW_DIR="/data/des30.a/data/annis/test/"
unset PYTHONPATH
export DESGW_DIR="/data/des60.a/data/alenon/desgw-map"
export DESGW_DATA_DIR="$DESGW_DIR/data/"
export PYTHONPATH="$DESGW_DIR/python/:$PYTHONPATH"

