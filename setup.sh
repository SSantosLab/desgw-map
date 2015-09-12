#!/bin/sh

# for use without installing desdw in the official pythonic way
# I prefer this method for research code, anyway
#

#export DESGW_DIR="/data/des30.a/data/annis/test/"
export DESGW_DIR="/home/s1/annis/daedalean/desgw-map"
export DESGW_DATA_DIR="$DESGW_DIR/data/"
PYTHONPATH="$PYTHONPATH:$DESGW_DIR/python/"

