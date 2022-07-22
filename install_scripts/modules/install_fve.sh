#!/bin/bash

CONDA_SOURCE=$1
FVE_ENV=$2
FVE_DIR=$3

source $CONDA_SOURCE
conda activate FVE_ENV

cd $FVE_DIR

javac -d bin src/*.java

conda deactivate