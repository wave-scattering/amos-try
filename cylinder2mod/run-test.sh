#!/bin/bash
ROOT_DIR=`pwd`
echo $ROOT_DIR
rm -rf test
mkdir test
cd test
git clone https://github.com/wave-scattering/amos-try.git
cd amos-try
cmake .
make
cd cylinder2mod
rm -rf *.dat
./cylinder2mod
./plot-ext.py
cd $ROOT_DIR
rm -rf test
