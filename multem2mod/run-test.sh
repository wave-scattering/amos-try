#!/bin/bash
rm -rf test
mkdir test
cd test
git clone https://github.com/wave-scattering/amos-try.git
cd amos-try
cmake .
make
cd multem2mod
./circle-BIC-spectra-plot.py
./go.sh
