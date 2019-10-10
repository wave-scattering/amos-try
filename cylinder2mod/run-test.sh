#!/bin/bash
rm -rf test
mkdir test
cd test
git clone https://github.com/wave-scattering/amos-try.git
cd amos-try
cmake .
make
cd cylinder2mod
./cylinder2mod && ./plot-ext.py
