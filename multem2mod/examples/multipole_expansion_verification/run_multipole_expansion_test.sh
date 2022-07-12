#!/bin/bash
cd ../../..
cmake .
make
cd multem2mod/examples/multipole_expansion_verification
cp ../../multem2 .
python3 DOI_10_1103_PhysRevB_95_195406_comparison.py

