#!/bin/bash
make
#cp input_sadrieva_band.txt fort.10 && ./multem2
cp input2_band.txt fort.10 && ./multem2
#cp input_sadrieva.txt fort.10 && ./multem2

#cp input1.txt fort.10 && ./multem2
#cp input_inoue.txt fort.10 && ./multem2
#cp input_andrey.txt fort.10 && ./multem2
#cp input_zarina.txt fort.10 && ./multem2

./plot_band.py
