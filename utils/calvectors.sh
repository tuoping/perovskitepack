#!/bin/bash
for d in annealT{52..63} annealT{65..105} annealT{110..135..5} annealT{136..145} annealT{147..150} annealT{152..154};
do
    echo $d;
    cd $d;
    # tail -6153 md4.dump > 4020000.lammpstrj;
    python ../perovskitecall.py *.lammpstrj;
    # python ../drawhist.py molecule_longaxis.dat
    # python ../drawhist.py molecule_polaraxis.dat
    cd ../;
done
