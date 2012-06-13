#!/bin/bash
$VFDIR/bin/SneddonCrack2D.py --lz .01 -o sneddon-gc1.dat
$VFDIR/bin/PressurePlot2.py 610480 610481 610482 609286 609287 sneddon-gc1.dat -o test15-newscaling-free-N355.pdf
$VFDIR/bin/CrackLengthPlot.py 610480 610481 610482 609286 609287 sneddon-gc1.dat -o test15-newscaling-free-N355-length.pdf

