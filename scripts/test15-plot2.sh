#!/bin/bash
$VFDIR/bin/SneddonCrack2D.py --lz .01 -o sneddon-gc1.dat
$VFDIR/bin/PressurePlot2.py   611699 611700 611701 611967 611968 612009 sneddon-gc1.dat -o sneddon-N355-pressure.pdf
$VFDIR/bin/CrackLengthPlot.py 611699 611700 611701 611967 611968 612009 sneddon-gc1.dat -o sneddon-N355-length.pdf

$VFDIR/bin/PressurePlot2.py   611701 611766 611771 sneddon-gc1.dat -o sneddon-scaling-pressure.pdf
$VFDIR/bin/CrackLengthPlot.py 611739 611767 611776  sneddon-gc1.dat -o sneddon-scaling-length.pdf
