#!/bin/bash
$VFDIR/bin/SneddonCrack2D.py --lz .01 -o sneddon-gc1.dat
$VFDIR/bin/PressurePlot2.py   611699 611700 611701 611967 611968 612009 sneddon-gc1.dat -o sneddon-N355-pressure.pdf
$VFDIR/bin/CrackLengthPlot.py 611699 611700 611701 611967 611968 612009 sneddon-gc1.dat -o sneddon-N355-length.pdf

$VFDIR/bin/PressurePlot2.py   611701 611766 611771 sneddon-gc1.dat --title 'Pressure v.s. injected volume, case 1' -o sneddon-scaling-pressure1.pdf
$VFDIR/bin/CrackLengthPlot.py   611701 611766 611771 sneddon-gc1.dat --title 'Crack increment v.s. injected volume, case 1' -o sneddon-scaling-length1.pdf
$VFDIR/bin/PressurePlot2.py 611739 611767 611776  sneddon-gc1.dat --title 'Pressure v.s. injected volume, case 6' -o sneddon-scaling-pressure2.pdf
$VFDIR/bin/CrackLengthPlot.py 611739 611767 611776  sneddon-gc1.dat --title 'Crack increment v.s. injected volume, case 6'  -o sneddon-scaling-length2.pdf
