#!/bin/bash
$VFDIR/bin/SneddonCrack2D.py --lz .01 -o sneddon-gc1.dat
$VFDIR/bin/SneddonCrack2D.py --lz .01 --gc 1.25 -o sneddon-gc1.25.dat

$VFDIR/bin/PressurePlot2.py 604397 605292 605062 604398 606674 606675 606657 606658 sneddon-gc1.25.dat -o test15-N251.pdf
$VFDIR/bin/PressurePlot2.py 605242 605243 605244 605245 605246 605247 605248 sneddon-gc1.25.dat -o test15-N355.pdf
$VFDIR/bin/PressurePlot2.py 605254 605255 605256 605257 605258 605259 605260 sneddon-gc1.25.dat -o test15-N501.pdf
$VFDIR/bin/PressurePlot2.py 605571 605572 605573 605574 605575 605576 605577 sneddon-gc1.25.dat -o test15-N707.pdf

$VFDIR/bin/PressurePlot2.py 605573 605256 605244 604398 sneddon-gc1.25.dat -o test15-ratio1.0.pdf
$VFDIR/bin/PressurePlot2.py 605574 605257 605245 606674 sneddon-gc1.25.dat -o test15-ratio1.5.pdf
$VFDIR/bin/PressurePlot2.py 605575 605258 605246 606675 sneddon-gc1.25.dat -o test15-ratio2.0.pdf
$VFDIR/bin/PressurePlot2.py 605576 605259 605247 606657 sneddon-gc1.25.dat -o test15-ratio2.5.pdf


$VFDIR/bin/CrackLengthPlot.py 604397 605292 605062 604398 606674 606675 606657 606658  sneddon-gc1.25.dat -o test15-length-N251.pdf
$VFDIR/bin/CrackLengthPlot.py 605242 605243 605244 605245 605246 605247 605248 sneddon-gc1.25.dat -o test15-length-N355.pdf
$VFDIR/bin/CrackLengthPlot.py 605254 605255 605256 605257 605258 605259 605260 sneddon-gc1.25.dat -o test15-length-N501.pdf
$VFDIR/bin/CrackLengthPlot.py 605571 605572 605573 605574 605575 605576 605577 sneddon-gc1.25.dat -o test15-length-N707.pdf

$VFDIR/bin/CrackLengthPlot.py 605573 605256 605244 604398 sneddon-gc1.25.dat -o test15-length-ratio1.0.pdf
$VFDIR/bin/CrackLengthPlot.py 605574 605257 605245 606674 sneddon-gc1.25.dat -o test15-length-ratio1.5.pdf
$VFDIR/bin/CrackLengthPlot.py 605575 605258 605246 606675 sneddon-gc1.25.dat -o test15-length-ratio2.0.pdf
$VFDIR/bin/CrackLengthPlot.py 605576 605259 605247 606657 sneddon-gc1.25.dat -o test15-length-ratio2.5.pdf

$VFDIR/bin/PressurePlot2.py 608976 608977 608978 sneddon-gc1.dat -o test15-newscaling-N251.pdf
$VFDIR/bin/PressurePlot2.py 609286 609287 609288 sneddon-gc1.dat -o test15-newscaling-N355.pdf


$VFDIR/bin/CrackLengthPlot.py 608976 608977 608978 sneddon-gc1.dat -o test15-newscaling-N251-length.pdf
$VFDIR/bin/CrackLengthPlot.py 609286 609287 609288 sneddon-gc1.dat -o test15-newscaling-N355-length.pdf

