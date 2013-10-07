#! /bin/bash

./convolve-mtss-rotatmers.py 	--resid 47 330\
 									--pdb peptso.gro\
 									--histogramBins 0 80 1\
 									--clashDistance 2.2\
 									--output "dat/peptso-xrd"\
									--dcdfilename "dcd/peptso-xrd-47-330"	
