#!/bin/bash

for ssp in 0 1 2 
do
	echo $ssp
		for ens in 0 1 2 3 4
		do
			echo $ens
				for rho in 790 850 910
				do
					echo $rho	
					python calibration_isimip.py $ens $ssp 40 $rho 
				done
		done
done
