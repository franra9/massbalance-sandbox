#!/bin/bash

#ssps=("0" "1" "2")
#enss=("0" "1" "2" "3" "4")
for ssp in 0 1 2
do
	echo $ssp
		for ens in 0 1 2 3 4
		do
			echo $ens	
			python calibration_isimip.py $ens $ssp 10 
		done
done
