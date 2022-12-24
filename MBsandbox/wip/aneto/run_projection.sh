#!/bin/bash

#ssps=("0" "1" "2")
#enss=("0" "1" "2" "3" "4")
for ssp in 0 1 2
do
	echo $ssp
		for ens in 0 1 2 3 4
		do
			echo $ens
				for rho in 850
				do
					echo $rho	
					python projection.py $ens $ssp 40 $rho 
				done
		done
done
