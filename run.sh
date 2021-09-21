#!/usr/bin/env bash

param="$1"

if [[ -f "$param" ]] && [[ $# == 1  ]]; then
	while IFS= read -r line;
	do
		if [[ "$line" == "#"* ]]; then
			true
		else 
			#echo "$line"
			IFS=" " read -r -a arr <<< "$line"
			if [[ ${#arr[@]} -ne 12 ]]; then
				echo "Please provide 1 text file containing 12 arguments: " 
				echo "pbar lab momentum, target mass, target charge, nuclear strength (real), nuclear strength (img), nuclear radius (real), nuclear radius (img), coulomb radius, diffusness (real), diffusness (img)"
				echo "Example: ./run.sh param.dat"
				exit 1
			fi
			#echo "${arr[@]:0:10}"
			#echo "${arr[@]:0:3}"
			./bin/antip "${arr[@]:0:12}" > /tmp/input.dat
			# echo "${arr[@]:0:3}" "${arr[11]}"
			./bin/scattering /tmp/input.dat "${arr[@]:0:3}" "${arr[11]}"
		fi
	done < "$param"
else
	echo "Please provide 1 text file: "\n 
	echo "Example: ./run.sh param.dat"
	exit 1
fi