#!/usr/bin/env bash

if [[ $# -ne 2 ]]; then
	echo "Please, provide 2 file as input: "
	echo "./execute.sh <input_file> <output_file>"
	echo "Example: ./execute.sh param.dat out_param.dat"
	exit 1
fi

echo "Compiling the code..."
./compile_2021.sh
echo "Executing the code..."
./run_2021.sh "$1" > "$2"
echo "Plotting the results..."
root -l 'plot.C("'"$2"'")'