#! /usr/bin/env bash

# SCRIPT TO SCAN ON MOMENTUM AND ANGLE IN ANTIP_SCAN.FOR

# Function to plot in gnuplot from the output file
# of antip_scan.for (option "mom" only)
function plotGraph()
{
	gnuplot <<- EOF
		set xlabel "plab (MeV/c)"
		set ylabel "{/Symbol s} (mb)"
		set title "Integral cross section (elastic and reaction)"
		set grid
		set term png
		set output "fig/${filename}.png"
		plot "${outfile}" using 1:2 title "elastic" with linespoints, \
		"${outfile}" using 1:3 title "reaction" with linespoints, \
		"${outfile}" using 1:4 title "elastic+reaction" with linespoints
	EOF
}

# input of script: parameter file (path)
inputfile=$1
# output filename and path
filename=int_cs
outfile=./results/${filename}.dat 

# Warning message/help option
if [ $# -ne 1 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	echo "Please, provide a parameter file: "
	echo "	Usage:	./p_scan.sh <parameter_file>"
	echo "		./p_scan.sh -h (or --help)"
	echo
	echo "Check if the file is present."
	echo "In case, create a text file and provides 12 parameters:"
	echo "plab(MeV/c), A target (amu), Z target, U0 (MeV), W0(MeV), r0r(fm), r0i(fm), r0c(fm), a0r(fm), a0i(fm), option (\"mom\" or \"ang\"), angle (deg)"
	echo "For comments, use # at the beginning of line."
	echo "If you use \"mom\" option, as convention put as angle 999."
	exit 1
fi

# Compile the antip_scan.for code and create the
# executable in bin/ directory
gfortran -std=gnu src/antip_scan.for -o bin/antip_scan

while IFS= read -r line;
	do
		if [[ "$line" == "#"* ]]; then 	#line with # are not read
			true
		else 
			# echo "$line"
			# Save the line in an array
			IFS=" " read -r -a arr <<< "$line"
		fi
	done < "$1"

if [ "${arr[10]}" = "mom" ]; then
	# Remove - if present - the old output
	rm $outfile 2>/dev/null

	# Number of step
	N=10
	# Value of the step (MeV/c)
	step=10

	echo "Executing the scan on momentum..."
	for ((i=0; i<"$N"+1; i++)); do
		# Create the array for momentum, starting
		# from the value given in the text file
		arrmom[i]=$(expr $((arr+i*step))) 
	done
	for mom in "${arrmom[@]:0:$N+1}"; do
		# Calculate the elastic and inelastic
		# cross sections
		./bin/antip_scan $mom.0 "${arr[@]:1:12}" 
		echo "plab momentum "$mom
	done
	echo "Text file with results saved in results/ directory."
	echo "Creating the plot..."
	# plot with gnuplot of the cross sections
	plotGraph
	echo "Plot saved in fig/ directory."
elif [ "${arr[10]}" = "ang" ]; then
	# Calculate the differential cross section
	# at specific angle and momentum
	./bin/antip_scan "${arr[@]:0:12}"
fi
