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


# option for input
opt=$1
# output filename and path
filename=int_cs
outfile=./results/${filename}.dat 

# Warning message/help option
if [ $# -ne 2 ] && [ $# -ne 13 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	echo "Please, provide an option and parameters or parameter file: "
	echo "	Usage:	./p_scan.sh <option> <parameters>"
	echo "		./p_scan.sh -h (or --help)"
	echo
	echo "If you use the file option, check if the file is present."
	exit 1
elif [ "$1" = "-f" ] && [ "$2" = "-h" ] || [ "$2" = "--help" ]; then
	echo "Usage of option -f (file mode):"
	echo "./p_scan.sh -f <filename>"
	echo
	echo "Check if the input file is present."
	echo "In case, create a text file and provides 12 parameters:"
	echo "plab(MeV/c), A target (amu), Z target, U0 (MeV), W0(MeV), r0r(fm), r0i(fm), r0c(fm), a0r(fm), a0i(fm), option (\"mom\" or \"ang\"), angle (deg)"
	echo "For comments, use # at the beginning of line."
	echo "If you use \"mom\" option, as convention put as angle 999."
	exit 1
elif [ "$1" = "-i" ] && [ "$2" = "-h" ] || [ "$2" = "--help" ]; then
	echo "Usage of option -i (interactive mode):"
	echo "./p_scan.sh -i <p1> <p2> ... <p12>"
	echo
	echo "The user must provide 12 parameters:"
	echo "plab(MeV/c), A target (amu), Z target, U0 (MeV), W0(MeV), r0r(fm), r0i(fm), r0c(fm), a0r(fm), a0i(fm), option (\"mom\" or \"ang\"), angle (deg)"
	echo "If you use \"mom\" option, as convention put as angle 999."
	exit 1
fi

# delcare the array for parameters
declare -a arr 

# Compile the antip_scan.for code and create the
# executable in bin/ directory
gfortran -std=legacy src/antip_scan.for -o bin/antip_scan

if [ "$opt" = "-f" ]; then

	# input of script: parameter file (path)
	inputfile=$2

	while IFS= read -r line;
		do
			if [[ "$line" == "#"* ]]; then 	#line with # are not read
				true
			else 
				# echo "$line"
				# Save the line in an array
				IFS=" " read -r -a arr <<< "$line"
			fi
		done < "$2"

	if [ ${#arr[@]} -ne 12 ]; then
		echo "ERROR!"
		echo "Please, provide a file with the correct number of parameter (12)."
		echo "Execute: "
		echo "./p_scan.sh -h (or --help)"
		echo "for more info."
		exit 1
	else
		true
	fi
elif [ "$opt" = "-i" ]; then
	# save the arguments as an array
	arr=( "${@:2:13}" )
fi

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
		# ${arr%.*} convert from float to int, then
		# again a float with .0
		arrmom[i]=$( expr $((${arr%.*}+i*step)) ).0 
	done
	for mom in "${arrmom[@]:0:$N+1}"; do
		# Calculate the elastic and inelastic
		# cross sections
		./bin/antip_scan $mom "${arr[@]:1:12}" 
		echo "plab momentum "$mom
	done
	echo "Text file with results saved in results/ directory with name ${filename}.dat"
	echo "Creating the plot..."
	# plot with gnuplot of the cross sections
	plotGraph
	echo "Plot saved in fig/ directory with name ${filename}.png"
elif [ "${arr[10]}" = "ang" ]; then
	# Calculate the differential cross section
	# at specific angle and momentum
	./bin/antip_scan "${arr[@]:0:12}"
else
	echo "Please, provide a valid option for scan (11th parameter)."
	exit 1
fi
