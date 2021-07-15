#! /usr/bin/env bash

filename=int_cs
outfile=./results/int_cs.dat 

gfortran -std=gnu src/antip_scan.for -o bin/antip_scan

opt=$1

if [ "$opt" = "mom" ]; then
	rm $outfile 2>/dev/null
	echo "Executing the scan on momentum..."
	for mom in {50..400..25} ; do
		./bin/antip_scan $mom.0 40.078 20.0 30.0 150.0 1.25 1.25 1.25 0.6 0.5 $opt 1.0 #> /tmp/file #2>/dev/null
		# awk -v m=$mom '{print m,$0}' /tmp/file >> $outfile
		echo "plab momentum "$mom
	done
	echo "Text file with results saved in results/ directory."
	echo "Creating the plot..."
	gnuplot <<-EOF
		set xlabel "plab (MeV/c)"
		set ylabel "{/Symbol s} (mb)"
		set title "Integral cross section (elastic and reaction)"
		set term png
		set output "fig/${filename}.png"
		plot "${outfile}" using 1:2 title "elastic" with linespoints, \
		"${outfile}" using 1:3 title "reaction" with linespoints, \
		"${outfile}" using 1:4 title "elastic+reaction" with linespoints
	EOF
	echo "Plot saved in fig/ directory."
elif [ "$opt" = "ang" ]; then
	mom=$2
	ang=$3
	./bin/antip_scan $mom 40.078 20.0 30.0 150.0 1.25 1.25 1.25 0.6 0.5 $opt $ang
else
	echo "Please, provide 1 or 3 arguments: "
	echo "	Usage: ./p_scan.sh <option> (<argument 1> <argument 2>) "
	echo
	echo "	option \"mom\" -> ./p_scan.sh mom : scan from 50 to 400 MeV/c (integral cross section)"
	echo "	option \"ang\" -> ./p_scan.sh ang <momentum(MeV/c)> <angle(deg)> : diff. cross section for specific momentum and angle"
fi

