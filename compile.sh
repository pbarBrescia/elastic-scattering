#!/usr/bin/env bash
if [ ! -d bin ]; then
  mkdir bin;
fi

gfortran src/antiprova.for -o bin/antip #2> /dev/null  ## suppressing a fortran compiling warning
gfortran src/antip_scan.for -o bin/antip_scan
g++ src/coulomb.cpp -o bin/scattering
echo "antip, antip_scan and scattering programs compiled (executable files in the bin directory)"

