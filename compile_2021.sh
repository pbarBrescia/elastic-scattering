#!/usr/bin/env bash
if [ ! -d bin ]; then
  mkdir bin;
fi

gfortran src/antiprova_2021.for -o bin/antip_new 2> /dev/null  ## suppressing a fortran compiling warning
g++ src/coulomb.cpp -o bin/scattering
echo "antip_new and scattering programs compiled (executable files in the bin directory)"

