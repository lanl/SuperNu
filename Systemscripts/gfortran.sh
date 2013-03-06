#!/bin/bash

#default values
FC="gfortran"
omp="-fopenmp"
optimize="-O3"
debug=""

#override default
for a in $@; do
 if [[ $a == 'ompOff' ]]; then
  omp=''
 elif [[ $a == 'mpi' ]]; then
  FC="mpif90"
 elif [[ $a == 'optOff' ]]; then
  optimize="-O0"
 elif [[ $a == 'debug' ]]; then
  optimize=""
  debug="-g -O0 -Wall -fcheck=all"
 else
  echo "option unknown: $a"
  exit
 fi
done

cat <<END >Makefile.compiler
#==== gfortran (current version: gcc 4.5) ====
FC = $FC
FCFLAGS = $optimize $debug
#FCFLAGS += -m64
FCFLAGS += $omp
FCFLAGS_FF = -ffixed-form
LDFLAGS = $omp

AR = ar
RANLIB=ranlib
END

if [[ $FC == 'gfortran' ]]; then
 ln -f mpimod_ser.f mpimod.f
else
 ln -f mpimod_mpi.f mpimod.f
fi

cat Makefile.compiler
