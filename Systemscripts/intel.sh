#!/bin/bash

#default values
FC="ifort"
omp="-openmp"
optimize="-O2 -msse3"
debug=""

#override default
for a in $@; do
 if [[ $a == 'ompOff' ]]; then
  omp=''
 elif [[ $a == 'mpi' ]]; then
  FC="mpif90"
 elif [[ $a == 'mpiftn' ]]; then
  FC="ftn"
 elif [[ $a == 'optOff' ]]; then
  optimize="-O0"
 elif [[ $a == 'debug' ]]; then
  optimize=""
  debug="-O0 -g -check -debug -traceback"
 else
  echo "option unknown: $a"
  exit
 fi
done

cat <<END >Makefile.compiler
#==== intel fortran compiler (current version: ifort 11.1 20091012) ====
FC = $FC
FCFLAGS = $optimize $debug
FCFLAGS += -m64 $omp
LDFLAGS = $omp

AR = ar
RANLIB=ranlib
END

if [[ $FC == 'ifort' ]]; then
 ln -f mpimod_ser.f mpimod.f
else
 ln -f mpimod_mpi.f mpimod.f
fi

cat Makefile.compiler
