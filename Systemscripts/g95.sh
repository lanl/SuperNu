#!/bin/bash

#default values
FC='g95'
optimize="-O3"
debug=""

#override default
for a in $@; do
 if [[ $a == 'ompOff' ]]; then omp=''
 elif [[ $a == 'mpi' ]]; then
  FC='mpif90'
 elif [[ $a == 'mpiftn' ]]; then
  FC="ftn"
 elif [[ $a == 'optOff' ]]; then optimize="-O0"
 elif [[ $a == 'debug' ]]; then
  optimize=""
  debug="-g -O0 -Wall -fstatic -fbounds-check -ftrace=full -fpointer=invalid -fwrapv"
 else
  echo "option unknown: $a"
  exit
 fi
done

cat <<END >Makefile.compiler
#==== g95 (current version: 0.92 - 24/jun/2009) ====
FC = $FC
FCFLAGS = $optimize $debug
#FCFLAGS += -m64
FCFLAGS_FF = -ffixed-form
LDFLAGS =

AR = ar
RANLIB=ranlib
END

if [[ $FC == 'g95' ]]; then
 ln -f mpimod_ser.f mpimod.f
else
 ln -f mpimod_mpi.f mpimod.f
fi

cat Makefile.compiler
