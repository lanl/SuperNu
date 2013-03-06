#!/bin/bash

#default values
FC='mpixlf2003_r'
optimize="-O3 -qfloat=fltint:rsqrt -qnohot -qinitauto=fff00000 -qalign=4k -qcache=auto -qdpc=e -qmaxmem=16384 -qnozerosize -qsigtrap -qstrict -qtune=auto -qhalt=S -qnolm"
debug=""
omp="-qsmp=omp"

#override default
for a in $@; do
 if [[ $a == 'ompOff' ]]; then
  FC=${FC/_r/}
  omp="-qsmp=noomp"
 elif [[ $a == 'mpiOff' ]]; then FC=${FC/mpi/}
 elif [[ $a == 'optOff' ]]; then optimize="-O0"
 elif [[ $a == 'debug' ]]; then
  optimize="-O0"
  debug="-C -g -qnoextchk -qflttrap=enab:ov:zero:inv"
 else
  echo "option unknown: $a"
  exit
 fi
done

cat <<END >Makefile.compiler
#==== xlf2003 (current version: 11.1 - 15/jun/2011) ====
FC = $FC
FCFLAGS = $omp $optimize $debug
FCFLAGS += #-qarch=450 #-qarch=450d
FCFLAGS_FF = -qfixed
LDFLAGS =

AR = ar
RANLIB=ranlib
END

if [[ $FC == ${FC/mpi/} ]]; then
 ln -f mpimod_ser.f mpimod.f
else
 ln -f mpimod_mpi.f mpimod.f
fi

cat Makefile.compiler
