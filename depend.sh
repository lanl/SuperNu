#!/bin/bash
#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
if [ $# != 0 ]; then
 files=$*
else
 files=*.f *.f90
fi

for i in $files; do
 fname=${i/%.o/.f}
 if [[ ! -f $fname ]]; then fname=${i/%.o/.f90}; fi
 [[ -s $fname ]] || continue
 dep=$(egrep '^ +(use |include )' $fname | tr '"' ' '|tr "'" ' '| sed -e "s/^cMPI//;s/ *,.*$//;s/ *!.*$//" | awk '{if(tolower($1) == "use"){print $2".mod"}else{print $2}}' | sort | uniq |tr "\n" " ")
 if [[ -n $dep ]]; then
  printf "${i/%.f*/.o}: $fname $dep\n"
 fi
done
