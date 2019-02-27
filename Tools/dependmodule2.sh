#!/bin/bash
#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.

#---
# SuperNu module dependency finder for Mac
# Disclaimer: All modules must have "mod" at end of name.
# Run ./dependmodule2.sh <directory> and
# copy to binary directory.
# Make sure mpimod.f exists.
#---

#-- specify module directory
if [ -z $1 ]; then
    echo "No module directory specified"
    exit 1
fi

echo
echo "Disclaimer: All modules must have \"mod\" at end of name."
echo

#-- remove existing Makefile.dependmod
if [ -f Makefile.dependmod ]; then
    rm Makefile.dependmod
fi

#-- find module files
declare mods
l=0
for i in $(ls -l $1 | awk '/mod\.f(90)? / {print $9}'); do
    mods[$l]=$(echo $i)
    ((l+=1))
done
n=$((${#mods[@]}-1))

#-- strip fortran file endings
declare modns
for l in $(seq 0 $n); do
    modns[$l]=$(echo ${mods[$l]} | sed 's/\.f$//g' | sed 's/\.f90$//g')
done

#-- find and add dependencies to Makefile.dependmod
for l in $(seq 0 $n); do
    supmods=$(egrep '(^ *use|^ *include)' "$1/${mods[$l]}" | \
		     awk '{print $2}' | sed 's/[,\!].*//g' \
		     |sort |uniq |sed 's/mod$/mod.o/g')
    echo "${modns[$l]}.o ${modns[$l]}.mod: ${mods[$l]}" $supmods \
	 >> Makefile.dependmod
done
