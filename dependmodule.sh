#!/bin/bash
#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
if [ $# != 0 ]; then
 files=$*
else
 return 1
fi

#
#-- collect the filenames belonging to each module
declare -A modnarr
for i in $files; do
 fname=${i/%.o/.f}
 if [[ ! -f $fname ]]; then fname=${i/%.o/.f90}; fi
 [[ -s $fname ]] || continue
 modnames=$(egrep '^ *module ' $fname |awk '{print $2}')
 for modname in $modnames; do
  modnarr[${modname}]=${i/%.f*/.o}
 done
done

#
for i in $files; do
#-- filename
 fname=${i/%.o/.f}
 if [[ ! -f $fname ]]; then fname=${i/%.o/.f90}; fi
 [[ -s $fname ]] || continue
#-- module name
 modname=$(egrep -m 1 '^ *module ' $fname |awk '{print $2}')
#-- extract dependencies
 depraw=($(egrep '^ +(use |include )' $fname | tr '"' ' '|tr "'" ' '| sed -e "s/ *,.*$//;s/ *!.*$//"))
#-- translate module dependencies to module filename
 deps=()
 for ((j=0;j<${#depraw[*]};j=j+2)); do
  typ=${depraw["$j"],,}  #-- lowercase
  dep=${depraw["$((j+1))"]}
#-- include file
  if [[ $typ == 'include' ]]; then
   deps+=($dep)
#-- module
  else
   dep2=${modnarr["$dep"]}
   [[ -z $dep2 ]] && printf "\nin file ($fname) unknown module: $dep\n" >&2 && exit 1
   deps+=($dep2)
  fi
 done
#-- uniq
 dep=$(echo ${deps[@]} |tr ' ' '\n' |sort |uniq |tr '\n' ' ')

 printf "${i/%.f*/.o} $modname.mod: $fname $dep\n"
done
