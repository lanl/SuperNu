#!/bin/bash
#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
########################################################################
# Symbollically mirror a directory tree with source supernu files.
# The mercurial repo directory .hg is excluded.
# (DR van Rossum, 2016/01/08)
########################################################################

#-- test input
[ $# -lt 1 ] && { echo "Usage: $0 src_dir [dst_dir]" ; exit 1;}

#-- source dir
[ ! -d "$1" ] && { echo "$1 is not a valid directory"; exit 1;}
src=$(readlink -f $1)

#-- destination dir
dst=${2:-.}
[ ! -d "$dst" ] && { echo "does not exist: $dst"; exit 1;}
dst=$(readlink -f $dst)

#-- create destination dir structure
while IFS= read -r -d '' line; do
  dir="$dst${line#$src}"
  mkdir -pv "$dir"
done < <(find "$src" -mindepth 1 -path "$src"/.hg -prune -o -type d -print0)

#-- link files in place
while IFS= read -r -d '' line; do
  file="$dst${line#$src}"
  ln -s "$line" "$file"
done < <(find "$src" -path "$src"/.hg\* -prune -o -type f -print0 -o -type l -print0)
