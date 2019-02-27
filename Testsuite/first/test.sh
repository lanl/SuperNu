#!/bin/bash
#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.

#-- test name
testname='first test'
#-- executable
binary='./supernu'

#-- path constants
testdir=$(cd $(dirname $0); pwd)
rundir=$1

#-- put files belonging to this test into place
cd $rundir || exit 1
ln -sf $testdir/* .

#-- run test
echo "TESTSUITE RUN: $testname"
./supernu || exit 1  #-- fail if unclean exit

#-- analyze output
for i in grd_temp flx_luminos tsp_time grd_fcoef; do
  diff -q "reference.$i" "output.$i" || exit 1
done
