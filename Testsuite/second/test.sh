#!/bin/bash

#-- test name
testname='static grid test'
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
for i in temp Lum tsp_time gas_fcoef; do
  diff -q "reference.$i" "output.$i" || exit 1
done
