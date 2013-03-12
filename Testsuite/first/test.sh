#!/bin/bash

#-- test name
testname='first test'

#-- executable
binary='./supernu'

#-- path constants
curdir=$(pwd)
rundir=$curdir/Run
testdir=$(cd $(dirname $0); pwd)

#-- put files belonging to this test into place
cd $rundir
ln -sf $testdir/* .

#-- run test
echo "TESTSUITE RUN: $testname"
./supernu || exit 1  #-- fail if unclean exit

#-- analyze output
for i in temp Lum tsp_time gas_fcoef; do
  diff -q "reference.$i" "output.$i" || exit 1
done
