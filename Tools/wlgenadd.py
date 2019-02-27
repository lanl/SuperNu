#!/usr/bin/python2
#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.

# open input.wlgrid file
f = open('input.wlgrid','r+')
print "File name: ", f.name

# obtain number of file lines
lnum = 0
while f.readline() != "":
    lnum+=1
## strip known header
lnum-=4
print lnum

# construct grid
## set min, max values and number of bins
ng = 10
wlmin = 1e-6
wlmax = 32e-5
wlgrid = []
for i in xrange(ng+1):
    wlgrid.append(wlmin*(wlmax/wlmin)**(float(i)/float(ng)))
print wlgrid
# convert grid values to string
str1 = str(wlgrid)[1:-1]+"\n"
str2 = str(lnum)
str3 = str2+" "+str(ng)+" "+str1

f.seek(-1,2)
f.write(str3)
f.write("#")

# close input.wlgrid file
f.close()
