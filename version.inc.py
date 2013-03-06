#!/usr/bin/python2
from version import revision_id, revision_number, date

fname = "version.inc"

line1 = "      data coderev_id, coderev_nr, coderev_date /\n"

try:
    f = open(fname,'r')
except:
    old_id = '-1'
else:
    dummy = f.readline()
    old_id = f.readline().split("'")[1]
    old_number = f.readline().split("'")[1]
    old_date = f.readline().split("'")[1]

if old_id != revision_id:
    print 'updating version.inc'
    with open(fname,'w') as f:
        f.write(line1)
        f.write("     &  '" + revision_id + "',\n")
        f.write("     &  '" + str(revision_number) + "',\n")
        f.write("     &  '" + date + "' /\n")
