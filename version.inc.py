#!/usr/bin/python2
from version import revision_id, revision_number, date

fname = "version.inc"

line1 = "      data coderev_id, coderev_nr, coderev_date, build_date /\n"

with open(fname,'w') as f:
    f.write(line1)
    f.write("     &  '" + revision_id + "',\n")
    f.write("     &  '" + str(revision_number) + "',\n")
    f.write("     &  '" + date + "',\n")
