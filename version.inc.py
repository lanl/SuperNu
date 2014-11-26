#!/usr/bin/python2
from version import revision_id, revision_number, date

fname = "version.inc"

with open(fname,'w') as f:
    f.write("      data coderev_id /'" + revision_id + "'/\n")
    f.write("      data coderev_nr /'" + str(revision_number) + "'/\n")
    f.write("      data coderev_date /'" + date + "'/\n")
