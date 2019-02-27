#!/usr/bin/python2
#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
from version import revision_id, revision_number, date

fname = "version.inc"

with open(fname,'w') as f:
    f.write("      coderev_id = '" + revision_id + "'\n")
    f.write("      coderev_nr = '" + str(revision_number) + "'\n")
    f.write("      coderev_date = '" + date + "'\n")
