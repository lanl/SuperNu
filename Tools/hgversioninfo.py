"""Write a version info .py file for your application to use at runtime.
"""

import os
import time

from mercurial import util


TEMPLATE = '''\
revision_id = %(revision_id)r
revision_number = %(revision_number)r
branch = %(branch)r
tags = %(tags)r
date = %(date)r
'''

def hook(ui, repo, node=None, **params):
    conf = _load_config(ui)
    changeset = repo[node]
    if not changeset.node():
        # this can happen when this function is called as an extension in a
        # working directory. in that case get the first parent changeset.
        changeset = changeset.parents()[0]
    versionfile = open(os.path.join(repo.root, conf['version_file_path']), 'w')
    _write_version_info(changeset, versionfile)

def _write_version_info(changeset, versionfile):
    template_vars = {
        'revision_id':changeset.node().encode('hex'),
        'revision_number':changeset.rev(),
        'branch':changeset.branch(),
        'tags':changeset.tags(),
        'date':time.ctime(changeset.date()[0]),
    }
    versionfile.write(TEMPLATE % template_vars)

def _load_config(ui):
    conf = dict(ui.configitems('hgversioninfo'))
    if not conf:
        raise util.Abort('You need a [hgversioninfo] section in your config')
    if not conf.get('version_file_path'):
        msg = ('You need to define the version_file_path in your '
                '[hgversioninfo] config section.')
        raise util.Abort(msg)
    return conf
