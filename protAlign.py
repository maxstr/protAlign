from fabric.api import local, env, settings
from itertools import chain


def tmAlign(file1, file2):
    with settings(capture=True):
        return local("%s %s %s" % (env.tmPath, file1, file2))



