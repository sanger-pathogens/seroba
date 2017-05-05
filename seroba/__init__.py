from pkg_resources import get_distribution

try:
    __version__ = get_distribution('seroba').version
except:
    __version__ = 'local'


__all__ = [
    'kmc',
    'common',
    'external_progs',
    'tasks',
]

from seroba import *
