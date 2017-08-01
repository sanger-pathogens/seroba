import argparse
from seroba import serotyping
import sys

def run(options):
    if (options.read1 == options.read2):
        print('Same file provided for forwards and reverse reads. Cannot continue', file=sys.stderr)
        sys.exit(1)

    if (options.read1.rsplit('_',1)[0] != options.read2.rsplit('_',1)[0]):
        print('Names for forwards and reverse reads does not match. Cannot continue', file=sys.stderr)
        sys.exit(1)
    sero = serotyping.Serotyping(options.databases,
    options.read1,
    options.read2,
    options.prefix,
    clean = (not options.noclean))
    cov = options.coverage
    sero.run()
