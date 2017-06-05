import argparse
from seroba import serotyping


def run(options):
    sero = serotyping.Serotyping(options.databases,
    options.read1,
    options.read2,
    options.prefix)
    sero.run()
