import argparse
from seroba import serotyping


def run(options):
    sero = serotyping.Serotyping(options.databases,
    options.fw_reads,
    options.bw_reads,
    options.prefix)
    sero.run()
