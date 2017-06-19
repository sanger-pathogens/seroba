import argparse
from seroba import summary


def run(options):
    summ = summary.summarise(options.out_dir)
    summ.run()
