import argparse
from seroba import summary

def run(options):
    summ = summary.Summary(options.out_dir)
    summ.run()
