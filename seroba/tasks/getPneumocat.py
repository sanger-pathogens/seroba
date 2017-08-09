import argparse
from seroba import get_pneumocat_data


def run(options):
    pneumo = get_pneumocat_data.GetPneumocatData(options.database_dir)
    pneumo.run()
