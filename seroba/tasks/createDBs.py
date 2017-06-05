import argparse
from seroba import ref_db_creator


def run(options):
    ref_db = ref_db_creator.RefDbCreator(
      options.out_dir,
      options.kmer_size
          )
    ref_db.run()
