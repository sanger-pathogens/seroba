import unittest
import os
import filecmp
from seroba import summary

modules_dir = os.path.dirname(os.path.abspath(summary.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSummary(unittest.TestCase):
    def test__summary(self):
        expected = os.path.join(data_dir,'exp_summary.tsv')
        got=summary.summarise(os.path.join(data_dir,'summ'))
        self.assertTrue(filecmp.cmp('summary.tsv',expected), 'files are equal')
        os.remove('summary.tsv')
