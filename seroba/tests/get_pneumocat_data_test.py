import unittest
import os
import filecmp
import shutil
import tempfile
from seroba import get_pneumocat_data

modules_dir = os.path.dirname(os.path.abspath(get_pneumocat_data.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestGetPneumocatData(unittest.TestCase):

    def test__pneumocat_db_2_tsv(self):
        serogroup_dir = os.path.join(data_dir,'streptococcus-pneumoniae-ctvdb')
        os.remove(os.path.join(serogroup_dir,'06A.fasta'))
        expected_pneumocat_meta_tsv = os.path.join(data_dir,'expected_pneumocat_meta.tsv')
        out_file = os.path.abspath('pneumocat_db_test.tsv')
        get_pneumocat_data.GetPneumocatData._pneumocat_db_2_tsv(serogroup_dir,out_file)
        self.assertTrue(filecmp.cmp(out_file ,expected_pneumocat_meta_tsv), 'files are not equal')
        os.remove(os.path.join(serogroup_dir,'07A.fasta'))
        os.remove(os.path.join(serogroup_dir,'07B.fasta'))
