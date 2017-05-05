import unittest
import os
from seroba import common

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestCommon(unittest.TestCase):
    def test_detect_sequence_format(self):
        '''detect_sequence_format'''
        expected_fa = 'fa'
        expected_fq = 'fq'
        expected_fq_gz = 'fq'
        fa_file = os.path.join(data_dir, 'test_detect_sequence_format.fa')
        fq_file = os.path.join(data_dir, 'test_detect_sequence_format.fq')
        fq_gz_file = os.path.join(data_dir, 'test_detect_sequence_format.fq.gz')
        wrong_file = os.path.join(data_dir, 'test_detect_wrong_format.txt')
        got_fa = common.detect_sequence_format(fa_file)
        got_fq = common.detect_sequence_format(fq_file)
        got_fq_gz = common.detect_sequence_format(fq_gz_file)
        with self.assertRaises(common.Error):
            wrong_format = common.detect_sequence_format(wrong_file)
        self.assertEqual(expected_fa,got_fa)
        self.assertEqual(expected_fq,got_fq)
        self.assertEqual(expected_fq_gz,got_fq_gz)
