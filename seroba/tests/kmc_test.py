import unittest
import os
import filecmp
import shutil
import tempfile
from seroba import external_progs, common, kmc

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
ext_progs = external_progs.ExternalProgs()
class TestKMC(unittest.TestCase):
    def test_run_kmc_fa(self):
        expected = os.path.join(data_dir,'expected_fa_database')
        sequence_file = os.path.join(data_dir,'test_build_kmc_database.fa')
        kmer_size = '51'
        tmp_dir = tempfile.mkdtemp(prefix = 'temp.kmc', dir=os.getcwd())
        got = kmc.run_kmc(sequence_file,kmer_size,tmp_dir,'count')
        trans_got = os.path.join(tmp_dir,'test_db')
        os.system(ext_progs.exe('kmc_tools')+' transform '+got+' histogram '+ trans_got)
        self.assertTrue(filecmp.cmp(trans_got,expected), 'files are equal')
        shutil.rmtree(tmp_dir)

    def test_run_kmc_fq(self):
        expected = os.path.join(data_dir,'test_build_kmc_read_hist')
        sequence_file = os.path.join(data_dir,'test_build_kmc_read.fq.gz')
        kmer_size = '51'
        tmp_dir = tempfile.mkdtemp(prefix = 'temp.kmc', dir=os.getcwd())
        got = kmc.run_kmc(sequence_file,kmer_size,tmp_dir,'count')
        trans_got = os.path.join(tmp_dir,'test_db')
        os.system(ext_progs.exe('kmc_tools')+' transform '+got+' histogram '+ trans_got)
        self.assertTrue(filecmp.cmp(trans_got,expected), 'files are equal')
        shutil.rmtree(tmp_dir)


    def test_run_kmc_intersect(self):
        expected = os.path.join(data_dir,'test_kmc_inter_hist')
        sequence_file = os.path.join(data_dir,'test_build_kmc_read')
        db = os.path.join(data_dir,'serotype_object','kmer_db','01','01')
        tmp_dir = (data_dir)
        got = kmc.run_kmc_intersect(sequence_file,tmp_dir,db)
        trans_got = os.path.join(tmp_dir,'test_db')
        os.system(ext_progs.exe('kmc_tools')+' transform '+got+' histogram '+ trans_got)
        self.assertTrue(filecmp.cmp(trans_got,expected), 'files are equal')
        os.remove(trans_got)


    def test_run_kmc_hist(self):
        expected = os.path.join(data_dir,'test_kmc_inter_hist')
        inter = os.path.join(data_dir,'inter')
        temp_dir = (data_dir)
        got = kmc.run_kmc_hist(inter,temp_dir)
        self.assertTrue(filecmp.cmp(got,expected), 'files are equal')
