import unittest
import os
import filecmp
import shutil
import tempfile
from seroba import ref_db_creator

modules_dir = os.path.dirname(os.path.abspath(ref_db_creator.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestRefDbCreator(unittest.TestCase):
    def test__read_meta_data_tsv(self):
        expected = {
            '6':
                {'allele':
                    {'wciN':
                        {'06A': 'wciN-1',
                        '06B': 'wciN-1',
                        '06C': 'wciN-2',
                        '06D': 'wciN-2',
                        '06E': 'wciN-3'}},
                'snps':
                        {'wciP':
                            {'583':
                                {'06A': 'AGT',
                                '06B': 'AAT',
                                '06C': 'AGT',
                                '06D': 'AAT' }}}},
            '7':
                {'snps':
                    {'wcwK':
                        {'46':
                            {'07B': 'GAT',
                            '07C': 'GGT',
                            '40': 'AAT',},
                         '145':
                            {'07B': 'CTT',
                            '07C': 'CTT',
                            '40': 'TTT' },
                        '385': {
                            '07B': 'TTT',
                            '07C': 'TGT',
                            '40': 'ACT'},
                        '487':
                            {'07B':'ACT',
                            '07C':'GCT',
                            '40':'ACT' },
                         '706':
                            {'07B': 'CAT',
                            '07C': 'CAT',
                            '40': 'TAT'},
                         '880':
                            {'07B': 'CTT',
                             '07C': 'CTT',
                             '40': 'TTT'},
                         '928':
                            {'07B': 'AAT',
                            '07C': 'AAT',
                            '40': 'AGT'},
                        '937':
                            {'07B': 'GCA',
                            '07C': 'GAA',
                            '40': 'GAA'},
                        '946':
                            {'07B': 'GGT',
                            '07C': 'GGT',
                            '40': 'GAT'}}}},
                '9':
                    {'genes':
                        {'09A':{'wcjD': '1'},
                         '09L':{'wcjD': '0'},
                         '09N':{'wcjD': '0'},
                         '09V':{'wcjD': '1'}},
                     'pseudo': {
                        '09A': {'wcjE': '1'},
                        '09V': {'wcjE': '0'}}}}
        self.meta_data_tsv =  os.path.join(data_dir,'test_dict_meta_data.tsv')
        got = ref_db_creator.RefDbCreator._read_meta_data_tsv(self.meta_data_tsv)
        self.assertEqual(expected,got)


    def test__split_meta_data2serogroup(self):
        ref_dir = os.path.join(data_dir,'test__split_meta_data2serogroup')
        kmer_size = '51'
        r = ref_db_creator.RefDbCreator(ref_dir,kmer_size)
        expected_temp_fasta = os.path.join(ref_dir,'expected_ref_fasta.fa')
        expected_temp_meta = os.path.join(ref_dir,'expected_meta_data.tsv')
        serogroup = '06A'
        r._split_meta_data2serogroup(serogroup)
        self.assertTrue(filecmp.cmp(r.temp_fasta,expected_temp_fasta), 'files are not equal')
        self.assertTrue(filecmp.cmp(r.temp_meta,expected_temp_meta), 'files are not equal')
        shutil.rmtree(r.temp_dir)

    def test__create_cdhit_cluster_file(self):
        expected_cluster = os.path.join(data_dir,'expected_cdhit_cluster.tsv')
        meta_data = os.path.join(data_dir,'meta.tsv')
        temp_dir = tempfile.mkdtemp(prefix = 'temp_test', dir=os.getcwd())
        ref_db_creator.RefDbCreator._create_cdhit_cluster_file(temp_dir,meta_data)
        shutil.rmtree(temp_dir)

    def test__check_meta_data(self):
        expected_meta_data =  os.path.join(data_dir,'test_check_meta_data.tsv')
        meta_data =  os.path.join(data_dir,'meta6.tsv')
        gene_fasta = os.path.join(data_dir,'06A.fasta')
        ref_db_creator.RefDbCreator._check_meta_data(meta_data,gene_fasta)
        self.assertTrue(filecmp.cmp(meta_data,expected_meta_data ), 'files are not equal')

    def test__create_ariba_db(self):

        fasta_file =  os.path.join(data_dir,'test__create_ariba_db','expected_ref_fasta.fa')
        cluster_meta_data = os.path.join(data_dir,'test__create_ariba_db','expected_meta_data.tsv')
        cdhit_clusters = os.path.join(data_dir,'test__create_ariba_db','cdhit_cluster.tsv')
        out_dir = 'temp'
        os.makedirs('temp/ariba_db')
        serogroup = '06A'
        ref_db_creator.RefDbCreator._create_ariba_db(fasta_file,cluster_meta_data,cdhit_clusters,serogroup,out_dir)
        db_files = os.listdir('temp/ariba_db/06A')
        for i in range(len(db_files)):
            with self.subTest (i= i):
                exp = os.path.join(data_dir,'ariba_db',db_files[i])
                print(exp)
                self.assertTrue(os.path.isfile(exp) , 'files does not exists')
        shutil.rmtree(out_dir)

    def test__create_complete_cdhit_cluster(self):
        out_dir =  tempfile.mkdtemp(prefix = 'temp_test', dir=os.getcwd())
        expected = os.path.join(data_dir,'expected_complete_cdhit_cluster')
        meta_data = os.path.join(data_dir,'test_complete_meta.tsv')
        got = ref_db_creator.RefDbCreator._create_complete_cdhit_cluster(meta_data,out_dir)
        self.assertTrue(filecmp.cmp(got,expected), 'files are not equal')
        shutil.rmtree(out_dir)

    def test___create_complete_ariba_db(self):
        pass
    def test__create_kmc_db(self):
        ref_dir= os.path.join(data_dir,'test__create_kmc_db')
        kmer_size = '51'
        r = ref_db_creator.RefDbCreator(ref_dir,kmer_size)
        r._create_kmc_db()

        self.assertTrue(os.path.isfile(os.path.join(ref_dir,'kmer_db','06A','06A.kmc_pre')) , 'files does not exists')
        self.assertTrue(os.path.isfile(os.path.join(ref_dir,'kmer_db','06A','06A.kmc_suf')) , 'files does not exists')
        self.assertTrue(os.path.isfile(os.path.join(ref_dir,'kmer_db','06B','06B.kmc_pre')) , 'files does not exists')
        self.assertTrue(os.path.isfile(os.path.join(ref_dir,'kmer_db','06B','06B.kmc_suf')) , 'files does not exists')
        shutil.rmtree(os.path.join(ref_dir,'kmer_db'))


    def test_run(self):

        pass
