import unittest
import os
import filecmp
import shutil
import tempfile
from seroba import serotyping

modules_dir = os.path.dirname(os.path.abspath(serotyping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSerotyping(unittest.TestCase):

    def test_serotype_2_cluster(self):
        expected1 = {
            'NT' : 'NT',
            '06A' : 'cluster',
            '06B' : 'cluster',
            '06C' : 'cluster',
            '06D' : 'cluster',
            '06E' : 'cluster',
            '07A' : 'cluster_1',
            '07F' : 'cluster_1'
            }
        expected2 = {
            'NT' : 1,
            'cluster' : 5,
            'cluster_1' : 2
            }
        expected3 = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F']
        }

        cd_hit_cluster_file = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        got1, got3, got2 =serotyping.Serotyping._serotype_2_cluster(cd_hit_cluster_file)
        self.assertEqual(expected1,got1)
        self.assertEqual(expected2,got2)
        self.assertEqual(expected3,got3)

    def test_run_kmc(self):
        expected = '06A'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'06B_1.fq.gz')
        bw_read = os.path.join(data_dir,'06B_2.fq.gz')
        meta_data = os.path.join(data_dir,'expected_pneumocat_meta.tsv')
        prefix = '06B_1'
        kmer_size = '51'
        kmer_db = 'kmer_db'
        ariba_cluster_db = 'ariba_db'
        reference_fasta =os.path.join(data_dir,'serotype_object','reference.fasta')
        s = serotyping.Serotyping( refs_dir,fw_read, bw_read, prefix)
        s._run_kmc()
        self.assertEqual(expected,s.best_serotype)

    def test_01_run_ariba_on_cluster(self):
        print('test ariba')
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'06B_1.fq.gz')
        bw_read = os.path.join(data_dir,'06B_2.fq.gz')
        meta_data = os.path.join(data_dir,'expected_pneumocat_meta.tsv')
        prefix = os.path.join(data_dir,'06B_1')
        reference_fasta =os.path.join(data_dir,'serotype_object','reference.fasta')
        cluster = 'cluster'
        s = serotyping.Serotyping( refs_dir,fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F']
        }
        s._run_ariba_on_cluster(cluster)

    def test__02_serotype6(self):
        print('test serotype 6')
        expected = '06B'
        assemblie_file = os.path.join(data_dir,'06B_1','assembled_genes.fa')
        cluster_db_path = os.path.join(data_dir,'ariba_cluster_db','06A')
        report_file = os.path.join(data_dir,'06B_1','report.tsv')
        got = serotyping.Serotyping.serotype6(assemblie_file,report_file)
        self.assertEqual(expected,got)


    def test__03_prediction_06(self):
        print('test pred 6')
        expected = '06B'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'06B_1.fq.gz')
        bw_read = os.path.join(data_dir,'06B_2.fq.gz')
        prefix = os.path.join(data_dir,'06B_1')
        cluster = 'cluster'
        assemblie_file = os.path.join(prefix,'assemblies.fa')
        s = serotyping.Serotyping(refs_dir,fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F']
        }
        s.best_serotype = '06A'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)
        shutil.rmtree(prefix)

    def test__prediction_09V(self):
        expected = '09V'

        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'09V')
        kmer_size = '51'
        cluster = 'cluster_2'
        assemblie_file = os.path.join(data_dir,'09V','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F'],
            'cluster_2':['09A','09L','09N','09V']        }
        s.best_serotype = '09N'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2,
        'cluster_2': 4
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)

    def test_mixed_sample_15B_C(self):
        expected = '15B/15C'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'15B_C','15B_C_1.fq')
        bw_read = os.path.join(data_dir,'15B_C','15B_C_2.fq')
        prefix = os.path.join(data_dir,'15B_C')
        cluster = 'cluster_3'
        assemblie_file = os.path.join(data_dir,'15B_C','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F'],
            'cluster_2':['10A','10B'],
            'cluster_3':['15A','15B','15C','15F']
                     }
        s.best_serotype = '15C'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2,
        'cluster_2': 2,
        'cluster_3': 4
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)

    def test_09N(self):
        expected = '09N'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'ERR1438851')
        kmer_size = '51'
        cluster = 'cluster_2'
        assemblie_file = os.path.join(data_dir,'09V','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F'],
            'cluster_2':['09A','09L','09N','09V']        }
        s.best_serotype = '09N'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2,
        'cluster_2': 4
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)

    def test_09L(self):
        expected = '09L'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'ERR1440275')
        kmer_size = '51'
        cluster = 'cluster_2'
        assemblie_file = os.path.join(data_dir,'ERR1440275','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F'],
            'cluster_2':['09A','09L','09N','09V']        }
        s.best_serotype = '09N'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2,
        'cluster_2': 4
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)

    def test_09N(self):
        expected = '09N'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'ERR1438851')
        kmer_size = '51'
        cluster = 'cluster_2'
        assemblie_file = os.path.join(data_dir,'ERR1438851','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F'],
            'cluster_2':['09A','09L','09N','09V']        }
        s.best_serotype = '09N'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2,
        'cluster_2': 4
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)

    def test_11A(self):
        expected = '11A'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'ERR1438910')
        kmer_size = '51'
        cluster = 'cluster_2'
        assemblie_file = os.path.join(data_dir,'ERR1438910','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F'],
            'cluster_2':['11A','11B','11C','11D','11F']        }
        s.best_serotype = '11A'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2,
        'cluster_2': 5
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)

    def test_33F(self):
        expected = '33F'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'ERR1438894')
        kmer_size = '51'
        cluster = 'cluster_2'
        assemblie_file = os.path.join(data_dir,'ERR1438894','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07A','07F'],
            'cluster_2':['33A','33F','37']        }
        s.best_serotype = '33F'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 2,
        'cluster_2': 3
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)
    def test_07C(self):
        expected = '07C'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'ERR1439287')
        kmer_size = '51'
        cluster = 'cluster_1'
        assemblie_file = os.path.join(data_dir,'ERR1439287','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07B','07C','40'],
            'cluster_2':['33A','33F','37']        }
        s.best_serotype = '07C'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 3,
        'cluster_2': 3
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)

    """def test_15B(self):
        expected = '15B'
        refs_dir = os.path.join(data_dir,'serotype_object')
        cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join(data_dir,'ERR1439297')
        kmer_size = '51'
        cluster = 'cluster_3'
        assemblie_file = os.path.join(data_dir,'ERR1439297','assemblies.fa')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.cluster_serotype_dict = {
            'cluster' : ['06A','06B','06C','06D','06E'],
            'NT' : ['NT'],
            'cluster_1' : ['07B','07C','40'],
            'cluster_2':['33A','33F','37'],
            'cluster_3':['15A','15B','15C','15F']         }
        s.best_serotype = '07C'
        s.cluster_count = {
        'cluster': 5,
        'NT': 1,
        'cluster_1': 3,
        'cluster_2': 3,
        'cluster_3': 4
        }
        s._prediction(assemblie_file,cluster)
        self.assertEqual(expected,s.sero)"""
    def test_11A(self):
         expected = '11A'
         refs_dir = os.path.join(data_dir,'serotype_object')
         cd_cluster = os.path.join(data_dir,'test_serotype_2_cluster.tsv')
         fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
         bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
         prefix = os.path.join(data_dir,'ERR1439321')
         kmer_size = '51'
         cluster = 'cluster_2'
         assemblie_file = os.path.join(data_dir,'ERR1439321','assemblies.fa')
         s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
         s.cluster_serotype_dict = {
             'cluster' : ['06A','06B','06C','06D','06E'],
             'NT' : ['NT'],
             'cluster_1' : ['07A','07F'],
             'cluster_2':['11A','11B','11C','11D','11F']        }
         s.best_serotype = '11A'
         s.cluster_count = {
         'cluster': 5,
         'NT': 1,
         'cluster_1': 2,
         'cluster_2': 3
         }
         s._prediction(assemblie_file,cluster)
         self.assertEqual(expected,s.sero)
    def full_run_09V(self):
        expected = '09V'
        refs_dir = os.path.join(data_dir,'serotype_object')
        fw_read = os.path.join(data_dir,'09V','09V_1.fq.gz')
        bw_read = os.path.join(data_dir,'09V','09V_2.fq.gz')
        prefix = os.path.join('09V')
        s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
        s.run()
        self.assertEqual(expected,s.sero)

    """def full_run_23F(self):
       expected = '23F'
       refs_dir = os.path.join(data_dir,'serotype_object_23F')
       fw_read = os.path.join(data_dir,'23F','23F_cover_22_1.fastq.gz')
       bw_read = os.path.join(data_dir,'23F','23F_cover_22_2.fastq.gz')
       prefix = os.path.join('23F')
       s = serotyping.Serotyping(refs_dir, fw_read, bw_read, prefix)
       s.run()
       self.assertEqual(expected,s.sero)"""	
