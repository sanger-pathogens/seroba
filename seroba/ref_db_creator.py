import ariba
from seroba import kmc
import tempfile
import os
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pprint
import pyfastaq
import shutil

class Error (Exception): pass

class RefDbCreator:
    def __init__(self,out_dir,kmer_size):
        self.ref_fasta = os.path.join(out_dir,'reference.fasta')
        self.meta_data_tsv = os.path.join(out_dir,'meta.tsv')
        self.kmer_size = kmer_size
        self.out_dir = out_dir
        """try:
            fd = os.open(ref_fasta,'r')
        except :
            raise Error('Reference fasta file does not exits: '+ref_fasta )
        try:
            fd = os.open(meta_data_tsv,'r')
        except :
            raise Error('Meta data file does not exits: '+meta_data_tsv )"""

    def _split_meta_data2serogroup(self,serogroup):

        print(serogroup)
        record_dict = SeqIO.to_dict(SeqIO.parse(self.ref_fasta, "fasta"))
        serogroup_fasta = os.path.join(os.path.dirname(self.out_dir),'streptococcus-pneumoniae-ctvdb',serogroup+'.fasta')
        print(serogroup_fasta)
        print(os.path.isfile(serogroup_fasta))
        if os.path.isfile(serogroup_fasta):
            print("file does exist at this time")
            record_serogroup_dict = SeqIO.to_dict(SeqIO.parse(serogroup_fasta, "fasta"))
        else:

            print("no such file")

        self.temp_dir = tempfile.mkdtemp(prefix = 'temp_ariba', dir=os.getcwd())
        self.temp_fasta = os.path.join(self.temp_dir,'temp_fasta.fasta')
        self.temp_meta = os.path.join(self.temp_dir,'temp_meta.tsv')

        #read fasta
        #grep all lines with serogroup
        with open(self.meta_data_tsv,'r') as fobj:
            tsvin = csv.reader(fobj, delimiter='\t')
            with open (self.temp_fasta,'a') as fasta:
                with open(self.temp_meta,'a') as meta:
                    prev_row = ''
                    for row in tsvin:
                        if serogroup in row and row[0] != prev_row:
                            prev_row = row[0]
                            if 'ref' in row:
                                records = SeqRecord(record_dict[row[0]].seq,id = row[0],description='')
                                meta.write(row[0]+'\t0\t0\t.\t.\t.\n')
                            else:
                                records = SeqRecord(record_serogroup_dict[row[0]].seq,id = row[0],description='')
                                meta.write(row[0]+'\t1\t0\t.\t.\t.\n')
                            SeqIO.write(records, fasta, "fasta")

        #copy all seqs to temp fasta
        #format lines to ariba valid format
    @staticmethod
    def _create_cdhit_cluster_file(temp_dir,meta_tsv):
        with open(meta_tsv,'r') as fobj:
            tsvin = csv.reader(fobj, delimiter='\t')
            temp_cdhit_cluster = os.path.join(temp_dir,'cdhit_cluster')
            with open(temp_cdhit_cluster,'w') as wobj:
                ref_seqs = []
                prefs = {}
                genes = []
                for row in tsvin:
                    if row[1] == '0':
                        ref_seqs.append(row[0])
                    else:
                        pre = row[0].split('-')[0].split('_')[0]
                        if pre not in prefs:
                            prefs.update({pre:len(genes)})
                            genes.append([row[0]])
                        else:
                            genes[prefs[pre]].append(row[0])
                wobj.write('\t'.join(ref_seqs)+'\n')
                for gene in genes:
                    wobj.write('\t'.join(gene)+'\n')
        return temp_cdhit_cluster

    @staticmethod
    def _check_meta_data(meta_tsv,gene_fasta):
        if os.path.isfile(gene_fasta):
            record_serogroup_dict = SeqIO.to_dict(SeqIO.parse(gene_fasta, "fasta"))
        meta_data = []
        with open(meta_tsv,'r') as meta:
            tsvin = csv.reader(meta, delimiter='\t')
            for row in tsvin:
                if row[1] == '1':
                    seq = pyfastaq.sequences.Fasta(row[0], record_serogroup_dict[row[0]].seq)
                    got = seq.looks_like_gene()

                    if got == False:

                        meta_data.append([row[0],'0']+row[2:])
                        print ([row[0],'0']+row[2:])
                    else:
                        meta_data.append(row)
                else:
                    meta_data.append(row)
        with open(meta_tsv,'w') as meta:
            for elem in meta_data:
                meta.write('\t'.join(elem)+'\n')


    @staticmethod
    def _create_ariba_db(fasta_file,cluster_meta_data,cdhit_clusters,serogroup,out_dir):
        ariba_dir = os.path.join(out_dir,'ariba_db',serogroup)
        command = ['ariba prepareref','-f',fasta_file,'-m',cluster_meta_data,'--cdhit_clusters',cdhit_clusters, ariba_dir]
        os.system(' '.join(command))

    @staticmethod
    def _create_complete_cdhit_cluster(meta_data_tsv,out_dir):
        serogroup_dict = {}
        with open(meta_data_tsv,'r') as fobj:
            tsvin = csv.reader(fobj, delimiter='\t')
            for row in tsvin:
                if row[4] == 'ref':
                    if row[1] not in serogroup_dict:
                        serogroup_dict[row[1]] = [row[2]]
                    else:
                        serogroup_dict[row[1]].append(row[2])

        temp_cdhit_cluster = os.path.join(out_dir,'cdhit_cluster')
        serogroups = sorted(list(serogroup_dict.keys()))
        with open(temp_cdhit_cluster,'w') as cd :
            for serogroup in serogroups:
                cd.write('\t'.join(sorted(serogroup_dict[serogroup]))+'\n')
        return temp_cdhit_cluster


    def _create_complete_ariba_db(self,complete_cdhit_cluster):
        temp_dir = tempfile.mkdtemp(prefix = 'temp_ariba', dir=os.getcwd())
        temp_dir = os.path.join(temp_dir,'ariba')
        command = ['ariba prepareref','-f',self.ref_fasta,'--all_coding no --cdhit_clusters',complete_cdhit_cluster, temp_dir]
        os.system(' '.join(command))
        shutil.copyfile(os.path.join(temp_dir,'02.cdhit.clusters.tsv'),os.path.join(self.out_dir,'cd_cluster.tsv'))
        shutil.rmtree(temp_dir)


    def _create_kmc_db(self):

        kmer_db_dir = os.path.join(self.out_dir,'kmer_db')
        record_dict = SeqIO.to_dict(SeqIO.parse(self.ref_fasta, "fasta"))
        with open(self.meta_data_tsv) as meta:
            tsvin = csv.reader(meta, delimiter='\t')
            for row in tsvin:
                if row[4] == 'ref':
                    if not os.path.exists(os.path.join(kmer_db_dir,row[0])):
                        os.makedirs(os.path.join(kmer_db_dir,row[0]))
                    serotype_fasta = os.path.join(kmer_db_dir,row[0],row[0]+'.fasta')
                    with open(serotype_fasta,'w') as fasta:
                        records = SeqRecord(record_dict[row[0]].seq,id = row[0],description='')
                        SeqIO.write(records, fasta, "fasta")
                    kmc.run_kmc(serotype_fasta,self.kmer_size,os.path.join(kmer_db_dir,row[0]),row[0])

    @staticmethod
    def _read_meta_data_tsv(meta_data_tsv):
        meta_data_dict = {}
        with open(meta_data_tsv,'r') as fobj:
            tsvin = csv.reader(fobj, delimiter='\t')

            for row in tsvin:
                if row[4] != 'ref':
                    if row[1] not in meta_data_dict:
                        sub_dict = {'allele':{},'snps':{},'genes':{},'pseudo':{}}
                        meta_data_dict[row[1]]= sub_dict

                    if row[4] == 'allele':

                        if row[0].split('-')[0] not in meta_data_dict[row[1]]['allele']:
                            meta_data_dict[row[1]]['allele'].update({row[0].split('-')[0]:{row[2]:row[0]}})
                        else:
                            meta_data_dict[row[1]]['allele'][row[0].split('-')[0]].update({row[2]:row[0]})

                    elif row[4] == 'snps':
                        if  row[0].split('_')[0] not in meta_data_dict[row[1]]['snps']:
                            meta_data_dict[row[1]]['snps'].update({row[0].split('_')[0]:{row[5]:{row[2]:row[6]}}})
                        elif row[5] not in meta_data_dict[row[1]]['snps'][row[0].split('_')[0]]:
                            meta_data_dict[row[1]]['snps'][row[0].split('_')[0]][row[5]]= {row[2]:row[6]}
                        else:
                            meta_data_dict[row[1]]['snps'][row[0].split('_')[0]][row[5]].update({row[2]:row[6]})

                    elif row[4] == 'pseudo':
                        meta_data_dict[row[1]]['pseudo'].update({row[2]:{row[0]:row[5]}})
                    elif row[4] == 'genes':
                        meta_data_dict[row[1]]['genes'].update({row[2]:{row[0]:row[5]}})
                elif row[4] == 'ref' and row[1] not in meta_data_dict:
                    meta_data_dict[row[1]] =  {'allele':{},'snps':{},'genes':{},'pseudo':{}}

        for serogroup in meta_data_dict:
            dels = []
            for genetic_inf in meta_data_dict[serogroup]:
                if meta_data_dict[serogroup][genetic_inf] == {}:
                    dels.append(genetic_inf)
            for key in dels:
                meta_data_dict[serogroup].pop(key)


        return meta_data_dict


    def run(self):
        self.meta_dict = self._read_meta_data_tsv(self.meta_data_tsv)
        os.makedirs(os.path.join(self.out_dir,'ariba_db'))
        cdhit_cluster = RefDbCreator._create_complete_cdhit_cluster(self.meta_data_tsv,self.out_dir)
        self._create_complete_ariba_db(cdhit_cluster)
        for serogroup in self.meta_dict:

            self._split_meta_data2serogroup(serogroup)
            self.cdhit_clusters= self._create_cdhit_cluster_file(self.temp_dir,self.temp_meta)
            gene_fasta = serogroup_fasta = os.path.join(os.path.dirname(self.out_dir),'streptococcus-pneumoniae-ctvdb',serogroup+'.fasta')
            RefDbCreator._check_meta_data(self.temp_meta,gene_fasta)
            RefDbCreator._create_ariba_db(self.temp_fasta,self.temp_meta,self.cdhit_clusters,serogroup,self.out_dir)
            shutil.rmtree(self.temp_dir)
        self._create_kmc_db()
        kmer_info = os.path.join(self.out_dir,'kmer_size.txt')
        with open(kmer_info,'w') as fobj:
            fobj.write(self.kmer_size)
        pprint.pprint(self.meta_dict)
