import os
import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import shutil
import tempfile
import tarfile
class Error (Exception): pass

class GetPneumocatData:
    def __init__(self,out_dir):
        self.out_dir = out_dir
        self.pneumocat_release_link = 'https://github.com/phe-bioinformatics/PneumoCaT/archive/v1.1.tar.gz'
        self.serogroup_dir = 'PneumoCaT-1.1/streptococcus-pneumoniae-ctvdb/'
        self.out_file = os.path.join(os.path.abspath(out_dir),'meta.tsv')



    def _get_pneumocat_db(self):

        os.system('wget '+self.pneumocat_release_link)
        t = os.listdir(self.tmpdir)[0]

        with tarfile.open('v1.1.tar.gz') as tar:
            subdir_and_files = [
            tarinfo for tarinfo in tar.getmembers()
            if tarinfo.name.startswith(self.serogroup_dir)
            ]
            tar.extractall()

    @staticmethod
    def _pneumocat_db_2_tsv(serogroup_dir,out_file):

        serogroup_info = os.listdir(serogroup_dir)

        info_list=[]
        serogroup_dict={}
        for subdir in serogroup_info:
                if subdir != 'reference.fasta':
                    info = []
                    serogroup = subdir.split('_')[0]
                    for serotype in subdir.split('_'):
                        serogroup_dict.update({serotype:serogroup})
                    serogroup_dict.update()
                    allele_snp=yaml.load( open( os.path.join(serogroup_dir,subdir,'mutationdb.yml'), "rb" ) )
                    if 'genes' in allele_snp:
                        for serotype in allele_snp['genes']:
                            gene = list(allele_snp['genes'][serotype].keys())
                            for i in range(len(gene)):
                                gene_inf = allele_snp['genes'][serotype][gene[i]]
                                info = [gene[i],serogroup,serotype,'0','genes',str(gene_inf)]
                                info_list.append(info)
                    if 'pseudo' in allele_snp:
                        for serotype in allele_snp['pseudo']:

                            gene = list(allele_snp['pseudo'][serotype].keys())[0]
                            gene_inf = allele_snp['pseudo'][serotype][gene][0]
                            info = [gene,serogroup,serotype,'1','pseudo',str(gene_inf)]
                            info_list.append(info)
                    if 'allele' in allele_snp:
                        for genes in allele_snp['allele']:
                            for serotype in allele_snp['allele'][genes]:

                                gene = allele_snp['allele'][genes][serotype]
                                info = [gene,serogroup,serotype,'0','allele']
                                info_list.append(info)
                    if 'snps' in allele_snp:
                        for genes in allele_snp['snps']:
                            for pos in allele_snp['snps'][genes]:
                                for serotype in  allele_snp['snps'][genes][pos]:
                                    gene = genes+'_'+serotype
                                    snp = allele_snp['snps'][genes][pos][serotype][0]
                                    info = [gene,serogroup,serotype,'0','snps',str(pos),snp]
                                    info_list.append(info)
                    with open(os.path.join(serogroup_dir,subdir,'reference.fasta'),'r') as fobj:
                       with open(os.path.join(serogroup_dir, serogroup+'.fasta'),'w') as wobj:
                           for line in fobj:
                              wobj.write(line.strip()+'\n')

        reference = os.path.join(serogroup_dir,'reference.fasta')
        record_dict = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
        for seqid in record_dict:
            if seqid  in serogroup_dict:
                info = [seqid,serogroup_dict[seqid],seqid,'1','ref']
                info_list.append(info)
            else:
                info = [seqid,seqid,seqid,'1','ref']
                info_list.append(info)
        info_list.sort()
        with open(out_file,'w') as fobj:
            for entry in info_list:

                fobj.write('\t'.join(entry)+'\n')



    def run(self):
        wd  = os.getcwd()
        os.makedirs(self.out_dir)
        self.tmpdir = tempfile.mkdtemp(prefix ='pneumo.temp', dir=wd)
        os.chdir(self.tmpdir)
        self._get_pneumocat_db()
        self._pneumocat_db_2_tsv(self.serogroup_dir,self.out_file)
        os.chdir(wd)
        pneumo_db = os.path.join(self.tmpdir,self.serogroup_dir)
        os.rename(pneumo_db,os.path.join(self.out_dir,self.serogroup_dir.split('/')[1]))
        shutil.copyfile(os.path.join(self.out_dir,self.serogroup_dir.split('/')[1],'reference.fasta'), os.path.join(self.out_dir,'reference.fasta'))
        shutil.rmtree(self.tmpdir)
