from seroba import ref_db_creator, kmc
import sys
from collections import Counter
import csv
import yaml
import os
import pymummer
import xml.etree.ElementTree as ET
import tempfile
import gzip
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import copy

class Error (Exception): pass

class Serotyping:
    def __init__(self,databases, fw_reads, bw_reads, prefix,clean=True,cov=20):

        self.pneumcat_refs = os.path.join(databases,'streptococcus-pneumoniae-ctvdb')
        self.cd_cluster =  os.path.join(databases,'cd_cluster.tsv')
        self.fw_read = fw_reads
        self.bw_read = bw_reads
        self.meta_data = os.path.join(databases,'meta.tsv')
        self.prefix = prefix
        self.kmer_size = open(os.path.join(databases,'kmer_size.txt'),'r').readline().strip()
        self.kmer_db = os.path.join(databases,'kmer_db')
        self.ariba_cluster_db = os.path.join(databases,'ariba_db')
        self.reference_fasta = os.path.join(databases,'reference.fasta')
        self.meta_data_dict = ref_db_creator.RefDbCreator._read_meta_data_tsv(self.meta_data)
        self.clean = clean
        self.cov = cov/100.0

    @staticmethod
    def _serotype_2_cluster(cd_cluster):
        with open (cd_cluster,'r') as fobj:
            tsvin = csv.reader(fobj, delimiter='\t')
            serotype_cluster_dict = {'NT':'NT'}
            cluster_count = {'NT':1}
            cluster_serotype_dict = {'NT' : ['NT']}
            for row in tsvin:
                cluster_count[row[0]] = len(row)-1
                cluster_serotype_dict[row[0]] = []
                for i in range(1 , len(row)):
                    serotype_cluster_dict[row[i]] = row[0]
                    cluster_serotype_dict[row[0]].append(row[i])

            return serotype_cluster_dict, cluster_serotype_dict ,cluster_count

    def _run_kmc(self):
    #kmc on fw_read
        temp_dir = tempfile.mkdtemp(prefix = 'temp.kmc', dir=os.getcwd())
        kmer_db_list = os.listdir(self.kmer_db)
        kmer_count = self.cov
        max_kmer_count = 0.0
        best_serotype = ''
        kmer_count_files = kmc.run_kmc(self.fw_read,self.kmer_size,temp_dir,self.prefix)
        record_dict = SeqIO.to_dict(SeqIO.parse(self.reference_fasta, "fasta"))
        for db in kmer_db_list:
            #kmc_tools intersect
            curr_db = os.path.join(self.kmer_db,db,db)
            temp_inter = kmc.run_kmc_intersect(kmer_count_files,temp_dir,curr_db)
            #kmc_tools get unique kmer count transform inter histogram
            temp_hist = kmc.run_kmc_hist(temp_inter,temp_dir)
            with open( temp_hist, 'r') as fobj:
                first_line = fobj.readline()
                unique_kmers = int(first_line.split('\t')[1])/float(len(record_dict[db].seq))
                print(unique_kmers)
                if unique_kmers > kmer_count:
                    kmer_count = unique_kmers
                    best_serotype = db
                if unique_kmers > max_kmer_count:
                   max_kmer_count = unique_kmers
        shutil.rmtree(temp_dir)
        if max_kmer_count < 0.01:
           self.best_serotype = 'coverage too low'
        elif best_serotype != '':
            self.best_serotype = best_serotype
        else:
            self.best_serotype = 'NT'

    def _run_ariba_on_cluster(self,cluster):
        os.makedirs(self.prefix)
        ref_dir = os.path.join(self.ariba_cluster_db ,self.cluster_serotype_dict[cluster][0]+'/','ref')
        command = ['ariba run ',ref_dir,self.fw_read,self.bw_read,os.path.join(self.prefix,'ref')]
        os.system(' '.join(command))
        if (os.path.isdir(os.path.join(self.ariba_cluster_db ,self.cluster_serotype_dict[cluster][0]+'/','genes'))):
           ref_dir = os.path.join(self.ariba_cluster_db ,self.cluster_serotype_dict[cluster][0]+'/','genes')
           command = ['ariba run ',ref_dir,self.fw_read,self.bw_read,os.path.join(self.prefix,'genes')]
           os.system(' '.join(command))
           shutil.copyfile(os.path.join(self.prefix,'genes','assembled_genes.fa.gz'),os.path.join(self.prefix,'assembled_genes.fa.gz'))
           os.system('cat '+ os.path.join(self.prefix,'ref','report.tsv')+' '+os.path.join(self.prefix,'genes','report.tsv')+' > ' + os.path.join(self.prefix,'report.tsv'))
           os.system('gzip -d '+os.path.join(self.prefix,'assembled_genes.fa.gz'))
           shutil.copyfile(os.path.join(self.prefix,'genes','assembled_genes.fa.gz'),os.path.join(self.prefix,'assembled_genes.fa.gz'))

        else:
            shutil.copyfile(os.path.join(self.prefix,'ref','report.tsv'),os.path.join(self.prefix,'report.tsv'))
        shutil.copyfile(os.path.join(self.prefix,'ref','assemblies.fa.gz'),os.path.join(self.prefix,'assemblies.fa.gz'))
        os.system('gzip -d '+os.path.join(self.prefix,'assemblies.fa.gz'))



    @staticmethod
    def serotype6(assemblie_file,report_file):
        #os.system('gzip -d '+assemblie_file)
        with open(report_file,'r') as fobj:
            next(fobj)
            first = fobj.readline().split('\t')[0]
        with open(report_file) as fobj:
            serotype = 'possible '+first+'\t but wciP gene might not be complete'
            tsvin = csv.reader(fobj, delimiter='\t')
            next(tsvin,None)
            row_dict = {}
            for row in tsvin:
                if row[0] not in row_dict:
                	row_dict[row[0]]=row[10]
            record= SeqIO.to_dict(SeqIO.parse(assemblie_file, "fasta"))
            if 'wciN_3'in row_dict:
                    for seq_id in row_dict:
                        if 'wciP' in seq_id:
                            for seq in record:
                                if row_dict[seq_id] in seq:
                                    snp = (record[seq].seq[583])
                                    if snp == 'G':
                                        serotype = '6E(6A)'
                                    elif snp =='A':
                                        serotype = '6E(6B)'

            elif 'wciN_1'in row_dict:
                for seq_id in row_dict:
                    if 'wciP' in seq_id:
                        for seq in record:
                            if row_dict[seq_id] in seq:
                                snp = (record[seq].seq[583])
                                if snp == 'G':
                                    serotype = '06A'
                                    for seq_id_2 in row_dict:
                                        if 'wciN_1' in seq_id_2:
                                            for seq_2 in record:
                                                if row_dict[seq_id_2] in seq_2:
                                                    if  (record[seq_2].seq[447]) == 'A':
                                                        serotype = '06F'
                                                    elif (record[seq_2].seq[112]) == 'A':
                                                        serotype = '06G'
                                elif snp =='A':
                                    serotype = '06B'

            elif 'wciN_2'in row_dict:
                for seq_id in row_dict:
                    if 'wciP' in seq_id:
                        for seq in record:
                            if row_dict[seq_id] in seq:
                                snp = (record[seq].seq[583])
                                if snp == 'G':
                                    serotype = '06C'
                                elif snp =='A':
                                    serotype = '06D'



        return serotype

    @staticmethod
    def _detect_mixed_samples(line,serotype_gene_dict):
        mixed_serotype = []
        if 'HET' in line[14]:
            for serotype in serotype_gene_dict:
                mixed_serotype.append(serotype)

            return '/'.join(sorted(mixed_serotype))



    @staticmethod
    def _get_present_genes(nucmer_hits,assemblie_file):
        record_dict = SeqIO.to_dict(SeqIO.parse(assemblie_file, "fasta"))
        present_genes ={}
        for seqid in record_dict:
            if seqid.replace('-','_') in nucmer_hits:
                present_genes[seqid.replace('-','_')] = '1'
            else:
                present_genes[seqid.replace('-','_')] = '0'

        return present_genes
    @staticmethod
    def _get_nucmer_snps(variants,genes):
        variant_dict = dict.fromkeys(genes)
        for v in variants:
            if variant_dict[v.ref_name.replace('-','_')]==None:
                variant_dict[v.ref_name.replace('-','_')]={v.ref_start: v.ref_base}
            else:
                variant_dict[v.ref_name.replace('-','_')].update({v.ref_start : v.ref_base})

        return variant_dict
    @staticmethod
    def _get_snp(snp_dict,serotype):

        triplets=[snp_dict[serotype]]
        bases = []
        for se in snp_dict:
            triplets.append(snp_dict[se])
        for i in range(3):
            base=[tripl[i] for tripl in triplets]
            if all(x == base[0] for x in base)== False :
               bases.append(base[0])

        return bases
    @staticmethod
    def _check_snps(yaml_dict,variant_dict,serotype_count,relevant_genetic_elements,prefix,report_file):
        record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(prefix,'assembled_genes.fa'), "fasta"))
        with open(report_file) as fobj:
            serotype = ''
            tsvin = csv.reader(fobj, delimiter='\t')
            next(tsvin,None)
            row_dict = {}
            for row in tsvin:
                if row[10] not in row_dict:
                	row_dict[row[10]]=row[0]
        record_dict2={}
        for seqid in row_dict:
            for seq in record_dict:
                if seqid in seq:
                    record_dict2[row_dict[seqid].split('_')[0]]=record_dict[seq]

        for gene in yaml_dict:
            for pos in yaml_dict[gene]:
                for serotype in serotype_count:
                    if serotype in yaml_dict[gene][pos]:
                        bases=Serotyping._get_snp(yaml_dict[gene][pos],serotype)

                        if int(pos)%3 == 0:
                            i = int(pos)
                            j = int(pos)+3
                        else:
                            i = int(pos)-1
                            j = int(pos)+2
                        if gene in record_dict2:
                            #for i in range(len(bases)):
                            if yaml_dict[gene][pos][serotype] == record_dict2[gene].seq[i:j]:
                                sub_dict = copy.deepcopy(relevant_genetic_elements[serotype])
                                sub_dict['snps'].append([gene,pos, yaml_dict[gene][pos][serotype]])
                                relevant_genetic_elements[serotype] = sub_dict
                                serotype_count[serotype]+=-1

                    #    else:
                    #        serotype_count[serotype]+=-0.5

        return serotype_count,relevant_genetic_elements

    @staticmethod
    def _find_serotype(assemblie_file,serogroup_fasta, serogroup_dict,serotypes,report_file,prefix):
        sub_dict = {'genes':[],'pseudo':[],'allele':[],'snps':[]}
        relevant_genetic_elements = dict.fromkeys(serotypes, sub_dict)
        allel_snp = serogroup_dict
        tmpdir = tempfile.mkdtemp(prefix = 'temp.nucmer', dir=os.getcwd())
        gene_ref = serogroup_fasta
        pymummer.nucmer.Runner(
            gene_ref,
            assemblie_file,
            os.path.join(tmpdir,'coords.txt'),
            min_id=90,
            min_length=200,
            maxmatch=True,
            show_snps=True,
            show_snps_C=False,
        ).run()
        hits = [[x.ref_name, x.ref_start, x.ref_end, x.ref_length,x.qry_length,x.ref_length] for x in pymummer.coords_file.reader(os.path.join(tmpdir,'coords.txt'))]
        serotype_count = dict.fromkeys(serotypes)
        variants = pymummer.snp_file.get_all_variants(os.path.join(tmpdir,'coords.txt.snps'))
        for key in serotype_count:
            serotype_count[key] = 0

        for i in range(len(hits)):
          if (int((hits[i][2])-int(hits[i][1]))/ float(hits[i][3])) < 0.5:
             hits[i]='0'
          elif int(hits[i][4])< 2000:
            hits[i]='0'
        for i in range(len(hits)):
           hits[i]=hits[i][0].replace('-','_')
        variants = pymummer.snp_file.get_all_variants(os.path.join(tmpdir,'coords.txt.snps'))
        hits = [ x for x in hits if '0' != x ]
        #check genes (present or not)
        gene_present_dict = Serotyping._get_present_genes(hits,gene_ref)
        for gene in gene_present_dict:
            for serotype in serotypes:
                if serotype in gene and gene_present_dict[gene] == 1:
                    serotype_count[serotype]+=-1

        variant_dict = Serotyping._get_nucmer_snps(variants, list(gene_present_dict.keys()) )
        mixed_serotype = None
        if "genes" in allel_snp:
            for serotype in serotypes:
                for gene in gene_present_dict:
                    if gene in allel_snp['genes'][serotype]:
                        if allel_snp['genes'][serotype][gene] != gene_present_dict[gene]:
                            serotype_count[serotype]+=4
                        if  gene_present_dict[gene] == '1':
                            sub_dict = copy.deepcopy(relevant_genetic_elements[serotype])
                            sub_dict['genes'].append(gene)
                            relevant_genetic_elements[serotype]=sub_dict


        print(serotype_count)
        if "pseudo" in allel_snp:

            for serotype in serotypes:
                sub_dict = copy.deepcopy(relevant_genetic_elements[serotype])
                print(serotype)
                print(sub_dict)

                if serotype in allel_snp['pseudo']:
                   for gene in allel_snp['pseudo'][serotype]:
                        if allel_snp['pseudo'][serotype][gene]=='1':
                            with open(report_file) as fobj:
                                tsvin = csv.reader(fobj, delimiter='\t')
                                for row in tsvin:
                                    if gene in row:
                                        mixed_serotype = Serotyping._detect_mixed_samples(row,allel_snp['pseudo'])
                                    if gene in row and ("FSHIFT" in row or 'TRUNC' in row):
                                        serotype_count[serotype]+=-2.5
                                        print(serotype)

                                        sub_dict['pseudo'].append(gene)



                        elif allel_snp['pseudo'][serotype][gene] == '0':
                            with open(report_file) as fobj:
                                tsvin = csv.reader(fobj, delimiter='\t')
                                next(tsvin,None)
                                count = 0
                                for row in tsvin:

                                       if gene in row and (float(row[8])/float(row[7]) >0.95):
                                           mixed_serotype = Serotyping._detect_mixed_samples(row,allel_snp['pseudo'])
                                       if gene in row and ("FSHIFT" in row or 'TRUNC' in row) and (float(row[8])/float(row[7]) > 0.95):
                                           count = 1

                                if count == 0:
                                   serotype_count[serotype] +=-2.5
                                   print(gene)
                                   sub_dict['pseudo'].append(gene)
                                   print(sub_dict)

                        relevant_genetic_elements[serotype] = sub_dict

                else:
                    serotype_count[serotype]+=-1
        print(serotype_count)

        if "allele" in allel_snp:
            for al in allel_snp['allele']:
                h = [[x.ref_name, x.ref_start, x.ref_end, x.ref_length,x.qry_length,x.ref_length,x.percent_identity] for x in pymummer.coords_file.reader(os.path.join(tmpdir,'coords.txt'))]
                for i in range(len(h)):
                     if (int((h[i][2])-int(h[i][1]))/ float(h[i][3])) < 0.5:
                        h[i][0]='0'
                     elif int(h[i][4])< 2000:
                       h[i][0]='0'
                best_al = ''
                score =0
                s =['']
                for serotype in allel_snp['allele'][al]:
                    for i in range(len(h)):
                        if h[i][0].replace('-','_') in allel_snp['allele'][al][serotype].replace('-','_') and float(h[i][6]) > score:
                           best_al = h[i][0].replace('-','_')
                           score = float(h[i][6])
                           s[0]= serotype
                        elif best_al == allel_snp['allele'][al][serotype].replace('-','_') and float(h[i][6]) == score and serotype not in s:
                              s.append(serotype)
                              sub_dict = copy.deepcopy( relevant_genetic_elements[serotype])
                              sub_dict['allele'].append(al)
                              relevant_genetic_elements[serotype]=sub_dict

                for se in s :
                   if se in serotype_count:
                       serotype_count[se]+=-1.5

        if "snps" in allel_snp:
            #check snps
            serotype_count,relevant_genetic_elements = Serotyping._check_snps(allel_snp['snps'],variant_dict,serotype_count,relevant_genetic_elements,prefix,
            report_file)

        print(relevant_genetic_elements)
        shutil.rmtree(tmpdir)
        min_value = min(serotype_count.values())
        min_keys = [k for k in serotype_count if serotype_count[k] == min_value]
        serotype = ''
        print(min_keys)
        if mixed_serotype != None:
            for key in min_keys:
                print(key)
                print(mixed_serotype)
                if key not in mixed_serotype:
                   mixed_serotype = None
        print(serotype_count)
        if  mixed_serotype!= None :
            serotype = mixed_serotype
        elif len(min_keys) > 1:
            with open(report_file) as fobj:
                tsvin = csv.reader(fobj, delimiter='\t')
                next(tsvin,None)
                first = next(tsvin)
                serotype = first[0]

        elif  min(serotype_count, key=serotype_count.get) =='33A':
                with open(report_file) as fobj:
                    tsvin = csv.reader(fobj, delimiter='\t')
                    next(tsvin,None)
                    first = next(tsvin)
                    if first[0] == min(serotype_count, key=serotype_count.get):
                        serotype = min(serotype_count, key=serotype_count.get)
                    else:
                        serotype = first[0]
        else :
            serotype =  min(serotype_count, key=serotype_count.get)
        return serotype , relevant_genetic_elements


    def _print_detailed_output(self,report_file,relevant_genetic_elements,serotype):
        print(relevant_genetic_elements)
        with open(report_file,'r') as fobj:
            next(fobj)
            first = fobj.readline().split('\t')
            with open(os.path.join(self.prefix,'detailed_serogroup_info.txt'),'w') as wobj:
                wobj.write('Predicted Serotype:\t'+ serotype+'\n')
                wobj.write('Serotype predicted by ariba\t:' +first[0]+'\n')
                wobj.write('assembly from ariba as an identiy of: \t'+ first[9]+'\t with this serotype\n')
                wobj.write('Serotype \t Genetic Variant\n')
                for serotype in relevant_genetic_elements:
                    for genetic_var in relevant_genetic_elements[serotype]:
                        for entry in relevant_genetic_elements[serotype][genetic_var]:
                            wobj.write(serotype+'\t'+genetic_var+'\t'+str(entry)+'\n')

    def _prediction(self,assemblie_file,cluster):
        sero = ''

        #db_path = os.path.join(self.pneumcat_refs,'_'.join(sorted(self.cluster_serotype_dict[cluster])))
        if self.cluster_count[cluster] == 1:
            self.sero = self.cluster_serotype_dict[cluster][0]
        elif '06' in self.best_serotype:
            report_file  = os.path.join(self.prefix,'report.tsv')
            assemblie_file = os.path.join(self.prefix,'assembled_genes.fa')
            self.sero = Serotyping.serotype6(assemblie_file, report_file)
        else:
            report_file = os.path.join(self.prefix,'report.tsv')
            serogroup = self.cluster_serotype_dict[cluster][0]
            serogroup_fasta = os.path.join(self.pneumcat_refs,serogroup+'.fasta')
            self.sero, self.imp = Serotyping._find_serotype(assemblie_file,serogroup_fasta,self.meta_data_dict[serogroup],\
                self.cluster_serotype_dict[cluster],report_file,self.prefix)
            self._print_detailed_output(report_file,self.imp,self.sero)

    def _clean(self):
        files = os.listdir(self.prefix)
        for f in files:
            if 'pred.tsv' != f and 'detailed_serogroup_info.txt' != f :
                path = os.path.join(self.prefix,f)
                os.remove(path)


    def run(self):
        self.serotype_cluster_dict, self.cluster_serotype_dict,\
        self.cluster_count = Serotyping._serotype_2_cluster(self.cd_cluster)
        assemblie_file = self.prefix+'/assemblies.fa'
        self._run_kmc()
        print(self.best_serotype)
        if self.best_serotype =='coverage too low':
           os.system('mkdir '+self.prefix)
           with open(self.prefix+'/pred.tsv', 'a') as fobj:
               fobj.write(self.prefix+'\t'+self.best_serotype+'\n')
        elif self.best_serotype == 'NT':
            os.system('mkdir '+self.prefix)
            with open(self.prefix+'/pred.tsv', 'a') as fobj:
                fobj.write(self.prefix+'\tuntypable\n')
        else:
            cluster = self.serotype_cluster_dict[self.best_serotype]
            self._run_ariba_on_cluster(cluster)
            self._prediction(assemblie_file,cluster)
            report_file = os.path.join(self.prefix,'report.tsv')
            flag = ''
            with open (report_file,'r') as report:
                for line in report:
                    if 'HET' in line:
                        flag = 'contamination'
            with open(self.prefix+'/pred.tsv', 'a') as fobj:
                if '24' in self.sero:
                    fobj.write(self.prefix+'\tserogroup 24\t'+flag+'\n')
                else:
                    fobj.write(self.prefix+'\t'+self.sero+'\t'+flag+'\n')
            shutil.rmtree(os.path.join(self.prefix,'ref'))
        if os.path.isdir(os.path.join(self.prefix,'genes')):
            shutil.rmtree(os.path.join(self.prefix,'genes'))
        if self.clean:
            self._clean()
