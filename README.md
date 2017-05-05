#SeroBA
SeroBA is a kmer based Pipeline to identify the Serotype from Illumina NGS reads for given references. You can use SeroBA to download references from (https://github.com/phe-bioinformatics/PneumoCaT) to do identify the capsular type of S.pneumoniae.
##Usage
usage: seroba  getPneumocat out_dir

Downlaods PneumoCat and build an tsv formated meta data file out of it

positional arguments:
  out_dir      directory to store the PneumoCaTs capsular type variant (CTV) database


usage: seroba createDBs  <database dir> <meta data tsv>

Creates a Database for kmc and ariba

positional arguments:
    out_dir     output directory for kmc and ariba Database
    kmer_size   kmer_size zou want to use for kmc , recommanded = 51

usage: seroba runSerotyping  <databases dict> <read1> <read2> <prefix>

indetify serotype of your input data

positional arguments:
    databases   path to database directory used for seroba createDBs
    read1       forward read file
    read2       backward read file
    prefix      unique prefix used for output directory
    
##Database
You can use the CTV von PneumoCat by using seroba  getPneumocat. It is also
possible so add new serotypes by adding the references sequence to the
"references.fasta" file in the database folder. Out of  the information provided
 by this database a TSV file is created while using seroba createDBs. You can
 easily put in additional genetic information for any of these serotypes in the
 given format.

##Installation
SeroBA has the following dependencies, which need to be installed:
  * [Python3][python] version >= 3.3.2
  * [KMC][kmc] version >= 3.0
  * [KMC_tools][kmc_tools] version >= 3.0
  * [MUMmer][mummer] version >= 3.1

If this dependencies are installed, you can download the latest release from this github repository,or clone the repository.
Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install
