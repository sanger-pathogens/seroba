# SeroBA
SeroBA is a k-mer based Pipeline to identify the Serotype from Illumina NGS reads for given references. You can use SeroBA to download references from (https://github.com/phe-bioinformatics/PneumoCaT) to do identify the capsular type of Streptococcus pneumoniae.

## Tutorial
A tutorial for SeroBA can be found here:

https://github.com/sanger-pathogens/pathogen-informatics-training

## Usage
Since SeroBA v0.1.3 an updated variant of the CTV from PneumoCat is provided in the SeroBA package. This includes the serotypes 6E, 6F, 11E, 10X, 39X and two NT references. It is not necessary to use SeroBA getPneumocat.

For SeroBA version 0.1.3 and greater, download the database provided within this git repository:
  * Install svn
```
svn checkout "https://github.com/sanger-pathogens/seroba/trunk/database"
```
Continue with Step 2.

############################

1.
For SeroBA version 0.1.2 and smaller:

```
usage: seroba  getPneumocat <database dir>

Downloads PneumoCat and build an tsv formatted meta data file out of it

positional arguments:
  database dir      directory to store the PneumoCats capsular type variant (CTV) database
```
###########################

2.
```

usage: seroba createDBs  <database dir> <kmer size>

Creates a Database for kmc and ariba

positional arguments:
    database dir     output directory for kmc and ariba Database
    kmer size   kmer_size you want to use for kmc , recommended = 71

    usage: seroba runSerotyping [options]  <databases directory> <read1> <read2> <prefix>

    Example : seroba createDBs my_database/ 71

Identify serotype of your input data

    positional arguments:
      database dir         path to database directory
      read1              forward read file
      read2              reverse read file
      prefix             unique prefix

    optional arguments:
      -h, --help         show this help message and exit

    Other options:
      --noclean NOCLEAN  Do not clean up intermediate files (assemblies, ariba
                         report)
      --coverage COVERAGE  threshold for k-mer coverage of the reference sequence (default = 20)                         



Summaries the output in one tsv file

usage: seroba summary  <output folder>

positional arguments:
  output folder   directory where the output directories from seroba runSerotyping are stored


```   

## Output
In the folder 'prefix' you will find a pred.tsv including your predicted serotype
as well as en file called detailed_serogroup_info.txt including information about
snps, genes, and alleles that are found in your reads.
After the use of "seroba summary" a tsv file called summary.tsv is created that
consists of three columns (sample Id , serotype, comments).
Serotypes that do not match any reference are marked as "untypable"(v0.1.3).

## Database
You can use the CTV von PneumoCat by using seroba  getPneumocat. It is also
possible so add new serotypes by adding the references sequence to the
"references.fasta" file in the database folder. Out of  the information provided
 by this database a TSV file is created while using seroba createDBs. You can
 easily put in additional genetic information for any of these serotypes in the
 given format.

## Installation

### CentOS 7
Ensure you have a development environment setup (you may have done this already):
```
yum -y update
yum -y groupinstall 'Development Tools'
yum -y install https://centos7.iuscommunity.org/ius-release.rpm
```

Install seroba and its dependancies:
```
yum -y install python36u python36u-pip python36u-devel zlib-devel wget which python36u-tkinter
ln -s $(which pip3.6) /usr/bin/pip3
bash <(curl -fsSL https://raw.githubusercontent.com/sanger-pathogens/seroba/master/install_dependencies.sh)
```
Make sure to add the PATHs outputted by this script to your .bashrc file (or equivalent). Finally install SeroBA:
```
pip3 install seroba
```

### Debian Testing/ Ubuntu 17.10

Install the dependancies:
```
sudo apt-get update
sudo apt-get install ariba python3-pip wget
```

Manually install [KMC version 3](https://github.com/refresh-bio/KMC/releases) (version 2 is the latest in Debian but is incompatible).
Add the binaries to your PATH (e.g. in your bash profile).
```
mkdir kmc && cd kmc
wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz
tar xvfz KMC3.linux.tar.gz
export PATH=$PWD:$PATH
```

Finally install SeroBA:
```
pip3 install seroba
```

### Ubuntu 16.04 (Xenial)


SeroBA has the following dependencies, which need to be installed:
  * Python3 version >= 3.3.2
  * KMC version >= 3.0
  * MUMmer version >= 3.1

Install the dependancies:
```
apt-get update
apt-get install --no-install-recommends -y build-essential cd-hit curl git libbz2-dev liblzma-dev mummer python python3-dev python3-setuptools python3-pip python3-tk python3-matplotlib unzip wget zlib1g-dev
wget -q https://raw.githubusercontent.com/sanger-pathogens/seroba/master/install_dependencies.sh && bash ./install_dependencies.sh
```

Once the dependencies are installed, install SeroBA using pip:
```
pip3 install seroba
```

# Linux/OSX/Windows/Cloud
## Docker
Install [Docker](https://www.docker.com/).  We have a docker container which gets automatically built from the latest version of SeroBA. To install it:

```
docker pull sangerpathogens/seroba
```
To use it you would use a command such as this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:
```
docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/seroba seroba runSerotyping seroba/database /data/read_1.fastq.gz /data/read_2.fastq.gz  /data/output_folder
```    
