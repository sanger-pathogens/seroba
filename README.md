# SeroBA
SeroBA is a k-mer based Pipeline to identify the Serotype from Illumina NGS reads for given references. You can use SeroBA to download references from (https://github.com/phe-bioinformatics/PneumoCaT) to do identify the capsular type of Streptococcus pneumoniae.

[![Build Status](https://travis-ci.org/sanger-pathogens/seroba.svg?branch=master)](https://travis-ci.org/sanger-pathogens/seroba)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/seroba/blob/master/LICENSE)   
[![status](https://img.shields.io/badge/MGEN-10.1099%2Fmgen.0.000056-brightgreen.svg)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000186)   
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/seroba/README.html)  
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/seroba)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/seroba.svg)](https://hub.docker.com/r/sangerpathogens/seroba)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/seroba.svg)](https://hub.docker.com/r/sangerpathogens/seroba)  
[![codecov](https://codecov.io/gh/sanger-pathogens/seroba/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/seroba)

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [conda](#conda)
    * [CentOS 7](#centos-7)
    * [Debian Testing/ Ubuntu 17\.10](#debian-testing-ubuntu-1710)
    * [Ubuntu 16\.04 (Xenial)](#ubuntu-1604-xenial)
    * [Docker](#docker)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
    * [Setting up the database](#setting-up-the-database)
    * [Usage](#usage-1)
    * [Output](#output)
  * [Troubleshooting](#troubleshooting)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)
  * [Further Information](#further-information)
    * [Tutorial](#tutorial)

## Introduction
SeroBA can predict serotypes, by identifying the cps locus, directly from raw whole genome sequencing read data with 98% concordance using a k-mer based method, can process 10,000 samples in just over 1 day using a standard server and can call serotypes at a coverage as low as 10x. SeroBA is implemented in Python3 and is freely available under an open source GPLv3

## Installation
SeroBA has the following dependencies:

### Required dependencies
* Python3 version >= 3.3.2
* KMC version >= 3.0
* MUMmer version >= 3.1
* Ariba

There are a number of ways to install SeroBA and details are provided below. If you encounter an issue when installing SeroBA please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/seroba/issues).

### conda
Set up bioconda channel:
```
conda config --add channels bioconda
```
Install SeroBA:
```
conda install -c bioconda seroba
```

### CentOS 7
Ensure you have a development environment setup (you may have done this already):
```
yum -y update
yum -y groupinstall 'Development Tools'
yum -y install https://centos7.iuscommunity.org/ius-release.rpm
```

Install SeroBA and its dependancies:
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

### Docker
Install [Docker](https://www.docker.com/).  We have a docker container which gets automatically built from the latest version of SeroBA. To install it:
```
docker pull sangerpathogens/seroba
```
To use it you would use a command such as this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:
```
docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/seroba seroba runSerotyping seroba/database /data/read_1.fastq.gz /data/read_2.fastq.gz  /data/output_folder
```    

### Running the tests
The test can be run from the top level directory:  

```
python setup.py test
```

## Usage
### Setting up the database
You can use the CTV of PneumoCaT by using seroba  getPneumocat. It is also possible to add new serotypes by adding the references sequence to the "references.fasta" file in the database folder. Out of the information provided by this database a TSV file is created while using seroba createDBs. You can easily put in additional genetic information for any of these serotypes in the given format.

Since SeroBA v0.1.3 an updated variant of the CTV from PneumoCaT is provided in the SeroBA package. This includes the serotypes 6E, 6F, 11E, 10X, 39X and two NT references. It is not necessary to use SeroBA getPneumocat.
For SeroBA version 0.1.3 and greater, download the database provided within this git repository:

__For git users__
```
Clone the git repository:
git clone https://github.com/sanger-pathogens/seroba.git
```

Copy the database to a directory:
```
cp -r seroba/database my_directory
```

Delete the git repository to clean up your system:
```
rm -r seroba
```

__For svn users__  
Install svn. Checkout the database directory:
```
svn checkout "https://github.com/sanger-pathogens/seroba/trunk/database"
```
Continue with Step 2.

__For SeroBA version 0.1.2 and smaller:__
```
usage: seroba  getPneumocat <database dir>

Downloads PneumoCat and build an tsv formatted meta data file out of it

positional arguments:
  database dir      directory to store the PneumoCats capsular type variant (CTV) database
```

### Usage
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

### Output
In the folder 'prefix' you will find a pred.tsv including your predicted serotype as well as a file called detailed_serogroup_info.txt including information about SNP, genes, and alleles that are found in your reads. After the use of "seroba summary" a tsv file called summary.tsv is created that consists of three columns (sample Id , serotype, comments). Serotypes that do not match any reference are marked as "untypable"(v0.1.3).

__detailed_serogroup_info example:__
```
Predicted Serotype:       23F
Serotype predicted by ariba:    23F
assembly from ariba has an identity of:   99.77    with this serotype

Serotype       Genetic Variant
23F            allele  wchA
```
In the detailed information you can see the finally predicted serotype as well as the serotypes that had the closest reference in that specific serogroup according to ARIBA. Furthermore you can see the sequence identity between the sequence assembly and the reference sequence.  

## Troubleshooting
* Case 1:
	* SeroBA predicts 'untypable'. An 'untypable' prediction can either be a
real 'untypable' strain or can be caused by different problems. Possible problems are:
bad quality of your input data, submission of a wrong species or to low coverage
of your sequenced reads. Please check your data again and run a quality control.

* Case 2:
	* 	Low alignment identity in the 'detailed_serogroup_info' file. This can
be a hint for a mosaic serotpye.
	* Possible solution: perform a blast search on the whole genome assembly

* Case 3:
	* The third column in the summary.tsv indicates "contamination". This means that
    at least one heterozygous SNP was detected in the read data with at least
    10% of the mapped reads at the specific position supporting the SNP.
	* Possible solution: please check the quality of your data and have a look
    for contamination within your reads

## License
SeroBA is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/seroba/blob/master/LICENSE)

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/seroba/issues).

## Citation
__SeroBA: rapid high-throughput serotyping of Streptococcus pneumoniae from whole genome sequence data__  
Epping L, van Tonder, AJ, Gladstone RA, GPS Consortium, Bentley SD, Page AJ, Keane JA, Microbial Genomics 2018, doi: [10.1099/mgen.0.000186](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000186)

## Further Information
### Tutorial
A tutorial for SeroBA can be found here:

https://github.com/sanger-pathogens/pathogen-informatics-training/
