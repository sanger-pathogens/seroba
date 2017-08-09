
# This container will install SeroBA from master
#
FROM debian:testing

# Install the dependancies
RUN  apt-get update -qq && apt-get install -y git bowtie2 cd-hit fastaq libc6 libfml0 libgcc1 libminimap0 libstdc++6 mummer python3 python3-setuptools python3-dev python3-pysam python3-pymummer     python3-dendropy gcc g++ zlib1g-dev
RUN apt-get -y install wget make unzip python3-pip python-pip python3-tk build-essential libbz2-dev liblzma-dev libfreetype6 libfreetype6-dev libpng-dev libxft-dev
# Get the latest code from github and install
RUN git clone https://github.com/sanger-pathogens/ariba.git && cd ariba && python3 setup.py install
RUN git clone https://github.com/eppinglen/seroba
RUN cd seroba  && ./install_dependencies.sh
env PATH /seroba/build:$PATH
RUN export PATH
RUN cd seroba && python3 setup.py install
RUN cd seroba && seroba createDBs database/ 71
