# This container will install SeroBA from master
#
FROM debian:testing

# Install the dependancies
RUN apt-get update -qq && apt-get install -y ariba python3-pip

# Get the latest code from github and install
RUN git clone https://github.com/sanger-pathogens/ariba.git && cd ariba && python3 setup.py install
RUN git clone https://github.com/sanger-pathogens/seroba
RUN cd seroba  && ./install_dependencies.sh
ENV PATH /seroba/build:$PATH
RUN export PATH
RUN cd seroba && python3 setup.py install
RUN cd seroba && seroba createDBs database/ 71
