#
# This container will install SeroBA from master
#
FROM ubuntu:latest



#
# Install the dependancies
#
RUN apt-get update -qq && apt-get install -y git libc6 libfml0 libgcc1 libminimap0 libstdc++6 mummer python3 python3-setuptools python3-dev python3-pysam python3-pymummer python3-dendropy

#
# Get the latest code from github and install
#
RUN git clone https://github.com/eppinglen/seroba && cd seroba
RUN bash install_dependencies.sh && python3 setup.py install
