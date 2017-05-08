#!/bin/bash
set -e
set -x

start_dir=$(pwd)

KMC_VERSION=3.0.0
MUMMER_VERSION=3.23

KMC3_DOWNLOAD_URL="https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC3.linux.tar.gz"
MUMMER_DOWNLOAD_URL="http://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz"


# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}



# --------------- KMC ------------------
cd $build_dir

download "${KMC3_DOWNLOAD_URL}" "KMC3.linux.tar.gz"
tar xzf KMC3.linux.tar.gz
chmod +x kmc
chmod +x kmc_tools
chmod +x kmc_dump

# --------------- mummer ------------------
cd $build_dir
download $MUMMER_DOWNLOAD_URL "MUMmer${MUMMER_VERSION}.tar.gz"
mummer_dir="$build_dir/MUMmer${MUMMER_VERSION}"
tar -zxf MUMmer${MUMMER_VERSION}.tar.gz
cd $mummer_dir
make


cd $start_dir

update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${build_dir}
update_path ${mummer_dir}

pip install ariba biopython pyyaml pynummer
