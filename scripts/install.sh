#!/usr/bin/env bash

# UCSC-Kent
if [[ "$*" == *"ucsc-kent"* ]] ## Check all the argument to the install.sh script
then
mkdir -p ext/kent/bin; cd ext/kent/bin
tools=(liftUp faSplit liftOver axtChain chainNet blat/blat chainSort faToTwoBit
twoBitInfo chainSplit chainMergeSort netChainSubset)
case "$(uname -s)" in
  Darwin*) ftp_dir="http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64";;
  Linux*) ftp_dir="http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64";;
esac
for tool in ${tools[@]}; do wget -c "${ftp_dir}/${tool}"; done
chmod +x *
cd -
fi

# GNU parallel
if [[ "$*" == *"parallel"* ]] ## Check all the argument to the install.sh script
then
cd ext
wget -c http://ftp.gnu.org/gnu/parallel/parallel-20150722.tar.bz2
tar xvf parallel-20150722.tar.bz2
rm parallel-20150722.tar.bz2
cd parallel-20150722
./configure
make
cd ../..
fi

# Genometools
if [[ "$*" == *"genometools"* ]] ## Check all the argument to the install.sh script
then
cd ext
wget -c https://github.com/genometools/genometools/archive/refs/tags/v1.6.2.tar.gz -O v1.6.2.tar.gz
tar xvf v1.6.2.tar.gz
rm v1.6.2.tar.gz
cd genometools-1.6.2
make cairo=no errorcheck=no
fi
