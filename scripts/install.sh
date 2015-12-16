# UCSC-Kent
mkdir -p ext/kent/bin
cd ext/kent/bin
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftUp"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/axtChain"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainNet"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainSort"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainSplit"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainMergeSort"
wget -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/netChainSubset"
chmod +x *
cd -

# GNU parallel
cd ext
wget -c http://ftp.gnu.org/gnu/parallel/parallel-20150722.tar.bz2
tar xvf parallel-20150722.tar.bz2
rm parallel-20150722.tar.bz2
cd parallel-20150722
./configure
make
cd ../..

# Genometools
cd ext
wget -c https://github.com/genometools/genometools/archive/v1.5.6.tar.gz -O v1.5.6.tar.gz
tar xvf v1.5.6.tar.gz
rm v1.5.6.tar.gz
cd genometools-1.5.6
make cairo=no errorcheck=no
