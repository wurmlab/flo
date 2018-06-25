# UCSC-Kent
mkdir -p ext/kent/bin; cd ext/kent/bin
tools=( liftUp faSplit liftOver axtChain chainNet blat chainSort faToTwoBit
twoBitInfo chainSplit chainMergeSort netChainSubset )
case "$(uname -s)" in
  Darwin*) ftp_dir="http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64";;
  Linux*) ftp_dir="http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64";;
esac
for tool in ${tools}; do wget -c "${ftp_dir}/${tool}"; done
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
