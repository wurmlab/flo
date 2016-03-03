# Copyright 2016 Anurag Priyam - MIT License
#
# Downloads NCBI's fire ant annotation on Si_gnG, process data and lift over
# annotations to same assembly: should lift over 100%.

set -ex

### Download and process data. ###

mkdir data && cd data

# Download genome. Merge chromosomal and mitochondrial genome while trimming
# long form NCBI sequence ids to accession numbers into Si_gnG.fa.
wget -c ftp://ftp.ncbi.nih.gov/genomes/Solenopsis_invicta/CHR_Un/13686_ref_Si_gnG_chrUn.fa.gz
gunzip 13686_ref_Si_gnG_chrUn.fa.gz
wget -c ftp://ftp.ncbi.nih.gov/genomes/Solenopsis_invicta/CHR_MT/13686_ref_Si_gnG_chrMT.fa.gz
gunzip 13686_ref_Si_gnG_chrMT.fa.gz
cat 13686_ref_Si_gnG_chrUn.fa 13686_ref_Si_gnG_chrMT.fa | ruby -pe 'gsub(/^>gi\|\d+\|\w+\|(.+)\| /) {">#{$1} "}' > Si_gnG.fa

# Download annotations. Alias as NCBIv100_Si_gnG.gff3.
wget -c ftp://ftp.ncbi.nih.gov/genomes/Solenopsis_invicta/GFF/ref_Si_gnG_scaffolds.gff3.gz
gunzip ref_Si_gnG_scaffolds.gff3.gz
ln -s ref_Si_gnG_scaffolds.gff3 NCBIv100_Si_gnG.gff3

# Select all mRNA and their children.
../scripts/select_mrna.rb NCBIv100_Si_gnG.gff3 > NCBIv100_Si_gnG.mRNA.gff3

# Get exonic, coding and protein sequences for all mRNA.
../scripts/extractfeat.rb Si_gnG.fa NCBIv100_Si_gnG.mRNA.gff3 NCBIv100_Si_gnG.mRNA.gff3.exonic.fa NCBIv100_Si_gnG.mRNA.gff3.coding.fa NCBIv100_Si_gnG.mRNA.gff3.protein.fa

# Filter out bad mRNA. See filter.rb for what constitutes bad mRNA.
../scripts/filter_mrna.rb NCBIv100_Si_gnG.mRNA.gff3 NCBIv100_Si_gnG.mRNA.gff3.exonic.fa NCBIv100_Si_gnG.mRNA.gff3.coding.fa > NCBIv100_Si_gnG.mRNA.filtered.gff3

# Select mRNA with longest CDS per loci, and select longest mRNA per loci.
../scripts/select_mrna.rb NCBIv100_Si_gnG.mRNA.filtered.gff3 NCBIv100_Si_gnG.mRNA.gff3.coding.fa > NCBIv100_Si_gnG.mRNA.filtered.longest_coding.gff3
../scripts/select_mrna.rb NCBIv100_Si_gnG.mRNA.filtered.gff3 longest > NCBIv100_Si_gnG.mRNA.filtered.longest.gff3

cd -

### gnG to gnG ###
cd test/gnG_to_gnG
ln -s ../../data
ln -s ../scripts
ln -s ../../ext
rake

### gnG to gnH ###
if [ -f data/Si_gnH.fa ]; then
  cd test/gnG_to_gnH
  ln -s ../../data
  ln -s ../scripts
  ln -s ../../ext
  rake
fi
