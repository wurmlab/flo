#!/usr/bin/env ruby

# Copyright 2016 Anurag Priyam - MIT License
#
# Extract exonic, coding and protein sequences from the given FASTA file,
# corresponding to the annotations given GFF3 file.

require_relative 'helpers'

# I/O: 1 FASTA, 1 GFF3 -> 3 FASTA (exonic, coding, protein).
fasta, gff3, exonic_fa, coding_fa, protein_fa = ARGV

# Get handle to a sorted and processed, temporary GFF3 file.
gt_gff3 = process_with_gt(gff3, gt_gff3)

unless exonic_fa.nil? || exonic_fa.empty?
  sh 'gt extractfeat -type exon -join -retainids -coords'                      \
     " -matchdescstart -seqfile #{fasta} #{gt_gff3.path} > #{exonic_fa}"
end

unless coding_fa.nil? || coding_fa.empty?
  sh 'gt extractfeat -type CDS -join -retainids -coords'                       \
     " -matchdescstart -seqfile #{fasta} #{gt_gff3.path} > #{coding_fa}"
end

unless protein_fa.nil? || protein_fa.empty?
  sh 'gt extractfeat -type CDS -join -translate -retainids -coords'            \
     " -matchdescstart -seqfile #{fasta} #{gt_gff3.path} > #{protein_fa}"
end

sh "rm #{fasta}.*"
