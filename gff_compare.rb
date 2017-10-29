#!/usr/bin/env ruby
# Copyright 2017 Anurag Priyam, Queen Mary University of London
#
# Compares two gff files and prints id of the transcripts that
# are not present in the second gff file and that do not yield
# identical cdna, cds, or protein sequence.
#
# Is used by flo to generate a list of ids of the transcripts
# that couldn't be lifted.
#
# cdna, cds, or protein sequences obtained from both the gff
# are read into memory for comparison :TODO optimise:

require 'rake'

type, source_fa, target_fa, source_gff, target_gff = ARGV

def extract_cdna(fas, gff)
  sh "gt gff3 -sort -retainids #{gff} | gt extractfeat" \
     " -type exon -join -retainids -seqfile #{fas} "    \
     " -matchdescstart - > #{gff.ext('.cdna.fa')}"
end

def extract_cds(fas, gff)
  sh "gt gff3 -sort -retainids #{gff} | gt extractfeat" \
     " -type CDS -join -retainids -seqfile #{fas}"      \
     " -matchdescstart - > #{gff.ext('.cds.fa')}"
end

def extract_pep(fas, gff)
  sh "gt gff3 -sort -retainids #{gff} | gt extractfeat" \
     " -type CDS -join -translate -retainids"           \
     " -seqfile #{fas} -matchdescstart - >"             \
     " #{gff.ext('.pep.fa')}"
end

def read_fa(file)
  seqs = {}
  IO.foreach(file).each_slice(2) do |id, seq|
    id = id.split[0][1..-1]
    seqs[id] = seq.chomp
  end
  seqs
end

send "extract_#{type}", source_fa, source_gff
send "extract_#{type}", target_fa, target_gff
inp_seqs = read_fa source_gff.ext(".#{type}.fa")
out_seqs = read_fa target_gff.ext(".#{type}.fa")

inp_seqs.each do |id, inp_seq|
  out_seq = out_seqs[id]
  if out_seq && inp_seq == out_seq
    # mapped - do nothing
  else
    # un-mapped
    puts id
  end
end
