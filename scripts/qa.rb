#!/usr/bin/env ruby

# blastn -task blastn -query cds.mapped.fa -db cds.orig.fa -outfmt '6 qseqid sseqid evalue qcovs pident' -num_threads 20 -evalue '1e-200' > out.tsv

lines = File.readlines 'qa.tsv'
lines.map!(&:split)

all = []
good = {}
lines.each do |line|
  qid = line[0]
  evl = line[2].to_f
  cov = line[3].to_f
  pid = line[4].to_f
  qlen = line[5].to_i
  slen = line[6].to_i
  alen = line[7].to_i

  all << qid
  good[qid] ||= 0
  good[qid] = good[qid] + 1 if pid == 100 && cov == 100 && qlen == slen
end

puts (all.uniq - Hash[good.select{|k, v| v > 0}].keys.uniq)
