#!/usr/bin/env ruby

require 'bio/db/gff'

# Obtain all annotations from the lifted GFF3 file.
records = Bio::GFF::GFF3.new(File.read('lifted.gff3')).records
records.reject! { |rec| rec.class != Bio::GFF::GFF3::Record }

# Group together sibling CDS and their mRNA parent.
records = records.group_by do |rec|
  key = 'ID'     if rec.feature_type == 'mRNA'
  key = 'Parent' if rec.feature_type == 'CDS'
  rec.attributes.assoc(key).last
end

a = 0
b = 0
records = records.map do |id, annots|
  a += 1 unless annots.map(&:seqname).uniq.length == 1    # mRNA on different scaffolds/contigs
  b += 1 unless annots.map(&:feature_type).include? 'CDS' # mRNA with no CDS
end

puts a
puts b
