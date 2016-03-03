#!/usr/bin/env ruby

# Copyright 2016 Anurag Priyam - MIT License
#
# Filters out mRNA that have N in their exonic sequence, or whose coding
# sequence is not a multiple of three, or have partial=true in its list of
# attributes, or have a Note indicating that the underlying sequence has been
# modified.

require_relative 'helpers'

gff3, exonic_fa, coding_fa = ARGV
mrna_subfeatures = parse_gff3 gff3

# filter_records returns [pass, fail]. Since rules below are constructued to
# pass for mRNA we don't want, we collect them as [rejected, selected].
rejected, _ = filter_records(mrna_subfeatures, exonic_fa, coding_fa,
  "'N' in exonic sequence"              => -> { exonic_seq =~ /N+/ },
  'CDS not a multiple of 3'             => -> { coding_seq.length % 3 != 0 },
  'partial=true in list of attributes'  => -> { [mrna].concat(subfeatures).any? { |f| attribute(f, 'partial') == 'true' } },
  "Note attributes mentions 'modified'" => -> { [mrna].concat(subfeatures).any? { |f| attribute(f, 'Note') =~ /modified/ } },
)

rejected.each do |rule, mrnas|
  $stderr.puts "#{rule}: #{mrnas.length}"
end

all_filtered = mrna_subfeatures.to_a.flatten - rejected.values.flatten.uniq
write_gff3(all_filtered.flatten)
