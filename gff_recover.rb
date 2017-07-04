#!/usr/bin/env ruby
# Copyright 2017 Anurag Priyam, Queen Mary University of London

require 'bio/db/gff'

# Addresses the following issues in the given GFF file.
#
# a. mRNA with subfeatures on different scaffolds. Such annotations are
#    removed.
# b. CDS with no mRNA. An mRNA is added around these CDS.
# c. mRNA with no CDS. Such mRNAs are removed.
def process_gff(gff)
  # Read the lifted gff file into memory and parse it.
  gff3 = Bio::GFF::GFF3.new File.read(gff)

  # Obtain transcripts and their children. genes and other features are not
  # processed.
  transcripts = Hash.new { |h, k| h[k] = [] }
  gff3.records.each do |record|
    # GFF file includes features, comments and directives. We are only
    # interested in "features".
    next unless record.respond_to?(:feature_type)

    # If the feature is a transcript, we consider its ID attribute.
    if record.feature_type =~ /gene|mRNA|transcript/
      transcripts[record.attributes.assoc('ID').last] << record
    end

    # If the feature is exon or CDS, we consider its Parent attribute.
    if record.feature_type =~ /exon|CDS/
      transcripts[record.attributes.assoc('Parent').last] << record
    end
  end

  transcripts = transcripts.map do |id, annots|
    next unless annots.map(&:seqname).uniq.length == 1    # mRNA on different scaffolds/contigs
    next unless annots.map(&:feature_type).include? 'CDS' # mRNA with no CDS

    # If there's a group of CDS without parent mRNA.
    if annots.map(&:feature_type).grep(/gene|mRNA|transcript/).empty?
      mrna = [
        annots.first.seqname,
        annots.first.source,
        'mRNA',
        annots.map(&:start).min,
        annots.map(&:end).max,
        nil,
        annots.first.strand,
        nil,
        [["ID", id]]
      ]
      mrna = Bio::GFF::GFF3::Record.new(*mrna)
      annots.unshift(mrna)
    end

    annots
  end.flatten.compact
end

puts process_gff(ARGV.pop).join
