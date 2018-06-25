#!/usr/bin/env ruby
# Copyright 2017 Anurag Priyam, Queen Mary University of London
#
# Retain the longest transcript for each gene. Gene annotations
# and Parent attribute of the longest mRNA retained are removed
# from the output.

require 'set'
begin
  require 'bio/db/gff'
rescue LoadError
  # Could not load BioRuby. Ask to install BioRuby and exit.
  puts <<MSG
Please install the bio gem first:

  sudo gem install bio
MSG
  exit!
end

# GFF file to process.
gff_file = ARGV.pop

# Print help and exit if gff_file not specified.
unless gff_file
  puts <<MSG
Usage: gff_longest_transcript.rb input_gff_file > output_gff_file
MSG
  exit!
end

# Read the gff file into memory and parse it.
gff3 = Bio::GFF::GFF3.new File.read(gff_file)

# Select all transcripts from the GFF file, grouped by their parent gene's id.
all_transcripts = Hash.new { |k,v| k[v] = []}
gff3.records.each do |record|
  # Transcript are usually annotated as 'mRNA' or 'transcript' in their 3rd
  # column.
  if record.respond_to?(:feature_type) &&
      record.feature_type =~ /mRNA|transcript/
    all_transcripts[record.attributes.assoc('Parent').last] << record
  end
end

# Build a set of ids of longest transcript for each gene. Set is used for
# quick lookups afterwards.
longest_transcript_ids = Set.new
all_transcripts.each do |gene, transcripts|
  # Select longest transcript.
  longest_transcript = transcripts.sort_by do |transcript|
    transcript.end - transcript.start
  end.last
  # Add its ID to the set.
  longest_transcript_ids << longest_transcript.attributes.assoc('ID').last
end

# Go through GFF again printing records whose ID or Parent is included in
# longest_transcript_ids Set.
gff3.records.each do |record|
  record_id = record.attributes.assoc('ID')
  record_pr = record.attributes.assoc('Parent')
  if record_id && longest_transcript_ids.include?(record_id.last) ||
      record_pr && longest_transcript_ids.include?(record_pr.last)
    puts record
  end
end
