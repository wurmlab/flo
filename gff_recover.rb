#!/usr/bin/env ruby
# Copyright 2017 Anurag Priyam, Queen Mary University of London
#
# Tries to reconstruct transcripts from a badly made gff files - such as
# liftOver's output.
#
# Transcripts and their exons and cds are first grouped together based
# on ID attribute of the transcript (if present) and Parent attribute
# of exons and CDS features. If the annotations in a group are not all
# on the same reference sequence, the group is eliminated. If a group
# of annotations lacks a transcript annotations, it is calculated.
#
# The gff file to be recovered is read into memory.

begin
  require 'bio/db/gff'
rescue LoadError
  puts <<MSG
Please install the bio gem first:

  sudo gem install bio
MSG
  exit!
end

# We can process 2-level features: transcripts and their subfeatures.
transcript_types = /gene|mRNA|transcript/
subfeature_types = /exon|CDS/

# Read the lifted gff file into memory and parse it.
gff3 = Bio::GFF::GFF3.new File.read(ARGV.pop)

# Obtain transcripts and their children. genes and other features are not
# processed.
record_groups = Hash.new { |h, k| h[k] = [] }
gff3.records.each do |record|
  # GFF file includes features, comments and directives. We are only
  # interested in "features".
  next unless record.respond_to?(:feature_type)

  # Consider ID attribute of transcripts, and Parent attribute of
  # exon and CDS.
  key = case record.feature_type
        when transcript_types
          'ID'
        when subfeature_types
          'Parent'
        end

  # If the annotation is neither a transcript nor a subfeature
  # type, key will not be set and this may cause assoc to
  # raise an error.
  val = record.attributes.assoc(key) if key

  # If the annotation is neither a transcript nor a subfeature
  # type, val will not be set. We print this feature to stderr.
  unless val
    $stderr.puts record
    next
  end

  # val is a 2-tuple: ['ID', value] or ['Parent', value]
  record_groups[val.last] << record
end

# How are transcripts annotated?
transcript_type = record_groups.each { |id, records|
  transcript = records.find { |record|
    record.feature_type =~ transcript_types
  }
  break transcript.feature_type if transcript
}

# Reject invalid groups. Print others.
record_groups.each do |id, records|
  # liftOver's output may contain sibling exons and CDS
  # mapped to different reference sequences. This is
  # simply because liftOver lifts features (even
  # subfeatures) one by one. Skip past those.
  next if records.map(&:seqname).uniq.length > 1

  # For single exon transcripts, it is possible to end
  # up with transcripts without CDS annotations. Skip
  # past those.
  next unless records.map(&:feature_type).include? 'CDS'

  # If transcript annotation wasn't lifted it's because
  # one or more exon and cds couldn't be lifted are were
  # lifted partially. Add a transcript annotation around
  # the exon and cds that could be lifted.
  transcript = records.find { |record|
    record.feature_type =~ transcript_types
  }
  unless transcript
    puts Bio::GFF::GFF3::Record.new(*[
      records.first.seqname,
      records.first.source,
      transcript_type,
      records.map(&:start).min,
      records.map(&:end).max,
      nil,
      records.first.strand,
      nil,
      [["ID", id]]
    ])
  end

  puts records
end
