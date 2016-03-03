# Copyright 2016 Anurag Priyam - MIT License
#
# Provides powerful helper functions for working with GFF3 and FASTA files
# using BioRuby's Bio::GFF::GFF3 and GenomeTools.

require 'rake'
require 'ostruct'
require 'tempfile'
require 'bio/db/gff'

def sequences(file)
  (@sequences ||= {})[file] ||= Hash[
    File.readlines(file).each_slice(2).
      map { |idl, seq| [ idl.split[0][1..-1], seq.chomp ] }
  ]
end

# Fetch sequence with the given id from the given FASTA file from the cache, or
# extract the sequence from the file, and cache and return it.
#
# NOTE: Assumes sequences aren't split across several lines.
def get_sequence(id, file)
  sequences(file).fetch(id)
rescue => e
  $stderr.puts "Error obtaining sequence: #{id}, #{file}"
  $stderr.puts "\t#{e}"
  exit
end

# Extract value of the given attribute for the given Bio::GFF:GFF3 feature
# record.
def attribute(feat, attr)
  val = feat.attributes.assoc(attr)
  val.last if val
end

# GT is buggy! Following function processes the given GFF3, and then processes
# the processed GFF3 and so on until the ouput equals input. At this point we
# can be sure GT has done everytying it can do. Processing is required for
# extractseq to work.
def process_with_gt(gff3, gt_gff3)
  tmp = Tempfile.new('gt', '.'); tmp.close
  sh "gt gff3 -tidy -sort -retainids -addids #{gff3} > #{tmp.path}"
  sh "diff #{tmp.path} #{gff3} &> /dev/null" do |ok, _|
    if !ok
      process_with_gt(tmp.path, gt_gff3)
    else
      tmp
    end
  end
end

# Returns mRNAs and their subfeatures from the given GFF3 file as:
#
#   {
#     mrna => [exons, cds],
#     ...
#   }
#
# NOTE: Assumes mRNA annotations will appear before its children (exon and
# CDS). One way to ensure that is to sort the GFF file with GenomeTools.
def parse_gff3(name)
  mrna_subfeatures = {}; mrnas = {}
  Bio::GFF::GFF3.new(File.read(name)).records.each do |rec|
    next unless rec.respond_to? :feature_type
    parent_id = attribute(rec, 'Parent')
    record_id = attribute(rec, 'ID')
    case rec.feature_type
    when 'mRNA'
      mrnas[record_id] = rec
      mrna_subfeatures[rec] = []
    when /exon|CDS/
      mrna = mrnas[parent_id]
      mrna_subfeatures[mrna] << rec if mrna
    end
  end
  mrna_subfeatures
end

# Runs a set of procs against each feature group. Returns two Hashes, the first
# containing feature groups for which the proc returned true and the second
# containing feature groups for which the proc returned false.
def filter_records(mrna_subfeatures, exonic_fa, coding_fa, filters)
  yay = Hash.new { |h, k| h[k] = [] }
  nay = Hash.new { |h, k| h[k] = [] }
  mrna_subfeatures.each do |mrna, subfeatures|
    mrna_id = attribute(mrna, 'ID')
    data = OpenStruct.new(mrna: mrna, subfeatures: subfeatures,
                          exonic_seq: get_sequence(mrna_id, exonic_fa),
                          coding_seq: get_sequence(mrna_id, coding_fa))
    filters.each do |name, filter|
      if data.instance_exec(&filter)
        yay[name] << [mrna, subfeatures]
      else
        nay[name] << [mrna, subfeatures]
      end
    end
  end
  [yay, nay]
end

# Write given GFF3 records to stdout. Writes the records to a temporary file,
# which is then processed with GT the output of which is written to stdout.
def write_gff3(recs)
  temp = Tempfile.open('temp', '.')
  temp.write("##gff-version 3\n\n")
  temp.write(recs.map(&:to_s).join)
  temp.close
  sh "gt gff3 -tidy -sort -retainids -addids #{temp.path} 2> /dev/null"
end
