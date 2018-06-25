#!/usr/bin/env ruby
# Copyright 2017 Anurag Priyam, Queen Mary University of London
#
# Removes features of the given type from the given gff file.
# Subfeature are not removed. Instead the subfeature's
# Parent atribute is unset.
#
# This script must be run by the user, if desired,
# before running flo - see README.

feat_type_to_remove, gff_file = ARGV

unless feat_type_to_remove and gff_file
  puts <<MSG
Usage: gff_remove_feats.rb feature_type_to_remove gff_file > output_file

e.g. gff_remove_feats.rb gene all_annotations.gff > transcripts_only.gff
MSG
  exit!
end

# Collect id of the features we want to remove.
feat_ids_to_remove = []
IO.foreach(gff_file) do |line|
  next if line[0] == '#'

  line = line.split
  if line[2] == feat_type_to_remove
    feat_id = line[8].match(/ID=(.+?);/)[1]
    feat_ids_to_remove << feat_id
  end
end

# Go over the file again. Skip over the lines containing
# ID=feat_ids_to_remove;. For the remaining lines,
# gsub Parent=feat_ids_to_remove;.with ''.
feat_ids_to_remove = feat_ids_to_remove.join('|')
pt_regex = /Parent=(#{feat_ids_to_remove});/
id_regex = /ID=(#{feat_ids_to_remove});/
IO.foreach(gff_file) do |line|
  next if line.match id_regex
  puts line.gsub pt_regex, ''
end
