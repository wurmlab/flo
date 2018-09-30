#!/usr/bin/env ruby
# Copyright 2017 Anurag Priyam, Queen Mary University of London
#
# Removes features of the given type from the given gff file.
# Subfeature are not removed. Instead the subfeature's
# Parent atribute is unset.
#
# This script must be run by the user, if desired,
# before running flo - see README.

require 'set'
require_relative 'gff_helpers'

# Take as input the feature type to remove and the GFF file to remove the
# features from.
feat_type_to_remove, gff_file = ARGV

# Print usage and exit if required input is not provided.
unless feat_type_to_remove and gff_file
  puts <<MSG
Usage: gff_remove_feats.rb feature_type_to_remove gff_file > output_file

e.g. gff_remove_feats.rb gene all_annotations.gff > transcripts_only.gff
MSG
  exit!
end

# Initialise a Set to hold id of the features we want to remove. We use Set
# for fast, subsequent lookup.
feat_ids_to_remove = Set[]

# Go over each line of the input GFF file and collect id of the features we
# want to remove (based on feature type).
IO.foreach(gff_file) do |line|
  # Skip comments.
  next if line[0] == '#'

  # Process feature lines.
  cols = line.split
  # If this feature is of the type we want to remove, add it to the Set.
  if cols[2] == feat_type_to_remove
    feat_ids_to_remove << get_attr(cols[8], 'ID')
  end
end

# Go over the input GFF file again. Output comments as it is. Eliminate
# features whose id we made a note of previously, and delete Parent
# attribute of the features whose parent feature was removed.
IO.foreach(gff_file) do |line|
  # Output comment lines as it is.
  if line[0] == '#'
    puts line
    next
  end

  # Process feature lines.
  cols = line.split

  # Get feature id and check if it was marked for removal. If so, skip to the
  # next line.
  feat_id = get_attr(cols[8], 'ID')
  next if feat_ids_to_remove.include? feat_id

  # Get id of feature's parent and check if it was marked for removal. If so,
  # delete the Parent attribute.
  feat_parent_id = get_attr(cols[8], 'Parent')
  if feat_ids_to_remove.include? feat_parent_id
    cols[8] = del_attr(cols[8], 'Parent')
  end

  # Print the line back.
  puts cols.join("\t")
end
