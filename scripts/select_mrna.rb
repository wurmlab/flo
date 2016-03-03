#!/usr/bin/env ruby

# Copyright 2015 Anurag Priyam - MIT License
#
# Prints all mRNA and there children from the given GFF3 file.  If 'longest' is
# specified as second argument, longest mRNA per gene is printed. If a FASTA of
# coding sequences is specified as the second argument, mRNA with longest
# coding sequence per gene is printed.
#
# As the parent gene of mRNA is not going to be in the output, Parent attribute
# of mRNA are changed to _Parent so that the output is valid GFF3, and still be
# able to track parent gene.

require_relative 'helpers'

gff3, coding_fa = ARGV
mrna_subfeatures = parse_gff3 gff3

# Group mRNA by there parent id, simultaneously deleting parent attribute
# from mRNA records (required for output GFF3 to be valid).
grouped_by_genes = mrna_subfeatures.group_by do |mrna, subfeatures|
  if mrna_parent_attr = mrna.attributes.assoc('Parent')
    mrna.attributes << ['_Parent', mrna_parent_attr[1]]
    mrna.attributes.delete(mrna_parent_attr)
  else
    attribute(mrna, '_Parent')
  end
end

if !coding_fa
  # All mRNA.
  write_gff3(mrna_subfeatures.to_a.flatten)

elsif coding_fa == 'longest'
  # Longest mRNA.
  longest_mrna = grouped_by_genes.map do |_, mrna_subfeatures|
    mrna_subfeatures.sort_by do |mrna, _|
      mrna.end - mrna.start
    end.last
  end
  write_gff3(longest_mrna.flatten)

elsif File.exist?(coding_fa)
  # mRNA with longest CDS.
  longest_coding_mrna = grouped_by_genes.map do |_, mrna_subfeatures|
    mrna_subfeatures.sort_by do |mrna, _|
      get_sequence(attribute(mrna, 'ID'), coding_fa).length
    end.last
  end
  write_gff3(longest_coding_mrna.flatten)
end
