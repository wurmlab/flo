# Remove annotations on the mitochondrial genome from the GFF.
grep -v 'NC_014672\.1' ref_Si_gnG_scaffolds.gff3 > ref_Si_gnG_scaffolds.no_mt.gff3
