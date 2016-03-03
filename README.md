# flo - same species annotations lift over pipeline
##### Anurag Priyam¹, Rodrigo Pracana¹, Yannick Wurm¹
##### ¹ Queen Mary University of London
##### Contact: a.priyam@qmul.ac.uk, y.wurm@qmul.ac.uk

## Introduction

"Lift over" is the process of transferring annotations from one genome assembly to another. Generally lift over is done to transfer annotations for a species from an old genome assembly to a new genome assembly. Lift over involves identifying blocks of highly similar sequences (synteny) between the two assemblies and using the synteny information to transform annotation coordinates to the new assembly. Using NCBI Remap \[[1](#ref1)\] and Ensemble API \[[2](#ref2)\] it is possible to lift over annotations between a few publicly available assemblies. `liftOver` \[[3](#ref3)\] or CrossMap \[[4](#ref4)\] can be used if synteny information is available as a chain file (which can be downloaded from the UCSC genome website for a few publicly available assemblies), but fail if annotations are in the popular GFF3 format: neither guarantees that the lifted annotations will make a valid gene model.

We developed flo, an end-to-end same-species lift over pipeline for annotations in GFF3 format. Lifted annotations are processed to ensure they are biologically meaningful and form a valid GFF3 file. Chain file created by flo can be used with other tools to lift over annotations in a variety of file format.

## Usage

To use `flo` you must have Ruby 2.0 or higher.

#### Download flo

    wget -c https://github.com/yeban/flo/archive/master.tar.gz -O flo.tar.gz
    tar xvf flo.tar.gz
    mv flo-master flo
    cd flo

#### Install flo's dependencies

1. The following utilities from [UCSC Genome site][1]: `faSplit`, `faToTwoBit`,
`twoBitInfo`, `blat`, `axtChain`, `chainSort`, `chainMergeSort`, `chainSplit`,
`chainNet`, `netChainSubset`, and `liftOver` (>= v315 / 29th April, 2015).
2. [GNU Parallel][2] (>= 20150422).
3. [GenomeTools][3] (>= 1.5.6).

On Linux, running `scripts/install.sh` will pull all of the above into `ext/`
directory.

#### Tell flo about your data

First, copy the example configuration file, `opts_example.yaml`, to
`opts.yaml`.

    cp opts_example.yaml opts.yaml

Then, edit `opts.yaml` to indicate: the location of utilities downloaded from
UCSC Genome site, GNU Parallel, and GenomeTools;location of source and target
assembly in FASTA format; location of GFF3 file(s) containing the annotations
on the source assembly and number of CPU cores to use.

If no GFF3 files are specified, `flo` stops after creating a chain file which
can be used with `liftOver` or `CrossMap` to lift over annotations in a
variety of file formats.

#### Run flo

Just invoke `rake` command (a utility bundled with Ruby).

    rake

#### Output

This will generate a directory corresponding to each GFF file:

```
.
├── input_gff_name-liftover-target_assembly_name
│   ├── input_gff_name-liftover-target_assembly_name.gff3   -> Final, processed, GFF3.
│   ├── ...                                                 -> Intermediate files.
│   └── summary.txt                                         -> Lift over summary.
|
├── ...               -> Same as above for additional GFF files.
|
└── run
    ├── liftover.chn  -> Chain file used for lift over.
    ├── ...           -> Intermediate files.
```

## Approach

`flo` is an implementation of the lift over approach as documented online \[[5](#ref5)\] by Jim W. Kent, with a few improvements. A high-level overview of the approach is provided below.

The target assembly is aligned to the source assembly using BLAT \[[6](#ref6)\].  Gapless blocks of aligned regions between pairs of target and source sequences are then grouped together to form "chains" using `axtChain` \[[7](#ref7)\]. Chains are like alignments, excpet they allow larger gaps between blocks and the gaps may occur simultaneously on both sequences. Subsequently, a hierarchy of chains, ordered by score, is built in which each chain lower down the hierarchy covers gap(s) in the chain above it. Parts of low-scoring chains that redundantly cover high-scoring chains are discarded. This is accomplished with `chainNet` \[[8](#ref8)\]. Finally, a subset of all chains that are represented in the net is derived using `netChainSubset` (distributed with UCSC genome utilities) and used for lift over.

Using `liftOver` \[[3](#ref3)\] and the chain file derived above, annotation coordinates (sequence id, start and stop coordinates; first three columns of a GFF3 file) are converted to the target assembly. Lifted annotations that do not form a valid gene model are eliminated and the resulting annotations processed with GenomeTools \[[8](#ref8)\] to correct the phase of CDS annotations and ensure the output is a well-formed, valid GFF3 file. Finally, a summary indicating the number of transcripts in the input and output set, and those that yield identical coding sequence is generated.

`flo` is implemented using the `rake` \[[9](#ref9)\] utility in the Ruby programming language \[[10](#ref10)\]. Lift over is computationally intensive. Thus the target assembly is split across several files, into chunks of 5000 nucleotides or smaller (required by `blat -fastMap`), and each file processed in parallel with GNU Parallel \[[11](#ref11)\]. Further, BLAT is run with `-fastMap -tileSize=12 -minIdentity=98` options to provide an optimal trade-off between sensitivity and performance, but the user may override BLAT options.

## Discussion

To evaluate `flo`, we considered that lifting annotations to the same assembly should yield annotations identical to the input set.  Secondly, to measure the success of `flo` in lifting annotations from one assembly to another, we considered that the input and output mRNA should yield identical coding sequence. We used NCBI's fire ant annotation release 100 \[[12](#ref12)\] on the genome assembly gnG and the fire ant genome assembly gnH (Riba-Grognuz et al., unpublished) as our test dataset. Indeed, lifting annotations from gnG to gnG yields the exact same annotations as the input set, while, lifting annotations from gnG to gnH lifted 94.52% of the annotations in the input set, 91.64% of which yielded identical coding sequences (86.61% of the input set).

| Assembly   | Input  | Lifted | High confidence |
| --------   | -----  | ------ | --------------- |
| gnG to gnG | 14,453 | 14,453 | 14,453          |
| gnG to gnH | 14,453 | 13,661 | 12,519          |

Due to differences in the sequence assemblies, parts of target assembly may not align or align weakly to the source assembly: annotations lying in such regions cannot be lifted.  Similarly, presence of 'NNN' and repeated regions in the assemblies can interrupt alignment and cause an annotation to be only partially covered by a chain: such annotations are not lifted. Further, duplicated regions in the assembly can cause parts of an mRNA to map to different scaffolds: such annotations are eliminated.

Lift over is computationally intensive: the ~350 Mbp fire ant genome took about 30 mins on a 40-core machine.

Identifying syntenic blocks between assemblies is independent of the annotations, and thus the chain file created by `flo` can be use to lift over annotations in format such as GTF, BED, GenePred, BAM, and VCF using `liftOver` \[[3](#ref3)\] and CrossMap \[[4](#ref4)\].

## References

<a id="ref1" href="http://www.ncbi.nlm.nih.gov/genome/tools/remap">1. NCBI. Genome remapping service. Accessed 7th May 2015. http://www.ncbi.nlm.nih.gov/genome/tools/remap.</a>

<a id="ref2" href="http://asia.ensembl.org/Homo_sapiens/Tools/AssemblyConverter">2. Ensembl. Assembly converter. Accessed 7th May 2015. http://asia.ensembl.org/Homo_sapiens/Tools/AssemblyConverter.</a>

<a id="ref3" href="http://bib.oxfordjournals.org/content/14/2/144.full">3. Kuhn RM, Haussler D, and Kent WJ. The UCSC genome browser and associated tools. Briefings in Bioinformatics. 2013; 14(32): 144–161.</a>

<a id="ref4" href="http://bioinformatics.oxfordjournals.org/content/30/7/1006.short">4. Zhao H, Sun Z, Wang J, Huang H, Kocher JP, and Wang L. CrossMap: a versatile tool for coordinate conversion between genome assemblies. Bioinformatics. 2014; 30(7): 1006-1007.</a>

<a id="ref5" href="http://hgwdev.cse.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt">5. Kent WJ. liftOver.txt. Accessed 7th May 2015. http://hgwdev.cse.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt.</a>

<a id="ref6" href="http://genome.cshlp.org/content/12/4/656">6. Kent WJ. BLAT - The BLAST-like alignment tool. Genome Research. 2002;12:656-664.</a>

<a id="ref7" href="http://www.pnas.org/content/100/20/11484.full">7. Kent WJ, Baertsch R, Hinrichs A, Miller W, Haussler D. Evolution's cauldron: Duplication, deletion, and rearrangement in the mouse and human genomes. Proceedings of the National Academy of Sciences of the United States of America. 2003; 100(20): 11484-11489.</a>

<a id="ref8" href="http://www.computer.org/csdl/trans/tb/2013/03/ttb2013030645-abs.html">8. Gremme G, Steinbiss S, and Kurtz S. GenomeTools: a comprehensive software library for efficient processing of structured genome annotations. IEEE/ACM Transactions on Computational Biology and Bioinformatics. 2013; 10(3): 645–656.</a>

<a id="ref9" href="https://github.com/ruby/rake">9. Weirich J et al. RAKE – Ruby Make. Accessed 7th May 2015. https://github.com/ruby/rake.</a>

<a id="ref9" href="https://www.ruby-lang.org/">10. Flanagan D, Matsumoto Y. The Ruby programming language. 1st ed. Beijing: O’Reilly; 2008.</a>

<a id="ref11" href="http://www.gnu.org/software/parallel/">11. Tange O. GNU Parallel - The Command-Line Power Tool. The USENIX Magazine. 2011; 36: 42-47</a>

<a id="ref12" href="http://www.ncbi.nlm.nih.gov/genome/annotation_euk/Solenopsis_invicta/100/">12. NCBI. Fire ant annotation release 100. Accessed 20th April 2015. http://www.ncbi.nlm.nih.gov/genome/annotation_euk/Solenopsis_invicta/100/.</a>

[1]: http://hgdownload.cse.ucsc.edu/admin/exe/
[2]: http://www.gnu.org/software/parallel/
[3]: http://genometools.org/
