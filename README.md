# flo - same species annotations lift over pipeline

Lift over is the process of transferring annotations from one genome assembly
to another. Usually lift over is done because there is a new, improved genome
assembly for the species but good quality annotations (maybe manually curated
or experimentally verified) are available on the old assembly.

The idea is simple: align the new assembly with the old one (e.g., with BLAT),
process the alignment data to define how a coordinate or coordinate range on
the old assembly should be transformed to the new assembly (e.g., as a chain
file), transform the coordinates (e.g., with liftOver).

Tools to lift over annotations are available from [UCSC-Kent toolkit][1]. flo
combines them together into a pipeline that is easy to run, tweak, and extend.
Additionally, flo validates and processes the lifted GFF file to ensure the
annotations are biologically meaningful, something liftOver doesn't do.

[Usage](#usage) | [Approach](#approach) | [Discussion](#discussion)

## Usage

To use `flo` you must have Ruby 2.0 or higher.

#### Download flo

    wget -c https://github.com/yeban/flo/archive/master.tar.gz -O flo.tar.gz
    tar xvf flo.tar.gz
    mv flo-master flo
    cd flo

#### Install flo's dependencies

1. [UCSC-Kent toolkit][1].
2. GNU Parallel - to easily parallelize different steps of flo.
3. genometools - to validate, post-process, and extract sequences from the
   lifted GFF file.
4. BioRuby gem - for parse GFF files for post-processing after lift over.

Running `scripts/install.sh` will install 1-3 into `ext/` directory. For
BioRuby gem, please run `gem install bio` or `sudo gem install bio`.

#### Tell flo about your data

First, copy the example configuration file (`opts_example.yaml`) to
`opts.yaml`.

    cp opts_example.yaml opts.yaml

Then, edit `opts.yaml` to indicate: the location of UCSC-Kent toolkit, GNU
Parallel, and genometools; location of source and target assembly in FASTA
format; location of GFF3 file(s) containing the annotations on the source
assembly and number of CPU cores to use.

#### Run flo

    rake

This will generate a directory corresponding to each GFF file:
"input_gff_name-liftover-target_assembly_name".

The above dir contains a similarly named GFF file which is the final output.
Additionally, if flo finds a input_gff_name.cdna/cds/pep.fa alongside the
input GFF file, flo will automatically generate corresponding FASTA from
the lifted annotations and print a comparison summary (how many sequences
in input, how many in output, how many identical).

## Approach

### Full genome alignment

This is done with BLAT. BLAT is fast, but doing a full genome can still take a
significant amount of time. So the target assembly is split into chunks of 5k
nucleotides. The chunks are then clubbed into "n" groups, depending on the
number of CPU cores on the system, and BLATed parallely against the source
assembly.

Because it was chunks of the new assembly that was BLATed against the old
assembly the alignment coordinates will have to be changed back to how it
would look if the entire new assembly were BLATed against the old assembly.
This process is called "lift up".

Lift up is done with the `liftUp` tool from UCSC-Kent toolkit and with the help
of a lift up (`.lft`) file. A lift up file is generated with the `-lift` option
of the `faSplit` tool used for chunking and grouping the target assembly.

### Chaining, netting, filtering

`.psl` files containing alignment data from BLAT are converted into a [chain
file][2]. This is done with the `axtChain` tool. `axtChain` finds pairwise
alignments in the psl files (one between the same target and query sequence)
and joins them together, if overlapping and combining them will result in a
higher scoring longer alignment.

The resulting chain files are then sorted, merged, and then split with new
assembly as the reference. This bit doesn't make sense to me either, but is
important.

The resulting chain file is then "netted" using `chainNet` tool. Netting
organizes the chains in a hierarchal collection with the highest-scoring
non-overlapping chains on top, and their gaps filled in where possible by
lower-scoring chains, which in turn may have gaps filled in by lower-level
chains and so on.

The chain file in the previous step is then filtered against the net file
obtained above to get a lift over file that can be used for coordinate
transformation.

### Lift over

With the above lift over file annotations on the old assembly can be
transferred to the new assembly using `liftOver` tool.

If lifting annotations in GFF format, one can run into issues like:

1. mRNA mapped to different scaffolds.
2. mRNA with no CDS.
3. CDS with no mRNA.
4. Duplicated CDS ids (not sure why this happens).

So the resulting GFF should be processed further.

## Discussion

`flo` in itself doesn't do anything new or novel. It merely combines existing
tools into a usable fashion.

We used `flo` to migrate to do same species lift over with ~350Mb genome size.
94% of the annotations were mapped to the final assembly. Of the mapped
annotations, 90% were good.

## What if annotations are not in GFF format

You could convert to GFF or tweak the "lift over" part of flow. If your
annotations are in bed or gp format, you should consider tweaking `flo`
instead of converting to GFF.

If you want to use BAM format check out [`CrossMap`][4].

## I have a big genome. It's taking a very long time to run BLAT.

Try creating an ooc file.

```
sh "blat ../../data/#{SOURCE}.fa /dev/null /dev/null" \
     " -tileSize=11 -makeOoc=11.ooc -repMatch=100"

# ...

"blat -noHead -fastMap -ooc=11.ooc -minScore=100 -minIdentity=98" \
```

You can do this in the `joblist` task.

## Can I use flo for different species lift over / creating chain files

Not as it is. But you perhaps reuse the framework and much of the alignment and
chaining-netting step. blastz or another may be more suitable instead of blat.

[1]: http://hgdownload.cse.ucsc.edu/admin/exe/
[2]: http://genome.ucsc.edu/goldenpath/help/chain.html
[3]: http://genometools.org/
[4]: http://crossmap.sourceforge.net/
