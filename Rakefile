# Copyright 2017 Anurag Priyam - MIT License
#
# Annotation lift over pipeline.
#
# Based on the lift over procedure deseribed at:
# http://genomewiki.ucsc.edu/index.php/LiftOver_Howto &
# http://hgwdev.cse.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt
#
# Additional references:
# http://genomewiki.ucsc.edu/index.php/Chains_Nets
# https://genome.ucsc.edu/goldenPath/help/net.html
# http://genome.ucsc.edu/goldenpath/help/chain.html
# http://asia.ensembl.org/info/website/upload/psl.html
#
# The pipeline depends on GNU parallel, genometools (> 1.5.5) and the following
# tools from UCSC-Kent tookit: faSplit, faToTwoBit, twoBitInfo, blat, axtChain,
# chainSort, chainMergeSort, chainSplit, chainNet, netChainSubset, and liftOver.

require 'yaml'
require 'tempfile'

# Reads config file. Runs task to create chain file, then runs liftOver.
task 'default' do
  # Check for presence of config file. Exit if not found.
  unless File.exist? 'flo_opts.yaml'
    puts "Config file not found. See README for how to use flo."
    exit!
  end

  # Read config file.
  CONFIG = YAML.load_file 'flo_opts.yaml'

  # Add dirs specified in config to PATH.
  Array(CONFIG[:add_to_path]).each { |path| add_to_PATH path }

  # Create chain file.
  task('run/liftover.chn').invoke

  # Lift over the given GFF3 files.
  Array(CONFIG[:lift]).each do |inp|
    outdir = "run/#{inp.ext}"
    mkdir outdir

    # Lift over the annotations from source assembly to target assembly.
    sh "liftOver -gff #{inp} run/liftover.chn #{outdir}/lifted.gff3" \
       " #{outdir}/unlifted.gff3"

    # Clean lifted annotations.
    sh "#{__dir__}/gff_recover.rb #{outdir}/lifted.gff3" \
       " | gt gff3 -tidy -sort -addids -retainids -"     \
       " > #{outdir}/lifted_cleaned.gff"

    # Symlink input gff to outdir and summarise / validate.
    sh "ln -s #{File.expand_path inp} #{outdir}/input.gff"
    summarize outdir
  end
end

# Task to create chain file.
file 'run/liftover.chn' do
  mkdir 'run'

  processes = CONFIG[:processes]
  blat_opts = CONFIG[:blat_opts]

  cp CONFIG[:source_fa], 'run/source.fa'
  cp CONFIG[:target_fa], 'run/target.fa'

  to_2bit 'run/source.fa'
  to_2bit 'run/target.fa'

  to_sizes 'run/source.2bit'
  to_sizes 'run/target.2bit'

  # Partition target assembly.
  sh "faSplit sequence run/target.fa #{processes} run/chunk_"

  parallel Dir['run/chunk_*.fa'],
    'faSplit -oneFile size %{this} 5000 %{this}.5k -lift=%{this}.lft &&'       \
    'mv %{this}.5k.fa %{this}'

  # BLAT each chunk of the target assembly to the source assembly.
  parallel Dir['run/chunk_*.fa'],
    "blat -noHead #{blat_opts} run/source.fa %{this} %{this}.psl"

  parallel Dir['run/chunk_*.fa'],
    "liftUp -type=.psl -pslQ -nohead"                                          \
    " %{this}.psl.lifted %{this}.lft warn %{this}.psl"

  # Derive a chain file each from BLAT's .psl output files.
  parallel Dir['run/chunk_*.psl.lifted'],
    'axtChain -psl -linearGap=medium'                                          \
    ' %{this} run/source.2bit run/target.2bit %{this}.chn'

  # Sort the chain files.
  parallel Dir["run/chunk_*.chn"],
    'chainSort %{this} %{this}.sorted'

  # Combine sorted chain files into a single sorted chain file.
  sh 'chainMergeSort run/*.chn.sorted | chainSplit run stdin -lump=1'
  mv 'run/000.chain', 'run/combined.chn.sorted'

  # Derive net file from combined, sorted chain file.
  sh 'chainNet'                                                                \
     ' run/combined.chn.sorted run/source.sizes run/target.sizes'              \
     ' run/combined.chn.sorted.net /dev/null'

  # Subset combined, sorted chain file.
  sh 'netChainSubset'                                                          \
     ' run/combined.chn.sorted.net run/combined.chn.sorted'                    \
     ' run/liftover.chn'
end


### Helpers ###

def add_to_PATH(path)
  return unless path
  return unless File.directory? path
  return if ENV['PATH'].split(':').include? path
  ENV['PATH'] = "#{path}:#{ENV['PATH']}"
end

def to_2bit(fas)
  sh "faToTwoBit #{fas} #{fas.ext('2bit')}"
end

def to_sizes(twobit)
  sh "twoBitInfo #{twobit} stdout | sort -k2nr > #{twobit.ext('sizes')}"
end

def extract_cdna(fas, gff)
  sh "gt extractfeat -type exon -join -retainids"  \
     " -seqfile #{fas} -matchdescstart"            \
     " #{gff} > #{gff.ext('.cdna.fa')}"
end

def extract_cds(fas, gff)
  sh "gt extractfeat -type CDS -join -retainids"   \
     " -seqfile #{fas} -matchdescstart"            \
     " #{gff} > #{gff.ext('.cds.fa')}"
end

def extract_pep(fas, gff)
  sh "gt extractfeat -type CDS -translate -join -retainids" \
     " -seqfile #{fas} -matchdescstart"                     \
     " #{gff} > #{gff.ext('.pep.fa')}"
end

def num_sequences(fas)
  `grep '>' #{fas} | wc -l`.strip
end

def num_exact(fas1, fas2)
  ENV['TMPDIR'] = Dir.pwd
  dir = Dir.mktmpdir
  system "grep -v '>' #{fas1} | sort > #{dir}/#{File.basename fas1}"
  system "grep -v '>' #{fas2} | sort > #{dir}/#{File.basename fas2}"
  `comm -12 #{dir}/#{File.basename fas1} \
   #{dir}/#{File.basename fas2} | wc -l`.strip
end

def summarize(outdir)
  inp_gff = "#{outdir}/input.gff"
  out_gff = "#{outdir}/lifted_cleaned.gff"

  extract_cdna('run/source.fa', inp_gff)
  extract_cds('run/source.fa', inp_gff)
  extract_pep('run/source.fa', inp_gff)

  extract_cdna('run/target.fa', out_gff)
  extract_cds('run/target.fa', out_gff)
  extract_pep('run/target.fa', out_gff)

  File.open("#{outdir}/summary.txt", 'w') do |file|
    %w(cdna cds pep).each do |tag|
      fas1 = inp_gff.ext(".#{tag}.fa")
      fas2 = out_gff.ext(".#{tag}.fa")
      next unless File.exist?(fas1) || File.exist?(fas2)

      file.puts tag.upcase
      file.puts "  source: #{num_sequences(fas1)}"
      file.puts "  target: #{num_sequences(fas2)}"
      file.puts "  exact:  #{num_exact(fas1, fas2)}"
    end
  end
end

def parallel(files, template)
  name = template.split.first
  jobs = files.map { |file| template % { :this => file } }
  joblst = "run/joblst.#{name}"
  joblog = "run/joblog.#{name}"
  File.write(joblst, jobs.join("\n"))
  sh "parallel --joblog #{joblog} -j #{jobs.length} -a #{joblst}"
end
