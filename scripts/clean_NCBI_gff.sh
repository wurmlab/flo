cat $1 | ruby -pe 'gsub(/end_range=\d+,\.;/, "")' | ruby -pe 'gsub(/start_range=\.,\d+;?/, "")' | ruby -pe 'gsub(/Gap=.+/, "")' | grep -v '##sequence-region' | grep -v '##species' > $2

# gt gff3 -tidy -sort -addids -retainids gnomon.gff > annot.gff
