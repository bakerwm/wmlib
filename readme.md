


## STAR

+ Get the fasta file from STAR index

Go to the index directory of STAR:

```
cat Genome | tr '\000\001\002\003\004' 'ACGTN' | tr -s '\005' '\n' | awk 'BEGIN {while (getline < "chrName.txt") {ii++;c[ii]=$1} } {if (NR>length(c)) exit; print ">" c[NR]; print}' > Genome.fasta

```
