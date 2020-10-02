# metadamage
Fast and efficient method for error and damage rate analysis. Possible modes for running. Program is utilizing mdz field of the aux part of reads and is therefore reference free.

1. Basic single genome analysis with one overall global estimate. Similar to mapdamage1.0 and mapdamage2.0.

2. Basic eDNA or metagenomic analyses. Output is a damage estimate for each referenceID.

3. Lowest common ancestor is used for determining specificity of alignments for the analyses.

For all analyses output is a binary '.bdamage.gz' file, that can be accessed with the 'metadamage print' functionality.


# Single genome analysis
`./metadamage getdamage -l 10 -p 5 --threads 8 ../data/subs.bam`

All options found below:

```
./metadamage getdamage

Usage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>

Example: ./metadamage getdamage -l 10 -p 5 --threads 8 ../data/subs.sam
Options:
  -f/--fasta	 is required with CRAM
  -l/--minlength	 reads shorter than minlength will be discarded
  -r/--runmode	runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.
		runmode 0 means one global estimate.
  -@/--threads	 Number of threads used for reading/writing
```

# Basic eDNA analysis
`./metadamage getdamage -l 10 -p 5 --threads 8 ../data/subs.bam -r 1`

All options found below:

```
./metadamage getdamage

Usage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>

Example: ./metadamage getdamage -l 10 -p 5 --threads 8 ../data/subs.sam
Options:
  -f/--fasta	 is required with CRAM
  -l/--minlength	 reads shorter than minlength will be discarded
  -r/--runmode	runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.
		runmode 0 means one global estimate.
  -@/--threads	 Number of threads used for reading/writing
```

# LCA analyses



./metadamage lca
 -names
 -nodes
 -acc2tax
 [-editdist[min/max] -simscore[low/high] -minmapq -discard]
 -bam
 -lca_rank #for damage estimate


./metadamage print file.bdamage.gz -bam /-names
-ctga #printer kun relevante mutationer
-r taxid #printer kun bestemte taxid
-names #NCBI names.dmp file - option that prints taxonomic names to output /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/names.dmp.gz
-bam file.bam #optional give it bamfile after 
#### Header in print: taxid nralign orientation position AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT 



./metadamage merge file.lca file.bdamage.gz [-names file.gz -bam file.bam -howmany 5 -nodes trestructure.gz]
-howmany #integer for many positions
-nodes #needs taxonomic paths to calculate damage higher than species level /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/nodes.dmp.gz
-names #NCBI names.dmp file - option that prints taxonomic names to output  


./metadamage getdamage file.bam

./metadamage mergedamage files.damage.*.gz

./metadamage index files.damage.gz

