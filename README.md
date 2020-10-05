# metaDamage2.0
Fast and efficient method for error and damage rate estimates. Possible modes for running. Program is utilizing mdz field of the aux part of reads and is therefore reference free.

1. Basic single genome analysis with one overall global estimate. Similar to mapdamage1.0 and mapdamage2.0.

2. Basic eDNA or metagenomic analyses. Output is a damage estimate for each referenceID.

3. Lowest common ancestor is used for determining specificity of alignments for the analyses.

For all analyses output is a binary '.bdamage.gz' file, that can be accessed with the 'metadamage print' functionality.


# Single genome analysis
`./metadamage getdamage -l 10 -p 5 --threads 8 input.bam`

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

# Multiple genome analysis
` ./metadamage getdamage -l 10 -p 5 --threads 8 input.bam -r 1 `

All options found below:

```
./metadamage getdamage

Usage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>

Example: ./metadamage getdamage -l 10 -p 5 --threads 8 input.bam
Options:
  -f/--fasta	 is required with CRAM
  -l/--minlength	 reads shorter than minlength will be discarded
  -r/--runmode	runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.
		runmode 0 means one global estimate.
  -@/--threads	 Number of threads used for reading/writing
```

# LCA analyses
` ./metadamage lca -names names.dmp.gz -nodes nodes.dmp.gz -acc2tax taxid_accssionNO.gz -simscorelow 0.95 -simscorehigh 1.0 -minmapq 30  -bam input.bam ´

All options found below:

```
./metadamage lca 

Usage: metadamage lca -names -nodes -acc2tax [-editdist[min/max] -simscore[low/high] -minmapq -discard] -bam

Example ./metadamage lca -names names.dmp.gz -nodes nodes.dmp.gz -acc2tax taxid_accssionNO.gz -simscorelow 0.95 -simscorehigh 1.0 -minmapq 30  -bam input.bam
Options:
-simscorelow 0.95 
-simscorehigh 1.0 
-names /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/names.dmp.gz 
-nodes /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/nodes.dmp.gz 
-acc2tax /willerslev/users-shared/science-snm-willerslev-zsb202/soft/ngslca/accession2taxid/combined_taxid_accssionNO_20200425.gz 
-bam UE1210-Mex-59-Lib-7.col.file.sort.bam 
-outnames output_prefix
-lca_rank family/genus/species
-discard
-minmapq #integer for minimum mapping quality
-howmany #integer for many positions
```


# ./metadamage print

` ./metadamage print 
./metadamage print file.bdamage.gz [-names file.gz -bam file.bam -ctga -countout]
infile: (null) inbam: (null) names: (null) search: -1 ctga: 0 `

```
 - -ctga ONLY print CT+ and GA- (the damage ones)
 - -countout print mismatch as counts and not as transition
 probabilites
 - -r taxid Only print for specific taxid
 - -names NCBI names.dmp file - option that prints taxonomic names to
 output
 - -bam print referencenames from bamfile, otherwise it prints integeroffset. 
 ```

#### Header in print: taxid,nralign,orientation,position,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT 

# ./metadamage merge 

```
./metadamage merge 
./metadamage merge file.lca file.bdamage.gz [-names file.gz -bam file.bam -howmany 5 -nodes trestructure.gz]
```

- -howmany #integer for many positions
- -nodes #needs taxonomic paths to calculate damage higher than species level /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/nodes.dmp.gz
- -names #NCBI names.dmp file - option that prints taxonomic names to output  


./metadamage getdamage file.bam

./metadamage mergedamage files.damage.*.gz

./metadamage index files.damage.gz

