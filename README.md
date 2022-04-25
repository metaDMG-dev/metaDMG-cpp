# metaDMG-cpp
## [![tests](https://github.com/ANGSD/metadamage/actions/workflows/make.yml/badge.svg)](https://github.com/ANGSD/metadamage/actions/workflows/make.yml) 

# Introduction

metaDMG-cpp is a fast and efficient method for estimating mutation and damage rates in ancient DNA data. It relies on a commonly used alignment files formats bam/sam/sam.gz and can calculate degree of damage from read data mapped against a single as well as multiple genomes. It is especially relevant for data where data have been mapped against multiple reference genomes or to speed up analysis for damage estimation for a single genome.

Possible modes for running. Program is utilizing mdz field of the aux part of reads and is therefore reference free.

1. Basic single genome analysis with one overall global estimate. Similar to mapdamage1.0 and mapdamage2.0.

2. Basic eDNA or metagenomic (e.g. multiple genome) analyses. Output is a damage estimate per reference, taxonomic name or accession no.

3. Integrating a Least Common Ancestor algorithm gives the opportunity to retrieve    specificity of alignments for the analyses.

For all analyses output is a binary '.bdamage.gz' file, that can be accessed with the 'metadamage print' functionality.

# Installation

### Dependencies
metaDMG-cpp requires HTSlib - a common library used for handling high-throughput sequencing data.
eigen3 and gsl.

On ubuntu these can be installed with:
```
sudo apt install libgsl-dev libeigen3-dev
```

To install HTSlib do:
```
git clone https://github.com/SAMtools/htslib
cd htslib
make
```
### Installing and compiling metaDMG
To install metaDMG-cpp do:
```
git clone https://github.com/metaDMG-dev/metaDMG-cpp.git
cd metadamage
make HTTSRC=../htslib
```
```

## Updating to latest version
For installing latest updates in the directory metaDMG-cpp do:
```
make clean
git pull https://github.com/metaDMG-dev/metaDMG-cpp.git
make HTTSRC=../htslib
```

# Taxonomic resource files
metaDMG-cpp can perform damage estimations on internal nodes within a taxonomy (e.g. species, genus and family level). In order to traverse up a taxonomic tree the program needs three files in NCBI taxonomy format. These can either be a custom taxonomy or simply rely on the NCBI taxonomy.
``` 
# Downloading resource files for program from NCBI
mkdir ncbi_tax_dmp;
cd ncbi_tax_dmp/;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip;
unzip new_taxdump.zip;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz;
gunzip nucl_gb.accession2taxid.gz;
```


# Single genome analysis
`./metaDMG-cpp getdamage -l 10 -p 5 --threads 8 input.bam`

All options found below:

```
./metaDMG-cpp getdamage

Usage: ./metaDMG-cpp getdamage [options] <in.bam>|<in.sam>|<in.cram>

Example: ./metaDMG-cpp getdamage -l 10 -p 5 --threads 8 ../data/subs.sam
Options:
  -f/--fasta	 is required with CRAM
  -l/--minlength	 reads shorter than minlength will be discarded
  -p/--printlength	number of base pairs from read termini to estimate damage (default: 5)
  -o/--outname	output prefix
  -r/--runmode	runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.
		runmode 0 means one global estimate.
  -@/--threads	 Number of threads used for reading/writing
```

# Multiple genome analysis
` ./metaDMG-cpp getdamage -l 10 -p 5 --threads 8 input.bam -r 1 `

All options found below:

```
./metaDMG-cpp getdamage

Usage: ./metaDMG-cpp getdamage [options] <in.bam>|<in.sam>|<in.cram>

Example: ./metaDMG-cpp getdamage -l 10 -p 5 --threads 8 input.bam
Options:
  -f/--fasta	 is required with CRAM
  -l/--minlength	 reads shorter than minlength will be discarded
  -p/--printlength	number of base pairs from read termini to estimate damage (default: 5)
  -o/--outname	output prefix
  -r/--runmode	runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.
		runmode 0 means one global estimate.
  -@/--threads	 Number of threads used for reading/writing
```

# LCA analyses

` ./metaDMG-cpp lca -names names.dmp.gz -nodes nodes.dmp.gz -acc2tax taxid_accssionNO.gz -simscorelow 0.95 -simscorehigh 1.0 -minmapq 30 -howmany 15 -bam input.bam`

All options found below:

```
./metaDMG-cpp lca 

Usage: ./metaDMG-cpp lca [options] -names -nodes -acc2tax [-editdist[min/max] -simscore[low/high] -minmapq -discard] -bam <in.bam>|<in.sam>|<in.sam.gz>

Example ./metaDMG-cpp lca -names names.dmp.gz -nodes nodes.dmp.gz -acc2tax taxid_accssionNO.gz -simscorelow 0.95 -simscorehigh 1.0 -minmapq 30 -howmany 15 -bam input.bam

Options:
  -simscorelow	number between 0-1
  -simscorehigh	number between 0-1
  -names 	/willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/names.dmp.gz 
  -nodes 	/willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/nodes.dmp.gz 
  -acc2tax 	/willerslev/users-shared/science-snm-willerslev-zsb202/soft/ngslca/accession2taxid/combined_taxid_accssionNO_20200425.gz 
  -bam		input alignment file in bam, sam or sam.gz format
  -outnames	output prefix
  -lca_rank	family/genus/species 
  -discard 
  -minmapq	integer for minimum mapping quality
  -howmany	integer for many positions
```


# ./metaDMG-cpp print

`./metaDMG-cpp print file.bdamage.gz [-names file.gz -bam file.bam -ctga -countout] infile: (null) inbam: (null) names: (null) search: -1 ctga: 0 `

All options found below:

```
 ./metaDMG-cpp print 
 
 Usage: ./metaDMG-cpp print 
 
Example ./metaDMG-cpp print file.bdamage.gz -names names.dmp.gz 
	./metaDMG-cpp print file.bdamage.gz -r 9639 -ctga
	./metaDMG-cpp print file.bdamage.gz -countout 
Options:
  -ctga		ONLY print CT+ and GA- (the damage ones)
  -countout	print mismatch as counts and not as transition probabilites
  -r taxid	Only print for specific taxid
  -names	NCBI names.dmp.gz file - option to print taxonomic names instead of NCBI TaxID
  -bam		print referencenames (Accession No.) from bamfile, otherwise it prints NCBI TaxId. 


#### Header in print: taxid,nralign,orientation,position,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT 
 ```
 
# ./metaDMG-cpp merge 

`./metaDMG-cpp merge file.lca file.bdamage.gz ` 


All options found below:

```
./metaDMG-cpp merge 

Usage: ./metaDMG-cpp merge file.lca file.bdamage.gz [-names names.dmp.gz -bam <in.bam>|<in.sam>|<in.sam.gz> -howmany 5 -nodes nodes.dmp.gz]

Example
Options:
-howmany #integer for many positions ?? Thorfinn is this also working for the merge module? Or only at the LCA module? 
-nodes #needs taxonomic paths to calculate damage higher than species level /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/nodes.dmp.gz
-names #NCBI names.dmp file - option that prints taxonomic names to output  
```

# ./metaDMG-cpp mergedamage files.damage.*.gz

# ./metaDMG-cpp index files.damage.gz

