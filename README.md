# metaDMG-cpp
## [![tests](https://github.com/metaDMG-dev/metaDMG-cpp/actions/workflows/make.yml/badge.svg)](https://github.com/metaDMG-dev/metaDMG-cpp/actions/workflows/make.yml) 

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
cd metaDMG-cpp
make HTTSRC=../htslib
```

## Updating to latest version
For installing latest updates in the directory metaDMG-cpp do:
```
make clean
git pull https://github.com/metaDMG-dev/metaDMG-cpp.git
```

## Installing metaDMG using Conda
```
conda create -n metaDMG htslib eigen cxx-compiler c-compiler gsl
conda activate metaDMG
git clone https://github.com/metaDMG-dev/metaDMG-cpp.git
cd metaDMG
make clean && make -j 8 HTSSRC=systemwide
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
`./metaDMG-cpp getdamage --threads 8 -l 10 -p 5 input.bam`

All options found below:

```
./metaDMG-cpp getdamage

Usage: ./metaDMG-cpp getdamage [options] <in.bam>|<in.sam>|<in.cram>

Example: ./metaDMG-cpp getdamage --threads 8 --minlength 10 --printlength 5 subs.sam
Options:
  -n/--threads	 Number of threads used for reading/writing
  -f/--fasta	 is required with CRAM
  -l/--minlength	 reads shorter than minlength will be discarded
  -p/--printlength	number of base pairs from read termini to estimate damage (default: 5)
  -o/--outname	output prefix
  -r/--runmode	runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.
		runmode 0 means one global estimate.
```

# Multiple genome analysis
`./metaDMG-cpp getdamage --threads 8 --minlength 10 --printlength 5 input.bam`

All options found below:

```
./metaDMG-cpp getdamage

Usage: ./metaDMG-cpp getdamage [options] <in.bam>|<in.sam>|<in.cram>

Example: ./metaDMG-cpp getdamage --threads 8 --minlength 10 --printlength 5 input.bam
Options:
  -n/--threads		Number of threads used for reading/writing
  -f/--fasta		is required with CRAM
  -l/--minlength	reads shorter than minlength will be discarded
  -p/--printlength	number of base pairs from read termini to estimate damage (default: 5)
  -o/--outname		output prefix
  -r/--runmode		runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.
			runmode 0 means one global estimate.
```

# LCA analyses

` ./metaDMG-cpp lca --threads 4 --bam input.bam --names names.dmp.gz --nodes nodes.dmp.gz --acc2tax taxid_accssionNO.gz --sim_score_low 0.95 --sim_score_high 1.0 --min_mapq 30 --how_many 15`

All options found below:

```
./metaDMG-cpp lca 

Usage: ./metaDMG-cpp lca --bam <in.bam>|<in.sam>|<in.sam.gz> --names --nodes --acc2tax [--edit_dist_[min/max] --sim_score_[low/high] --min_mapq --discard]

Example ./metaDMG-cpp lca --bam input.bam --names names.dmp.gz --nodes nodes.dmp.gz --acc2tax taxid_accssionNO.gz --sim_score_low 0.95 --sim_score_high 1.0 --min_mapq 30 --how_many 15

Options:
  --bam			input alignment file in bam, sam or sam.gz format
  --names 		names.dmp.gz 
  --nodes 		nodes.dmp.gz 
  --acc2tax 		combined_taxid_accssionNO_20200425.gz 
  --sim_score_low	number between 0-1
  --sim_score_high	number between 0-1
  --lca_rank		family/genus/species 
  --discard 
  --min_mapq		integer for minimum mapping quality
  --how_many		integer for many positions
  --out 		output prefix
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

