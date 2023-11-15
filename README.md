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
`metaDMG-cpp` requires `HTSlib`, a common library used for handling high-throughput sequencing data, `eigen3` and `gsl`.

On ubuntu these can be installed with:
```
sudo apt install libgsl-dev libeigen3-dev
```

### Installing and compiling metaDMG
To install metaDMG-cpp do:
```
git clone https://github.com/metaDMG-dev/metaDMG-cpp.git
cd metaDMG-cpp
make HTTSRC=../htslib
```

### To install HTSlib:
```
git submodule update --init --recursive
cd htslib
make
```

## Updating to latest version
For installing latest updates:
```
make clean
git pull https://github.com/metaDMG-dev/metaDMG-cpp.git
```

## Installing metaDMG using Conda
```
conda create -c bioconda -n metaDMG metadmg
conda activate metaDMG
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


# Damage analysis
```
./metaDMG-cpp getdamage [options] <in.bam>|<in.sam>|<in.cram>

Options:
  -n/--threads	        number of threads used for reading/writing (default: 4)
  -f/--fasta	        reference genome (required with CRAM)
  -l/--min_length	reads shorter than minlength will be discarded (default: 35)
  -p/--print_length	number of base pairs from read termini to estimate damage (default: 5)
  -r/--run_mode	        0: global estimate (default)
                        1: damage patterns will be calculated for each chr/scaffold contig.
  -i/--ignore_errors    continue analyses even if there are errors
  -o/--out_prefix	output prefix (default: meta)
```

# LCA analyses
```
./metaDMG-cpp lca [options]

Options:
  --threads             number of threads used for reading/writing (default: 4)
  --bam			input alignment file SAM/SAM.gz/BAM
  --names 		names.dmp.gz
  --nodes 		nodes.dmp.gz
  --acc2tax 		accesion to taxid table
  --edit_dist_min	minimum read edit distance
  --edit_dist_max	maximum read edit distance
  --min_mapq		minimum mapping quality
  --min_length		minimum read length
  --sim_score_low	number between 0-1
  --sim_score_high	number between 0-1
  --fix_ncbi		
  --discard
  --how_many		integer for many positions
  --lca_rank		family/genus/species
  --used_reads
  --no_rank2species
  --skip_no_rank
  --weight_type
  -i/--ignore_errors       continue analyses even if there are errors 1 or stop when error 0 (default)
  --temp                temp prefix
  -o/--out_prefix 		output prefix
```

# metaDMG print

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
 
# metaDMG merge 

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

# metaDMG mergedamage

# metaDMG index

