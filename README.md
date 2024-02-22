# metaDMG-cpp
## [![tests](https://github.com/metaDMG-dev/metaDMG-cpp/actions/workflows/make.yml/badge.svg)](https://github.com/metaDMG-dev/metaDMG-cpp/actions/workflows/make.yml) 

# Introduction

metaDMG-cpp is a fast and efficient method for estimating mutation and damage rates in ancient DNA data, especially from ancient metagenomes. It relies on commonly used alignment file formats bam/sam/sam.gz and can calculate the degree of damage from read data mapped against a single as well as multiple genomes. It is especially relevant for data where data have been mapped against multiple reference genomes or to speed up analysis for damage estimation for a single genome.

There are three possible modes for running metaDMG. In all cases, the program utilises the MD:Z field of the aux part of reads and is therefore reference-free.

1. Basic single genome analysis with one overall global estimate. Similar to the mapdamage2.0 programme (./metaDMG-cpp getdamage --run_mode 0).

2. Basic metagenomes (e.g. multiple genome alignments) analyses. Output is a damage estimate per reference, taxonomic name or accession no that includes all alignments without an LCA analysis (./metaDMG-cpp getdamage --run_mode 1).

3. Integrating a Least Common Ancestor algorithm (ngsLCA) which allows retrieving damage estimates for alignments classified to the given taxonomic levels (./metaDMG-cpp lca).

For all analyses, the output is a binary '.bdamage.gz' file, which contains a substitution matrix that can be accessed and read with the (./metaDMG print)' functionality.

# Installation

### Dependencies
`metaDMG-cpp` requires `HTSlib`, a common library used for handling high-throughput sequencing data, `eigen3` and `gsl`.

On Ubuntu these can be installed with:
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
./metaDMG-cpp lca counts substitutions between read and reference on internal nodes within a taxonomy (e.g. species, genus and family level). To traverse up a taxonomic tree the program needs three files in NCBI taxonomy format. These can either be a custom taxonomy built as the NCBI taxonomy or simply rely on the NCBI taxonomy, or it can  be a combination. NOTE the taxonomy file shall reflect the version of the database you are using. 

**Downloading resource files for the program from NCBI**
``` 
mkdir ncbi_tax_dmp;
cd ncbi_tax_dmp/;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip;
unzip new_taxdump.zip;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz;
gunzip nucl_gb.accession2taxid.gz;
```


# Damage analysis
Calculations of substitution matrices (without any LCA analysis), either at a:
	- global mode (--run_mode 0) e.g. one matrix for the whole alignment file. Which is useful for single taxa analysis or if you want a global estimate for a metagenome. 
 	- local mode  (--run_mode 1) e.g. a matrix for each reference with alignments. Which can be useful for microbial analysis and simulation of ancient metagenomes. 
```
./metaDMG-cpp getdamage [options] <in.bam>|<in.sam>|<in.cram>

Options:
  -n/--threads	        number of threads used for reading/writing (default: 4)
  -f/--fasta	        reference genome (required with CRAM)
  -l/--min_length	reads shorter than minlength will be discarded (default: 35)
  -p/--print_length	number of positions along the read termini that are used to estimate the damage (default: 5)
  -r/--run_mode	        0: **global**  (default)
                        1: **local** damage patterns will be calculated for each chr/scaffold contig.
  -i/--ignore_errors    continue analyses even if there are errors
  -o/--out_prefix	output prefix (default: meta)
```

# LCA analyses
Calculations, where each read, is classified to a given taxonomic level (based on the similarity of the alignments as specified below, we recommend --sim_score_low 0.95 --sim_score_high 1.0).  
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
  --how_many		integer for number of positions OBS rasmus will change this to --print_length 
  --lca_rank		such as family/genus/species
  --used_reads
  --no_rank2species
  --skip_no_rank
  --weight_type
  -i/--ignore_errors       continue analyses even if there are errors 1 or stop when error 0 (default)
  --temp                temp prefix
  -o/--out_prefix 		output prefix
```


# metaDMG aggregate  
Aggregating the lca statistics when transversing through the tree structure, creating files with prefix .aggregate.stat.txt.gz

`./metaDMG-cpp aggregate file.bdamage.gz --lcastat lca.stat --names names.dmp --nodes nodes.dmp --out file` 

```
	-> ./metaDMG-cpp aggregate -h 
    Aggregation of lca produced statistics (mean length, variance length, mean GC, variance GC) when transversing up the nodes of the tree structure
		./metaDMG-cpp aggregate file.bdamage.gz --names file.gz --nodes trestructure.gz --lcastat file.stat --out filename
--help 		 Print extended help page to see all options.
--names 	 names.dmp.gz
--nodes 	 nodes.dmp.gz
--lca 		 lcaout.stat lca produced statistics
--out 		 Suffix of outputname with the predetermined prefix (.aggregate.stat.txt.gz)
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
-howmany #integer for many positions.
-nodes #needs taxonomic paths to calculate damage higher than species level /ncbi_tax_dump_files/nodes.dmp.gz
-names #NCBI names.dmp file - option that prints taxonomic names to output  
```

 
# metaDMG dfit  
Performing numerical optimization of the deamination frequencies based on the mismatch matrix (.bdamage.gz) to estimate four parameters: A,q,c,phi. Either applying a beta-binomial distribution or binomial distribution as the choice for likelihood model.

A: Amplitude of damage on position one.

q: Relative decrease of damage per position.

c: Background noise, equivalent to sequencing errors.

phi: signifies the uncertainty of beta-binomial model, with larger values indicating the probability of deamination is reduced to binomial distribution.

The optimization can be performed in three modes, all of which depends on the input parameters and information stored within the .bdamage.gz file, 

global: with one damage estimate for entire BAM file. 

local: damage estimate for all chromosomes/contigs/scaffolds present in BAM header.

lca: damage estimates of the nodes within the lca tree structure.

`./metaDMG-cpp dfit file.bdamage.gz ` 

A simple help page is printed with 

`./metaDMG-cpp dfit -h ` 

While an extended and full helppage is printed with --help and shown below.

```
	-> ./metaDMG-cpp dfit --help 
Extended helppage damage estimation using numerical optimization
Estimates damage in either local, global or lca mode depending on the bdamage input format
Dfit command performs either optimization of beta-binomial (default --nbootstrap = 0) or binomial model (--nbootstrap > 2)

Local mode: 	  damage estimated for each reference/chromosome in the BAM file
Global mode: 	  one damage estimate for whole BAM file
lca mode: 	  damage estimated over the lca tree at different ranks


--help 				 Print extended help page to see all options.

./metaDMG-cpp dfit file.bdamage.gz --names file.gz --nodes trestructure.gz --lcastat fil.gz --bam file.bam --showfits int --nopt int --nbootstrap int --seed int --doCI int --CI float --lib <ds,ss> --out file

------------ Required ------------- 
./metaDMG-cpp dfit file.bdamage.gz 	 bdamage file contains the misincorporation matrix, in global mode from getdamage command or local mode from lca command

------------ Optional ------------- 

--nopt 				 Number of optimization calls (default: 10).

--seed 				 Seed value for random number generators (default: computer time).

--out 				 Prefix for output name.

--lib 				 double stranded (ds) use C>T (forward) and G>A (reverse); single stranded (ss) use C>T for both forward and reverse (default ds)

--nbootstrap 			 number of bootstrap iterations. default: 1 -> use Beta-binomial model, -nbootstrap >1 use Binomial model 
--bam 				In local mode - convert the internal id numbering from bdamage.gz to the reference in the bam header
--nodes       #needs taxonomic paths to calculate damage higher than species level
--names       #NCBI names.dmp file - option that prints taxonomic names to output  
--acc2tax 		accesion to taxid table

---- Optimization model specific ---- 

------  Beta-binomial model ------

--showfits: 			 Verbose parameter, stored in the .dfit.txt.gz, default = 0, see documentation for full description of columns
	--showfits 0		 [id;A;q;c;phi;llh;ncall;sigmaD;Zfit] 
				 id:contig; A:amplitute of damage; q:decrease in damage; c:offset (background noise); phi:variance between Beta-binomial and binomial model
				 llh:minimized negative log-likelihood; ncall:optimization calls; sigmaD: standard deviation; Zfit: significance score

	--showfits 1		 [id;A;q;c;phi;llh;ncall;sigmaD;Zfit;
				 fwdxi;fwdxConfi;..;fwdxn;fwdxConfn..;
					 bwdx;dxConfn;..;bwdxn;bwdxConfn]
				 fwdx: damage estimates forward strand; fwdxConf: confidence interval; fbwdx: reverse strand; bwdxConfn: confidence interval.
 					 i position 1 to position n (zero-index)

	--showfits 2		 [id;A;q;c;phi;llh;ncall;sigmaD;Zfit;
				 fwKi;fwNi;fwdxi;fwf0;fwdxConfi;..;fwKi;fwNi;fwdxi;fwfi;fwdxConfi..;
					 bwKi;bwNi;bwdxi;bwfi;bwdxConfi;..;bwKn;bwNn;bwdxn;fwfn;bwdxConfn]
				 fwKi: Number of transitions forward strand (C>T); fwNi: Number reference counts forward strand (C); fwfi: C>T frequency based on forward strand from bdamage file
				 bwKi: Number of transitions reverse strand (G>A); bwNi: Number reference counts reverse strand (G); bwfi: C>T frequency based on reverse strand from bdamage file


---------- Binomial model ----------

--showfits: 			 Verbose parameter, stored in the .dfit.txt.gz, default = 0
	--showfits 0		 [id;A;q;c;phi;llh;ncall;sigmaD;Zfit;A_b;q_b;c_b;phi_b;A_CI_l;A_CI_h;q_CI_l;q_CI_h;c_CI_l;c_CI_h;phi_CI_l;phi_CI_h] 
					 id;A;q;c;phi;llh;ncall;sigmaD;Zfit estimates from beta-binomial (see above)
					 A_b;q_b;c_b;phi_b;A_CI_l;A_CI_h;q_CI_l;q_CI_h;c_CI_l;c_CI_h;phi_CI_l;phi_CI_h estimates binomial model with bootstrap method (*_b), confidence interval (*_CI_l) and  (*_CI_h)

	--showfits 1		 similar columns added as described above, using the bootstrap estimated parameters;

	--showfits 2		 similar columns added as described above, using the bootstrap estimated parameters;
--nbootstrap: 			 Number of bootstrap iterations, default = 0
--doboot: 			 Store all bootstrap iterations of [id,A_b,q_b,c_b,phi_b] in seperate file, suffix: .bdamage.gz.boot.stat.txt.gz
--sigtype 			 Determines the bootstrap method <1,2,3>, default = 1
	0: 				 Sample with replacement from the C>T frequency from the .bdamage format
	1: 				 Sample with replacement from the K and N column from the .bdamage format and calculates f column
	1: 				 Sample with replacement from binomial distribution defined from K and N column from the .bdamage format
--CI: 				 Desired confidence interval, default = 0.95
--doCI: 			 Confidence interval type <0,1>
	0: 				 Mean and Standard Deviation-Based CI Calculation
	1: 				 Percentile CI Calculation

---------- Examples ----------
	 Local mode binomial:
 		 ./metaDMG-cpp getdamage Pitch6.bam -l 10 -p 15 -r 1 -o Pitch6getDMG
 		 ./metaDMG-cpp dfit Pitch6getDMG.bdamage.gz --showfits 1
	 Global mode binomial:
 		 ./metaDMG-cpp getdamage Pitch6.bam -l 10 -p 15 -r 0 -o Pitch6getDMG
 		 ./metaDMG-cpp dfit Pitch6getDMG.bdamage.gz --nbootstrap 2 --showfits 2
	 Lca mode beta-binomial:
 		 ./metaDMG-cpp lca --names names.dmp --nodes nodes.dmp --acc2tax acc2taxid.map.gz --weight_type 1 --fix_ncbi 0 --bam Pitch6.bam --out Pitch6lcatest 
 		 ./metaDMG-cpp dfit Pitch6lcatest.bdamage.gz --names names.dmp --nodes nodes.dmp --showfits 0
```
