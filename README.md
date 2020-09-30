# metadamage

/willerslev/scratch/mwp/metaDamage_test/metadamage/metadamage lca
 -names
 -nodes
 -acc2tax
 [-editdist[min/max] -simscore[low/high] -minmapq -discard]
 -bam
 -lca_rank #for damage estimate


/willerslev/scratch/mwp/metaDamage_test/metadamage/metadamage print file.bdamage.gz -bam /-names
-ctga #printer kun relevante mutationer
-r taxid #printer kun bestemte taxid
-names #NCBI names.dmp file - option that prints taxonomic names to output /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/names.dmp.gz
-bam file.bam #optional give it bamfile after 
#### Header in print: taxid nralign orientation position AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT 



/willerslev/scratch/mwp/metaDamage_test/metadamage/metadamage merge file.lca file.bdamage.gz [-names file.gz -bam file.bam -howmany 5 -nodes trestructure.gz]
-howmany #integer for many positions
-nodes #needs taxonomic paths to calculate damage higher than species level /willerslev/users-shared/science-snm-willerslev-npl206/ngsLCA/ngsLCA/ncbi_tax_dump_files/nodes.dmp.gz
-names #NCBI names.dmp file - option that prints taxonomic names to output  


./metadamage getdamage file.bam

./metadamage mergedamage files.damage.*.gz

./metadamage index files.damage.gz

