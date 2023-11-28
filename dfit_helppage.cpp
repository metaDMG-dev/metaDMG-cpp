#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <math.h>
#include "dfit_helppage.h"

int HelpPageSimple(FILE *fp){
  fprintf(fp,"Simple helppage damage estimation using numerical optimization\n");
  fprintf(fp,"Estimates damage in either local, global or lca mode depending on the bdamage input format\n");
  fprintf(fp,"Dfit command performs either optimization of beta-binomial (default --nbootstrap = 0) or binomial model (--nbootstrap > 2)\n\n");

  fprintf(fp,"Local mode: \t  damage estimated for each reference/chromosome in the BAM file\n");
  fprintf(fp,"Global mode: \t  one damage estimate for whole BAM file\n");
  fprintf(fp,"lca mode: \t  damage estimated over the lca tree at different ranks\n\n");

  fprintf(fp,"\n--help \t\t\t\t Print extended help page to see all options.\n\n");

  fprintf(stderr, "./metaDMG-cpp dfit file.bdamage.gz --names file.gz --nodes trestructure.gz --lcastat fil.gz --bam file.bam --showfits int --nopt int --nbootstrap int --seed int --doCI int --CI float --lib <ds,ss> --out file\n");

  fprintf(fp,"\n------------ Required ------------- \n");
  fprintf(fp,"file.bdamage.gz contains the misincorporation matrix, in global mode from getdamage command either global or local mode or from lca command\n");
  fprintf(fp,"e..g ./metaDMG-cpp dfit file.bdamage.gz\n");

  fprintf(fp,"\n------------ Optional ------------- \n");
  fprintf(fp,"--showfits: \t\t\t Verbose parameter, stored in the .dfit.txt.gz, default = 0, see documentation or extendend helppage --help, for full description of columns\n");
  fprintf(fp,"--nopt \t\t\t\t Number of optimization calls (default: 10).\n");
  fprintf(fp,"--lib \t\t\t\t double stranded (ds) use C>T (forward) and G>A (reverse); single stranded (ss) use C>T for both forward and reverse (default ds)\n");
  fprintf(fp,"--nbootstrap: \t\t\t Number of bootstrap iterations, default = 0, i.e. optimization of beta-binomial model\n");
  fprintf(fp,"--nthreads: \t\t\t Number of threads, default = 1, i.e. no threading\n");

  fprintf(stderr, "\n---------- Examples ----------\n");
  fprintf(stderr, "\t Local mode binomial:\n \t\t ./metaDMG-cpp getdamage Pitch6.bam -l 10 -p 15 -r 1 -o Pitch6getDMG\n \t\t ./metaDMG-cpp dfit Pitch6getDMG.bdamage.gz --showfits 1\n");
  fprintf(stderr, "\t Global mode binomial:\n \t\t ./metaDMG-cpp getdamage Pitch6.bam -l 10 -p 15 -r 0 -o Pitch6getDMG\n \t\t ./metaDMG-cpp dfit Pitch6getDMG.bdamage.gz --nbootstrap 2 --showfits 2\n");
  fprintf(stderr, "\t Lca mode beta-binomial:\n \t\t ./metaDMG-cpp lca --names names.dmp --nodes nodes.dmp --acc2tax acc2taxid.map.gz --weight_type 1 --fix_ncbi 0 --bam Pitch6.bam --out Pitch6lcatest \n \t\t ./metaDMG-cpp dfit Pitch6lcatest.bdamage.gz --names names.dmp --nodes nodes.dmp --lcastat Pitch6lcatest.stat --showfits 0\n");

  exit(1);
  return 0;
}

int HelpPage(FILE *fp){
  fprintf(fp,"Extended helppage damage estimation using numerical optimization\n");
  fprintf(fp,"Estimates damage in either local, global or lca mode depending on the bdamage input format\n");
  fprintf(fp,"Dfit command performs either optimization of beta-binomial (default --nbootstrap = 0) or binomial model (--nbootstrap > 2)\n\n");

  fprintf(fp,"Local mode: \t  damage estimated for each reference/chromosome in the BAM file\n");
  fprintf(fp,"Global mode: \t  one damage estimate for whole BAM file\n");
  fprintf(fp,"lca mode: \t  damage estimated over the lca tree at different ranks\n\n");

  fprintf(fp,"\n--help \t\t\t\t Print extended help page to see all options.\n\n");

  fprintf(stderr, "./metaDMG-cpp dfit file.bdamage.gz --names file.gz --nodes trestructure.gz --lcastat fil.gz --bam file.bam --showfits int --nopt int --nbootstrap int --seed int --doCI int --CI float --lib <ds,ss> --out file\n");
  
  fprintf(fp,"\n------------ Required ------------- \n");
  fprintf(fp,"./metaDMG-cpp dfit file.bdamage.gz \t bdamage file contains the misincorporation matrix, in global mode from getdamage command or local mode from lca command\n");
  fprintf(fp,"\n------------ Optional ------------- \n");
  fprintf(fp,"\n--nopt \t\t\t\t Number of optimization calls (default: 10).\n");
  fprintf(fp,"\n--seed \t\t\t\t Seed value for random number generators (default: computer time).\n");
  fprintf(fp,"\n--out \t\t\t\t Prefix for output name.\n");
  fprintf(fp,"\n--lib \t\t\t\t double stranded (ds) use C>T (forward) and G>A (reverse); single stranded (ss) use C>T for both forward and reverse (default ds)\n");
  fprintf(fp,"\n--nbootstrap \t\t\t number of bootstrap iterations. default: 1 -> use Beta-binomial model, -nbootstrap >1 use Binomial model ");
  fprintf(fp,"\n--bam \t\t\t\t In local mode - convert the internal id numbering from bdamage.gz to the reference in the bam header\n");
  fprintf(fp,"\n-rng | --rand: \t\t\t Pseudo-random number generator, OS specific\n");
  fprintf(fp,"\t\t <0,1,2,3> \n"); 
  fprintf(fp,"\t\t 0 :  \t\t\t drand48_r, default for linux or unix, not available for MacOS.\n"); 
  fprintf(fp,"\t\t 1 :  \t\t\t std::uniform_int_distribution\n"); 
  fprintf(fp,"\t\t 2 :  \t\t\t rand_r\n"); 
  fprintf(fp,"\t\t 3 :  \t\t\t erand48, default for MacOS.\n"); 

  fprintf(fp,"\n\n---- Optimization model specific ---- \n");
  fprintf(stderr, "\n\n------  Beta-binomial model ------\n\n");
  fprintf(fp,"--showfits: \t\t\t Verbose parameter, stored in the .dfit.txt.gz, default = 0, see documentation for full description of columns\n");
  fprintf(fp,"\t--showfits 0\t\t [id;A;q;c;phi;llh;ncall;sigmaD;Zfit] \n");
  fprintf(fp,"\t\t\t\t id:contig; A:amplitute of damage; q:decrease in damage; c:offset (background noise); phi:variance between Beta-binomial and binomial model\n");
  fprintf(fp,"\t\t\t\t llh:minimized negative log-likelihood; ncall:optimization calls; sigmaD: standard deviation; Zfit: significance score\n");
  //showfits 1
  fprintf(fp,"\n\t--showfits 1\t\t [id;A;q;c;phi;llh;ncall;sigmaD;Zfit;\n");
  fprintf(fp,"\t\t\t\t fwdxi;fwdxConfi;..;fwdxn;fwdxConfn..;\n\t\t\t\t\t bwdx;dxConfn;..;bwdxn;bwdxConfn]\n");
  fprintf(fp,"\t\t\t\t fwdx: damage estimates forward strand; fwdxConf: confidence interval; fbwdx: reverse strand; bwdxConfn: confidence interval.\n \t\t\t\t\t i position 1 to position n (zero-index)\n");
  //showfits 2
  fprintf(fp,"\n\t--showfits 2\t\t [id;A;q;c;phi;llh;ncall;sigmaD;Zfit;\n");
  fprintf(fp,"\t\t\t\t fwKi;fwNi;fwdxi;fwf0;fwdxConfi;..;fwKi;fwNi;fwdxi;fwfi;fwdxConfi..;\n\t\t\t\t\t bwKi;bwNi;bwdxi;bwfi;bwdxConfi;..;bwKn;bwNn;bwdxn;fwfn;bwdxConfn]\n");
  fprintf(fp,"\t\t\t\t fwKi: Number of transitions forward strand (C>T); fwNi: Number reference counts forward strand (C); fwfi: C>T frequency based on forward strand from bdamage file\n");
  fprintf(fp,"\t\t\t\t bwKi: Number of transitions reverse strand (G>A); bwNi: Number reference counts reverse strand (G); bwfi: C>T frequency based on reverse strand from bdamage file\n");
  
  fprintf(stderr, "\n\n---------- Binomial model ----------\n\n");
  fprintf(fp,"--showfits: \t\t\t Verbose parameter, stored in the .dfit.txt.gz, default = 0\n");
  fprintf(fp,"\t--showfits 0\t\t [id;A;q;c;phi;llh;ncall;sigmaD;Zfit;A_b;q_b;c_b;phi_b;A_CI_l;A_CI_h;q_CI_l;q_CI_h;c_CI_l;c_CI_h;phi_CI_l;phi_CI_h] \n");
  fprintf(fp,"\t\t\t\t\t id;A;q;c;phi;llh;ncall;sigmaD;Zfit estimates from beta-binomial (see above)\n");
  fprintf(fp,"\t\t\t\t\t A_b;q_b;c_b;phi_b;A_CI_l;A_CI_h;q_CI_l;q_CI_h;c_CI_l;c_CI_h;phi_CI_l;phi_CI_h estimates binomial model with bootstrap method (*_b), confidence interval (*_CI_l) and  (*_CI_h)\n");
  fprintf(fp,"\n\t--showfits 1\t\t similar columns added as described above, using the bootstrap estimated parameters;\n");
  fprintf(fp,"\n\t--showfits 2\t\t similar columns added as described above, using the bootstrap estimated parameters;\n");
  fprintf(fp,"--nbootstrap: \t\t\t Number of bootstrap iterations, default = 0\n");
  fprintf(fp,"--doboot: \t\t\t Store all bootstrap iterations of [id,A_b,q_b,c_b,phi_b] in seperate file, suffix: .bdamage.gz.boot.stat.txt.gz\n");
  fprintf(fp,"--sigtype \t\t\t Determines the bootstrap method <1,2,3>, default = 1\n");
  fprintf(fp,"\t0: \t\t\t\t Sample with replacement from the C>T frequency from the .bdamage format\n");
  fprintf(fp,"\t1: \t\t\t\t Sample with replacement from the K and N column from the .bdamage format and calculates f column\n");
  fprintf(fp,"\t1: \t\t\t\t Sample with replacement from binomial distribution defined from K and N column from the .bdamage format\n");

  fprintf(fp,"--CI: \t\t\t\t Desired confidence interval, default = 0.95\n");
  fprintf(fp,"--doCI: \t\t\t Confidence interval type <0,1>\n");
  fprintf(fp,"\t0: \t\t\t\t Mean and Standard Deviation-Based CI Calculation\n");
  fprintf(fp,"\t1: \t\t\t\t Percentile CI Calculation\n");
  fprintf(fp,"--showfits: \t\t\t Verbose parameter, stored in the .dfit.txt.gz, default = 0\n");

  fprintf(stderr, "\n---------- Examples ----------\n");
  fprintf(stderr, "\t Local mode binomial:\n \t\t ./metaDMG-cpp getdamage Pitch6.bam -l 10 -p 15 -r 1 -o Pitch6getDMG\n \t\t ./metaDMG-cpp dfit Pitch6getDMG.bdamage.gz --showfits 1\n");
  fprintf(stderr, "\t Global mode binomial:\n \t\t ./metaDMG-cpp getdamage Pitch6.bam -l 10 -p 15 -r 0 -o Pitch6getDMG\n \t\t ./metaDMG-cpp dfit Pitch6getDMG.bdamage.gz --nbootstrap 2 --showfits 2\n");
  fprintf(stderr, "\t Lca mode beta-binomial:\n \t\t ./metaDMG-cpp lca --names names.dmp --nodes nodes.dmp --acc2tax acc2taxid.map.gz --weight_type 1 --fix_ncbi 0 --bam Pitch6.bam --out Pitch6lcatest \n \t\t ./metaDMG-cpp dfit Pitch6lcatest.bdamage.gz --names names.dmp --nodes nodes.dmp --lcastat Pitch6lcatest.stat --showfits 0\n");

    
  
  exit(1);
  return 0;
}
/*
    fprintf(stderr, "./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz -bam file.bam -showfits int -nopt int -nbootstrap int -seed int -doCI int -CI float -lib <ds,ss> -out file\n");
    fprintf(stderr, "Estimate damage patterns with beta-binomial model\n");
    fprintf(stderr, "\tEstimate damage patterns for each chr/scaffold contig (local mode), using lca stats\n");
    fprintf(stderr, "\t\t./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz -bam file.bam -showfits int -nopt int -nbootstrap int -seed int -doCI int -CI float -lib <ds,ss> -out file\n");
    fprintf(stderr, "\tEstimate one global damage pattern \n");
    fprintf(stderr, "\t\t./metaDMG-cpp dfit metaDMG-cpp/metaDMG-cpp dfit Pitch6getDMG.bdamage.gz -doboot 1 -nbootstrap 5 -nopt 10 -showfits 0\n");
    fprintf(stderr, "Estimate damage patterns with binomial model\n");
    fprintf(stderr, "\tEstimate damage patterns for each chr/scaffold contig (local mode), using lca stats\n");
    fprintf(stderr, "\t\t./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz -bam file.bam -showfits int -nopt int -nbootstrap int -seed int -doCI int -CI float -lib <ds,ss> -out file\n");
    fprintf(stderr, "\tEstimate one global damage pattern \n");
    fprintf(stderr, "\t\t./metaDMG-cpp dfit metaDMG-cpp/metaDMG-cpp dfit Pitch6getDMG.bdamage.gz -doboot 1 -nbootstrap 5 -nopt 10 -showfits 0\n");
*/
