// gpl thorfinn@binf.ku.dk
#include <getopt.h>       // for optarg, getopt_long, optind, option
#include <htslib/bgzf.h>  // for bgzf_read, bgzf_close, bgzf_open, BGZF
#include <htslib/hts.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/sam.h>   // for sam_hdr_read, bam_hdr_destroy, sam_hdr_t...
#include <strings.h>      // for strcasecmp
#include <time.h>         // for clock, time, clock_t, time_t
#include <zlib.h>         // for gzprintf, gzclose, gzgets, gzopen, Z_NULL
#include <climits>

#include <cstdio>   // for fprintf, NULL, stderr, stdout, fopen
#include <cstdlib>  // for atoi, exit, free, atof, calloc, malloc
#include <cstring>  // for strdup, strcmp, strtok, strlen, strncpy
#include <map>      // for __map_iterator, map, operator!=, map<>::...
#include <vector>   // for vector

#include "main_pmd.h"    // for main_pmd
#include "main_dfit.h"   // for main_pmd
#include "Aggregate_stat.h"   // for main_pmd
#include "ngsLCA.h"      // for mean, var, gccontent, print_chain
#include "profile.h"     // for mydataD, mydata2, load_bdamage3, (anonym...
#include "shared.h"      // for parse_names, parse_nodes
#include "types.h"       // for int2intvec, int2int
#include "version.h"     // for METADAMAGE_VERSION
#include "main_print.h"

#ifdef __REGRESSION__
#include "regression.h"  // for main_regression
#endif


typedef std::map<int, char *> int2char;
int usage_getdamage(FILE *fp) {
    fprintf(fp, "\nUsage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>\n");
    fprintf(fp, "\nExample: ./metaDMG-cpp getdamage -l 10 -p 5 --threads 8 ../data/subs.sam\nOptions:\n");
    fprintf(fp, "  -n/--threads\t\t number of threads used for reading/writing (default: 4)\n");
    fprintf(fp, "  -f/--fasta\t\t reference genome (required with CRAM)\n");
    fprintf(fp, "  -l/--min_length\t minimum read length (default: 35)\n");
    fprintf(fp, "  -p/--print_length\t number of base pairs from read termini to estimate damage (default: 5)\n");
    fprintf(fp, "  -r/--run_mode\t\t 0: global estimate (default)\n\t\t\t 1: damage patterns will be calculated for each chr/scaffold contig\n");
    fprintf(fp, "  -i/--ignore_errors\t continue analyses even if there are errors.\n");
    fprintf(fp, "  -o/--out_prefix\t output prefix (default: meta)\n");
    fprintf(fp, "  -z/--rlens_flat_out\t make flat output of bins. Nice for computers\n");
    return 1;
}

int main_getdamage(int argc, char **argv) {
    if (argc == 1)
        return usage_getdamage(stderr);

    //  int MAXLENGTH = 256;
    int minLength = 35;
    int printLength = 5;
    char *refName = NULL;
    char *fname = NULL;
    int runmode = 0;  // this means one species, runmode=1 means multi species
    htsFile *fp = NULL;
    char *onam = strdup("meta");
    int nthreads = 4;
    int ignore_errors = 0;
    int rlens_flat_out = 0;
    htsFormat *dingding2 = (htsFormat *)calloc(1, sizeof(htsFormat));
    // fix thesepro
    static struct option lopts[] = {
        {"threads", required_argument, 0, 'n'},
	{"rlens_flat_out", required_argument, 0, 'z'},
        {"fasta", required_argument, 0, 'f'},
        {"min_length", required_argument, 0, 'l'},
        {"print_length", required_argument, 0, 'p'},
        {"run_mode", required_argument, 0, 'r'},
	{"ignore_errors", no_argument, &ignore_errors, 1},
        {"out_prefix", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}};

    int c;
    while ((c = getopt_long(argc, argv,
                            "iz:n:f:l:p:r:o:h",
                            lopts, NULL)) >= 0) {
      switch (c) {
      case 'i':
        ignore_errors = 1;
        break;
      case 'n':
	nthreads = atoi(optarg);
	break;
      case 'f':
	refName = strdup(optarg);
	break;
      case 'z':
	rlens_flat_out = atoi(optarg);
	break;
      case 'l':
	minLength = atoi(optarg);
	break;
      case 'p':
	printLength = atoi(optarg);
	break;
      case 'r':
	runmode = atoi(optarg);
	break;
      case 'o':
	free(onam);
	onam = strdup(optarg);
	break;
      case 'h':
	return usage_getdamage(stdout);
      default:
	fprintf(stderr, "Never here: %s\n", optarg);
	break;
      }
    }
    if (optind < argc)
        fname = strdup(argv[optind]);
    fprintf(stderr, "\t-> ./metaDMG-cpp refName: %s min_length: %d print_length: %d run_mode: %d out_prefix: %s nthreads: %d ignore_errors: %d rlens_flat_out: %d\n", refName, minLength, printLength, runmode, onam, nthreads, ignore_errors,rlens_flat_out);
    if (fname == NULL) {
        usage_getdamage(stderr);
        return 0;
    }
    if (refName) {
      int lenlen = 10 + strlen(refName) + 1;
      char *ref = (char *)malloc(lenlen);
      snprintf(ref,lenlen, "reference=%s", refName);
        hts_opt_add((hts_opt **)&dingding2->specific, ref);
        free(ref);
    }

    if ((fp = sam_open_format(fname, "r", dingding2)) == NULL) {
        fprintf(stderr, "[%s] nonexistant file: %s\n", __FUNCTION__, fname);
        exit(1);
    }

    bam1_t *b = bam_init1();
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if(hdr==NULL){
      fprintf(stderr,"\t-> Hello Doctor! im afraid there is an error reading the header\n");
      exit(1);
    }
    int checkIfSorted(char *str);
    if(checkIfSorted(hdr->text)){
      fprintf(stderr, "Input alignment file is not sorted.");
      if(!ignore_errors)
	return 1;
    }
    int ret;
    damage *dmg = new damage(printLength, nthreads, 0);
    int skipper[4] = {3, 3, 3, 3};
    std::map<int, std::vector<float> > gcconts;
    std::map<int, std::vector<float> > seqlens;
    while (((ret = sam_read1(fp, hdr, b))) >= 0) {
        if (bam_is_unmapped(b)) {
            if (skipper[0])
                fprintf(stderr, "skipping: %s unmapped, this msg is printed: %d times more\n", bam_get_qname(b), --skipper[0]);
            continue;
        }
        if (bam_is_failed(b)) {
            if (skipper[1])
                fprintf(stderr, "skipping: %s failed: flags=%d, this msg is printed: %d times more\n", bam_get_qname(b), b->core.flag, --skipper[1]);
            continue;
        }
        if (b->core.l_qseq < minLength) {
            if (skipper[2])
                fprintf(stderr, "skipping: %s too short, this msg is printed %d times more \n", bam_get_qname(b), --skipper[2]);
            continue;
        }

        dmg->damage_analysis(b, runmode != 0 ? b->core.tid : 0, 1);

        float mygc = gccontent(b);
        float mylen = b->core.l_qseq;

        int whichref = 0;
        if (runmode == 1)//<- runmode 0 means local one means global
            whichref = b->core.tid;
        std::map<int, std::vector<float> >::iterator it = gcconts.find(whichref);
        if (it == gcconts.end()) {
	  std::vector<float> tmp1;
	  tmp1.push_back(mygc);
	  gcconts[whichref] = tmp1;
	  std::vector<float> tmp2;
	  tmp2.push_back(mylen);
	  seqlens[whichref] = tmp2;
        } else {
	  if(it->second.size()<1000000)//<- we dont need more than a million values do we?
            it->second.push_back(mygc);
	  it = seqlens.find(whichref);
	  if (it == seqlens.end()) {
	    fprintf(stderr, "\t-> Error: iterator reached end in seqlens, will exit\n");
	    exit(1);
	  }
	  if(it->second.size()<1000000)
            it->second.push_back(mylen);
        }
    }

    dmg->printit(stdout, printLength);
    dmg->write(onam, runmode == 1 ? hdr : NULL);
    dmg->bwrite(onam,rlens_flat_out);

    // write stat
    char buf[1024];
    snprintf(buf, 1024, "%s.stat.gz", onam);
    fprintf(stderr, "\t-> Outputting overall statistic in file: \"%s\"\n", buf);

    gzFile fpstat = NULL;
    if((fpstat = gzopen(buf, "wb")) == NULL){
      fprintf(stderr,"\t-> Error problem opening file %s will exit\n",buf);
      exit(1);
    }
    gzprintf(fpstat,"taxid\tnreads\tmean_len\tvar_len\tmean_gc\tvar_gc\tlca\trank\n");
    for (std::map<int, std::vector<float> >::iterator it = gcconts.begin(); it != gcconts.end(); it++) {
        std::map<int, triple>::iterator it2 = dmg->assoc.find(it->first);
	if (it2 == dmg->assoc.end()) {
	  fprintf(stderr, "\t-> Error: iterator reached end in dmg->assoc, will exit\n");
	  exit(1);
	}
        std::map<int, std::vector<float> >::iterator it3 = seqlens.find(it->first);
	if (it3 == seqlens.end()) {
	  fprintf(stderr, "\t-> Error: iterator reached end in seqlens (it3), will exit\n");
	  exit(1);
	}
        if (0)
            gzprintf(fpstat, "%d\t%lu\t%f\t%f\t%f\t%f\tNA\tNA\n", it->first, it2->second.nreads, mean(it3->second), var(it3->second), mean(it->second), var(it->second));
        else{
	  if(runmode==1)
	    gzprintf(fpstat, "%s\t%lu\t%f\t%f\t%f\t%f\tNA\tNA\n", sam_hdr_tid2name(hdr, it->first), it2->second.nreads, mean(it3->second), var(it3->second), mean(it->second), var(it->second));
	  else
	    gzprintf(fpstat, "global\t%lu\t%f\t%f\t%f\t%f\tNA\tNA\n", it2->second.nreads, mean(it3->second), var(it3->second), mean(it->second), var(it->second));
	    
	}
    }
    gzclose(fpstat);

    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    sam_close(fp);
    destroy_damage(dmg);
    free(fname);
    free(onam);
    if(refName) free(refName);
    free(dingding2);
    return 0;
}
// from ngsLCA.cpp
int main_lca(int argc, char **argv);
int main(int argc, char **argv) {
  fprintf(stderr,"\t-> metaDMG version: %s (htslib: %s) build(%s %s)\n",METADAMAGE_VERSION,hts_version(),__DATE__,__TIME__); 
    clock_t t = clock();
    time_t t2 = time(NULL);

    if (argc == 1) {
#ifdef __REGRESSION__
        fprintf(stderr, "./metaDMG-cpp regression [other options]\n");
#endif
        fprintf(stderr, "./metaDMG-cpp pmd [other options]\n");
        fprintf(stderr, "./metaDMG-cpp getdamage file.bam\n");
        fprintf(stderr, "./metaDMG-cpp index files.damage.gz\n");
        fprintf(stderr, "./metaDMG-cpp lca [many options]\n");
        fprintf(stderr, "./metaDMG-cpp print bdamage.gz\n");
        fprintf(stderr, "./metaDMG-cpp print2 [many options] bdamage.gz\n");
        fprintf(stderr, "./metaDMG-cpp print_all [many options] bdamage.gz\n");
        fprintf(stderr, "./metaDMG-cpp print_ugly [many options] bdamage.gz\n");
	    fprintf(stderr, "./metaDMG-cpp dfit [many options] bdamage.gz\n");
        fprintf(stderr, "./metaDMG-cpp aggregate [many options] bdamage.gz\n");

        return 0;
    }
    fprintf(stderr, "\t-> ");
    for (int i = 0; i < argc; i++)
        fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, "\n");
    fflush(stderr);
    argc--;
    ++argv;
#ifdef __REGRESSION__
    if (!strcmp(argv[0], "regression"))
        return main_regression(argc, argv);
#endif
    if (!strcmp(argv[0], "pmd"))
        main_pmd(argc, argv);
    else if (!strcmp(argv[0], "getdamage"))
        main_getdamage(argc, argv);
    else if (!strcmp(argv[0], "print"))
        main_print(argc, argv);
    else if (!strcmp(argv[0], "print_all"))
        main_print_all(argc, argv);
    else if (!strcmp(argv[0], "print_ugly"))
        main_print_ugly(argc, argv);
    else if (!strcmp(argv[0], "dfit"))
      main_dfit(argc, argv);
    else if (!strcmp(argv[0], "aggregate"))
      main_aggregate(argc, argv);
    else if (!strcmp(argv[0], "print2"))
        main_print2(argc, argv);
    else if (!strcmp(argv[0], "lca"))
        main_lca(argc, argv);
    else{
      fprintf(stderr, "\t-> Unknown command: %s will exit error\n", argv[0]);
      return 1;
    }
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
    return 0;
}
