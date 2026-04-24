// gpl thorfinn@binf.ku.dk
#include <time.h>         // for clock, time, clock_t, time_t

#include <cstdio>   // for fprintf, stderr
#include <cstring>  // for strcmp

#include "Aggregate_stat.h"  // for main_aggregate
#include "main_dfit.h"       // for main_dfit
#include "main_getdamage.h"
#include "main_pmd.h"        // for main_pmd
#include "main_print.h"
#include "version.h"         // for METADAMAGE_VERSION

#ifdef __REGRESSION__
#include "regression.h"  // for main_regression
#endif

// from ngsLCA.cpp
int main_lca(int argc, char **argv);
int main(int argc, char **argv) {
  fprintf(stderr,"\t-> metaDMG version: %s (htslib: %s) build(%s %s)\n",METADAMAGE_VERSION,hts_version(),__DATE__,__TIME__); 
    clock_t t = clock();
    time_t t2 = time(NULL);
    int rc = 0;

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
        rc = main_regression(argc, argv);
    else
#endif
    if (!strcmp(argv[0], "pmd"))
        rc = main_pmd(argc, argv);
    else if (!strcmp(argv[0], "getdamage"))
        rc = main_getdamage(argc, argv);
    else if (!strcmp(argv[0], "print"))
        rc = main_print(argc, argv);
    else if (!strcmp(argv[0], "print_all"))
        rc = main_print_all(argc, argv);
    else if (!strcmp(argv[0], "print_ugly"))
        rc = main_print_ugly(argc, argv);
    else if (!strcmp(argv[0], "dfit"))
      rc = main_dfit(argc, argv);
    else if (!strcmp(argv[0], "aggregate"))
      rc = main_aggregate(argc, argv);
    else if (!strcmp(argv[0], "print2"))
        rc = main_print2(argc, argv);
    else if (!strcmp(argv[0], "lca"))
        rc = main_lca(argc, argv);
    else{
      fprintf(stderr, "\t-> Unknown command: %s will exit error\n", argv[0]);
      rc = 1;
    }
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
    return rc;
}
