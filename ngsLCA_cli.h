#pragma once
#include <htslib/bgzf.h>  // for BGZF
#include <htslib/sam.h>   // for bam_hdr_t, samFile
#include <stdio.h>        // for FILE
#include <zlib.h>         // for gzFile

#include "types.h"  // for char2int

typedef struct {
    // filenames
    char *htsfile;  // bam,cram,sam
    char *nodesfile;
    char *namesfile;
    char *acc2taxfile;
    char *filteredAcc2taxfile;
    // hts strucutures
    samFile *hts;
    bam_hdr_t *header;
    // parameters for filtering reads
    double simscoreLow;
    double simscoreHigh;
    int editdistMin;
    int editdistMax;
    int skipnorank;
    char *outnames;
    gzFile fp1;
    gzFile fp_lcadist;
    gzFile fp2;
  //    FILE *fp3; //this is the logfile that fgv thinkgs sholld be removed
    int minmapq;
    int discard;  // or bitoperation with the flag of the read
    int minlength;
    char2int *charref2taxid;
    char *lca_rank;
    int norank2species;
    int howmany;
    char *usedreads_sam;
    int fixdb;  // used for disabling mod_db function
    int nthreads;
    int weighttype;
    char *tempfolder;
    int ignore_errors;
  int reallyDump;
  long maxreads;
} pars;

pars *get_pars(int argc, char **argv);
void print_pars(FILE *fp, pars *p);
void pars_free(pars *p);
int fexists(const char *str);
int fexists2(const char *str1, const char *str2);
int fexists3(const char *str1, const char *str2, const char *str3);
BGZF *getbgzf(const char *str1, const char *mode, int nthreads);
BGZF *getbgzf2(const char *str1, const char *str2, const char *mode, int nthreads);
BGZF *getbgzf3(const char *str1, const char *str2, const char *str3, const char *mode, int nthreads);
