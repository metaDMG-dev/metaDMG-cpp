#pragma once
#include <zlib.h>
#include <htslib/sam.h>
#include <cstring>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <pthread.h>
#include "types.h"

typedef struct{
  //filenames
  char *htsfile;//bam,cram,sam
  char *nodesfile;
  char *namesfile;
  char *acc2taxfile;
  //hts strucutures
  samFile *hts;
  bam_hdr_t *header;
  //parameters for filtering reads
  double simscoreLow;
  double simscoreHigh;
  int editdistMin;
  int editdistMax;
  char *outnames;
  FILE *fp1;
  FILE *fp2;
  FILE *fp3;
  int minmapq;
  int discard; //or bitoperation with the flag of the read
  int minlength;
  char2int *charref2taxid;
  char *lca_rank;
  int norank2species;
}pars;

pars *get_pars(int argc,char **argv);
void print_pars(FILE *fp,pars *p);
void pars_free(pars *p);
int fexists(const char* str);
int fexists2(const char*str1,const char* str2);
int fexists3(const char*str1,const char* str2,const char *str3);
BGZF *getbgzf(const char*str1,const char *mode,int nthreads);
BGZF *getbgzf2(const char*str1,const char *str2,const char *mode,int nthreads);
BGZF *getbgzf3(const char*str1,const char *str2,const char *str3,const char *mode,int nthreads);

