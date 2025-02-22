#include <cstdio>

#include "profile.h"


int main_mergedamage(int argc, char **argv){
  fprintf(stderr,"argc: %d\n",argc);
  argc--;argv++;
  for(int i=0;i<argc;i++){
    fprintf(stderr,"argc: %d val: %s\n",i,argv[i]);
    int howmany;
    char *fname = argv[i];
    std::map<int, mydataD> retmap = load_bdamage_full(fname, howmany);
    fprintf(stderr, "\t-> Number of entries in damage pattern file: %lu printlength(howmany):%d\n", retmap.size(), howmany);
  }

  return 0;
}
