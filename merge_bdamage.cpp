#include <cstdio>

#include "profile.h"


int main_mergedamage(int argc, char **argv){

  for(int i=0;i<argc;i++){
    fprintf(stderr,"argc: %d val: %s\n",i,argv[i]);
  }

  return 0;
}
