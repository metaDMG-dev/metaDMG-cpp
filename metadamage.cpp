 //gpl thorfinn@binf.ku.dk
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctype.h>

#include "profile.h"

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

int main(int argc, char **argv){
  int MAXLENGTH = 1000;
  int minLength = 30;
  char *refName = NULL;
  char *fname = argv[1];
  htsFile *fp = NULL;

  if(argc>2)
    refName = argv[2];
  if(refName){
    char *ref =(char*) malloc(10 + strlen(refName) + 1);
    sprintf(ref, "reference=%s", refName);
    hts_opt_add((hts_opt **)&dingding2->specific,ref);
    free(ref);
  }
  if((fp=sam_open_format(fname,"r",dingding2))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(0);
  }
  
  bam1_t *b = bam_init1();
  bam_hdr_t  *hdr = sam_hdr_read(fp);
  int ret;
  damage *dmg = init_damage(MAXLENGTH);
  while(((ret=sam_read1(fp,hdr,b)))>0){
      if(bam_is_unmapped(b) ){
      	fprintf(stderr,"skipping: %s unmapped \n",bam_get_qname(b));
      continue;
    }
    if(bam_is_failed(b) ){
      fprintf(stderr,"skipping: %s failed: flags=%d \n",bam_get_qname(b),b->core.flag);
      continue;
    }
    if(b->core.l_qseq < minLength){
      fprintf(stderr,"skipping: %s too short \n");
      continue;
    }
    if(bam_is_paired(b)){
      fprintf(stderr,"skipping: %s  is paired (can be considered using the -paired flag\n",bam_get_qname(b));
      continue;
    }
    dmg->damage_analysis(b,0);
  }


  dmg->printit(stdout,5);
  dmg->write(NULL,NULL);
    
  sam_hdr_destroy(hdr);
  bam_destroy1(b);
  sam_close(fp);
  destroy_damage(dmg);
  return 0;
}

