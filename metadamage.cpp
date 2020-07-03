 //gpl thorfinn@binf.ku.dk
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctype.h>
htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

size_t ***fw = NULL;
size_t ***rv = NULL;

size_t ***init(int length){
  size_t ***ret=(size_t ***) new size_t**[length];
  for(int at=0;at<length;at++){
    ret[at] = new size_t*[5];
    for(int i=0;i<5;i++){
      ret[at][i] = new size_t[5];
      for(int j=0;j<5;j++)
	ret[at][i][j] =0;
    }
  }
  return ret;
}

#define LENS 10

void merge_md_cig(int ncig,uint32_t *cigs,uint8_t *mdz){
  fprintf(stderr,"%s\n",mdz);
  if(*mdz) mdz++;//skip z
  fprintf(stderr,"%s\n",mdz);
  uint32_t ret[1024];
  int val_predel =0;
  int at=0;
  while(*mdz){
    fprintf(stderr,"mdzbegin: \'%c\'\n",*mdz);
    //skip deletion part
    if(*mdz=='^'){
      fprintf(stderr,"Deletion part\n");
      while(!isdigit(*mdz))
	mdz++;
      fprintf(stderr,"Done Deletion delpart\n");
      continue;
    }
    
    int val=0;
    while(isdigit(*mdz)){
      val = val_predel+ val*10+(int)*mdz-'0';
      val_predel =0;
      mdz++;
    }
    fprintf(stderr,"Value part after loop: %d\n",val);

    char ch = *mdz;
    if(ch!='^'){
      ret[at++] =0;
      continue;
    }else
      val_predel = 0;
    fprintf(stderr,"nextchar: mdz: %c\n",ch);
    mdz++;
  }
  //  for(int i=0;i<)
  fprintf(stderr,"nop in mdz:%d\n",at);
  exit(0);
  
}



int main(int argc, char **argv){
  char *refName = NULL;
  char *fname = argv[1];
  htsFile *fp = NULL;
  fw = init(LENS);
  rv = init(LENS);
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
  while(((ret=sam_read1(fp,hdr,b)))>0){
    uint8_t *mdz = bam_aux_get(b, "MD");

    uint32_t *cigs = bam_get_cigar(b);
    int nCig = b->core.n_cigar;

    fprintf(stdout,"ret:%d b.core->tid:%d pos:%d mdz:%s flags:%d qnam:%s\n",ret,b->core.tid,(int)b->core.pos,(char *)mdz,b->core.flag,(char *)bam_get_qname(b));
    merge_md_cig(nCig,cigs,mdz);
    return 0;
    for(int i=0;0&&i<nCig;i++) {
      int opCode = cigs[i]&BAM_CIGAR_MASK; //what to do
      int opLen = cigs[i]>>BAM_CIGAR_SHIFT; //length of what to do
      fprintf(stderr,"opCode=%d opLen=%d\n",opCode,opLen);
      if(opCode==BAM_CINS||opCode==BAM_CDEL) 
	fprintf(stderr,"INS OR DEL\n");
      else if(opCode==BAM_CSOFT_CLIP)
	fprintf(stderr,"SOFTCLIP\n");
      else if(opCode==BAM_CMATCH||opCode==BAM_CEQUAL||opCode==BAM_CDIFF) 
	fprintf(stderr,"MATCH EQUAL OR DIFF\n");
      else if(opCode==BAM_CREF_SKIP) 
	fprintf(stderr,"REF SKIP\n");
      else if(opCode==BAM_CPAD||opCode==BAM_CHARD_CLIP)
	fprintf(stderr,"PAD OR HARDCLIP\n");
      else
	fprintf(stderr,"Problem with unsupported CIGAR opCode=%d\n",opCode);
    }
  }
  return 0;
}

