#include <cstdio>
#include <zlib.h>
#include <map>
#include <cstring>
#include <cstdlib>
#include <htslib/sam.h>
#include <cassert>
#include <ctype.h>
#include "../shared.h"

int2int getrefused(samFile *htsfp,bam_hdr_t *hdr){
   //now mainloop
  bam1_t *aln = bam_init1();
  int2int hit;
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    int chr = aln->core.tid ; //contig name (chromosome)
    assert(chr!=-1);
    int2int::iterator it = hit.find(chr);
    if(it==hit.end())
      hit[chr] = 1;
    else
      it->second = it->second+1;
    
    chr =aln->core.mtid; 
    if(chr!=-1){
      int2int::iterator it = hit.find(chr);
      if(it==hit.end())
	hit[chr] = 1;
      else
	it->second = it->second+1;
    }
  }
  fprintf(stderr,"\t-> Done reading list of refids to keep: %lu\n",hit.size());
  bam_destroy1(aln);
  return hit;
}

void writemod(samFile *htsfp,bam_hdr_t *hdr,int2int &keep,samFile *outhts){
  sam_hdr_t *newhdr = sam_hdr_dup(hdr);
  fprintf(stderr,"sizeof header before: %d headerlen: %d\n",sam_hdr_nref(newhdr),sam_hdr_length(newhdr));
  for(int i=0;i<sam_hdr_nref(hdr);i++){
    int2int::iterator it = keep.find(i);
    if(it==keep.end())
      sam_hdr_remove_line_id(newhdr, "SQ", "SN", sam_hdr_tid2name(hdr,i)); 

  }
  fprintf(stderr,"sizeof header before: %d headerlen: %d\n",sam_hdr_nref(newhdr),sam_hdr_length(newhdr));

  for(int2int::iterator it=keep.begin();it!=keep.end();it++){
    //    fprintf(stderr,"pre: %d ",it->first);
    it->second = sam_hdr_name2tid(newhdr,sam_hdr_tid2name(hdr,it->first));
    //fprintf(stderr,"post: %d \n",it->second);
  }
  
  int ret = sam_hdr_write(outhts,newhdr);
  //now mainloop
  bam1_t *aln = bam_init1();
  int at =1;
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    int oldchr = aln->core.tid ; //contig name (chromosome)
    aln->core.tid = keep[oldchr];
    if(aln->core.mtid!=-1)
      aln->core.mtid = keep[aln->core.mtid];
    assert(sam_write1(outhts, newhdr,aln)>=0);
  }
  
  sam_hdr_destroy(newhdr);
  bam_destroy1(aln);
}

int main(int argc,char**argv){
  
  if(argc==1){
    fprintf(stderr,"./compressbam -hts -ref -out -type [Sam,Cram,Bam]\n");
    return 0;
  }
  argv++;
  char *hts = NULL;
  char *names = NULL;
  char *ref = NULL;
  char *outfile = "tmp.sam";

  char out_mode[5]="ws";
  
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    //fprintf(stderr,"key: %s val: %s\n",key,val);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-out",key)) outfile=strdup(val);
    else if(!strcasecmp("-ref",key)) ref=strdup(val);
    else if(!strcasecmp("-type",key)){
      char c = tolower(val[0]);
      if(c!='b'&&c!='s'&&c!='c'){
	fprintf(stderr,"\t-> Problem understanding output format, bam,cram or sam val: %s val[0]:%c\n",val,val[0]);
	return 0;
      }
      out_mode[1] = c;
    }
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      return 0;
    }
    ++argv;
  }
  
  fprintf(stderr,"\t-> hts: %s out: %s type: %s ref: %s\n",hts,outfile,out_mode,ref);

  //open inputfile and parse header
  samFile *htsfp = hts_open(hts,"r");
  bam_hdr_t *hdr = sam_hdr_read(htsfp); 

  htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));
  samFile *outhts = NULL;
  if ((outhts = sam_open_format(outfile,out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing: %s\n",outfile);
    exit(0);
  }

  int2int usedref = getrefused(htsfp,hdr);
  sam_close(htsfp);sam_hdr_destroy(hdr);
  htsfp=hts_open(hts,"r" );
  hdr = sam_hdr_read(htsfp); 
  writemod(htsfp,hdr,usedref,outhts);
  sam_close(htsfp);

  if(hts) free(hts);
  if(outfile) free(outfile);
  if(ref) free(ref);
  sam_close(outhts);
  free(dingding2);
  sam_hdr_destroy(hdr);

  return 0;
}
