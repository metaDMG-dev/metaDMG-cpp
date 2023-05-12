// g++ compressbam.cpp -O3 -o compressbam ../../htslib/libhts.a -lz -lbz2 -llzma -lpthread -lcurl -ggdb                                                                                                                                 
// ./compressbam -hts small.sam -out small.out.bam2 -type bam -@ 16

#include <cstdio>
#include <zlib.h>
#include <map>
#include <cstring>
#include <cstdlib>
#include <iostream>//for printing time
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/bgzf.h>
#include <cassert>
#include <ctype.h>
#include <htslib/kstring.h>
#include <htslib/hfile.h>

//#include "../shared.h"

#define VERSION "0.1"

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));
htsThreadPool p = {NULL, 0};
int nthreads = 8;
char out_mode[5]="ws";

//get ref used not get refused
size_t getrefused(samFile *htsfp,bam_hdr_t *hdr,int *keeplist,int &nkeep){
   //now mainloop
  bam1_t *aln = bam_init1();

  size_t at =0;
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    if((at++ %100000)==0  )
      fprintf(stderr,"\r\t->  Now at read:     %lu ",at);
    int chr = aln->core.tid ; //contig name (chromosome)
    assert(chr!=-1);
    keeplist[chr] = keeplist[chr]+1;
    chr =aln->core.mtid; 
    if(chr!=-1)
      keeplist[chr] = keeplist[chr]+1;
  }
  for(int i=0;i<sam_hdr_nref(hdr);i++)
    if(keeplist[i]!=-1)
      nkeep++;
  fprintf(stderr,"\n\t-> Done reading list of refids: to keep %d\n",nkeep);
  bam_destroy1(aln);
  return at;
}

void writemod(const char *outfile ,bam_hdr_t *hdr,int *keeplist,samFile *htsfp,char *mycl){
  BGZF *fp = NULL;
  fp=bgzf_open(outfile,"w5");//write bgzf with compression level5. Maybe this should be changed.
  bgzf_mt(fp,nthreads,256);
  assert(fp!=NULL);
  
  kstring_t kstmp ={0,0,NULL};
  kstring_t newhdr_ks ={0,0,NULL};


  sam_hdr_find_hd(hdr,&kstmp);
  kputs(kstmp.s,&newhdr_ks);
  kputc('\n',&newhdr_ks);
  //  for(int2int::iterator it=keep.begin();it!=keep.end();it++){
  for(int i=0;i<sam_hdr_nref(hdr);i++){
    if(keeplist[i]==-1)
      continue;
    //    fprintf(stderr,"SQ: it->first: %d\n",it->first);
    kstmp.l =0;
    assert(sam_hdr_find_line_pos(hdr,"SQ",i,&kstmp)==0);
    kputs(kstmp.s,&newhdr_ks);
    kputc('\n',&newhdr_ks);
    if(newhdr_ks.l>10000000){
      assert(bgzf_write(fp,newhdr_ks.s,newhdr_ks.l)==newhdr_ks.l);
      newhdr_ks.l =0;
    }
    
  }

  for(int i=0;i<sam_hdr_count_lines(hdr,"RG");i++){
    kstmp.l =0;
    assert(sam_hdr_find_line_pos(hdr,"RG",i,&kstmp)==0);
    kputs(kstmp.s,&newhdr_ks);
    kputc('\n',&newhdr_ks);
    if(newhdr_ks.l>10000000){
      assert(bgzf_write(fp,newhdr_ks.s,newhdr_ks.l)==newhdr_ks.l);
      newhdr_ks.l =0;
    }
  }
  for(int i=0;i<sam_hdr_count_lines(hdr,"PG");i++){
    kstmp.l =0;
    assert(sam_hdr_find_line_pos(hdr,"PG",i,&kstmp)==0);
    kputs(kstmp.s,&newhdr_ks);
    kputc('\n',&newhdr_ks);
  }
  for(int i=0;i<sam_hdr_count_lines(hdr,"CO");i++){
    kstmp.l =0;
    assert(sam_hdr_find_line_pos(hdr,"CO",i,&kstmp)==0);
    kputs(kstmp.s,&newhdr_ks);
    kputc('\n',&newhdr_ks);
  }


  assert(bgzf_write(fp,newhdr_ks.s,newhdr_ks.l)==newhdr_ks.l);
  free(newhdr_ks.s);
  free(kstmp.s);
  bgzf_close(fp);
  fprintf(stderr,"\t-> Done writing new header as text\n");fflush(stderr);

  samFile *tmpfp = NULL;
  tmpfp = hts_open(outfile,"r");
  assert(tmpfp!=NULL);
  sam_hdr_t *newhdr = NULL;
  newhdr = sam_hdr_read(tmpfp);
  assert(newhdr!=NULL);
  sam_close(tmpfp);

  //fprintf(stderr,"newhdr: %s\n",sam_hdr_str(newhdr));
  fprintf(stderr,"\t-> header info before: nref: %d bytesize: %lu\n",sam_hdr_nref(hdr),sam_hdr_length(hdr));
  fprintf(stderr,"\t-> header info  after: nref: %d bytesize: %lu\n",sam_hdr_nref(newhdr),sam_hdr_length(newhdr));
  fprintf(stderr,"\t-> header info reduction nref: %f bytesize: %f\n",(double)sam_hdr_nref(newhdr)/sam_hdr_nref(hdr),(double)sam_hdr_length(newhdr)/sam_hdr_length(hdr));
  fflush(stderr);
  //write outputfile
  samFile *outhts = NULL;
  if ((outhts = sam_open_format(outfile,out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing: %s\n",outfile);
    exit(0);
  }
  if(nthreads>1){
    if (!(p.pool = hts_tpool_init(nthreads))) {
      fprintf(stderr, "Error creating thread pool\n");
      exit(0);
    }
    hts_set_opt(outhts,  HTS_OPT_THREAD_POOL, &p);
  }
  for(int i=0;i<sam_hdr_nref(hdr);i++){
    if(keeplist[i]==-1)
      continue;

    keeplist[i] = sam_hdr_name2tid(newhdr,sam_hdr_tid2name(hdr,i));
    assert(keeplist[i]>=0);
  }
  sam_hdr_add_pg(newhdr,"compressbam","VN",VERSION,"CL",mycl,NULL);

  int ret = sam_hdr_write(outhts,newhdr);
  fprintf(stderr,"\t-> Done writing new header as binary\n");fflush(stderr);
  //now mainloop
  bam1_t *aln = bam_init1();
  size_t at =1;
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    if((at++ %100000)==0  ){
      fprintf(stderr,"\r\t->  Now at read:     %lu ",at);
      sam_flush(htsfp);
      fflush(stderr);
    }

    aln->core.tid = keeplist[aln->core.tid];
    if(aln->core.mtid!=-1)
      aln->core.mtid = keeplist[aln->core.mtid];
    assert(sam_write1(outhts, newhdr,aln)>=0);
  }
  sam_flush(htsfp);
  sam_hdr_destroy(newhdr);
  bam_destroy1(aln);
  sam_close(outhts);
  hts_tpool_destroy(p.pool);

}

int main(int argc,char**argv){

  //print time
  clock_t t=clock();
  time_t t2=time(NULL);

  if(argc==1){
    fprintf(stderr,"./compressbam -hts -ref -out -type [Sam,Cram,Bam] -@ nthreads\n");
    return 0;
  }
  char *mycl = stringify_argv(argc,argv);
  fprintf(stderr,"\t-> compressbam: (%s;%s;%s): \'%s\'\n",__FILE__,__DATE__,__TIME__,mycl);  
  
  argv++;
  char *hts = NULL;
  char *names = NULL;
  char *ref = NULL;
  char *outfile = strdup("tmp.sam");

  //  char out_mode[5]="ws";
  //  int nthreads = 4; //this is now global
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    //fprintf(stderr,"key: %s val: %s\n",key,val);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-out",key))
      {
	free(outfile);
	outfile=strdup(val);
      }
    else if(!strcasecmp("-ref",key)) ref=strdup(val);
    else if(!strcasecmp("-@",key)) nthreads=atoi(val);
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
  
  fprintf(stderr,"\t-> hts: %s out: %s type: %s ref: %s nthreads: %d\n",hts,outfile,out_mode,ref,nthreads);
  
  //open inputfile and parse header
  samFile *htsfp = hts_open(hts,"r");
  bam_hdr_t *hdr = sam_hdr_read(htsfp);
  //  int64_t record_begin = htell(htsfp);
  //  int64_t ret = hts_tell_func(htsfp->fp); //this should work at some piont
  fprintf(stderr,"\t-> Header has now been read. Will now start list of refIDs to use\n");fflush(stderr);
  int nkeep = 0;
  int *keeplist = new int[sam_hdr_nref(hdr)];
  for(int i=0;i<sam_hdr_nref(hdr);i++)
    keeplist[i] = -1;
  size_t nalign = getrefused(htsfp,hdr,keeplist,nkeep);
  fprintf(stderr,"\t-> Number of alignments parsed: %lu\n",nalign);
  sam_close(htsfp);

  sam_hdr_destroy(hdr);
  htsfp=hts_open(hts,"r" );
  hdr = sam_hdr_read(htsfp);
  writemod(outfile,hdr,keeplist,htsfp,mycl);
  fprintf(stderr,"\t-> Done writing file: \'%s\'\n",outfile);
  free(mycl);
  sam_close(htsfp);

  if(hts) free(hts);
  if(outfile) free(outfile);
  if(ref) free(ref);

  free(dingding2);
  sam_hdr_destroy(hdr);
  delete [] keeplist;

  fprintf(stderr, "\t-> [ALL done] cpu-time used =  %.2f sec walltime used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));
  return 0;
}
