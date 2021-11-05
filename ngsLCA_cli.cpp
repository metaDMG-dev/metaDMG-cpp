#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <cassert>
#include <htslib/bgzf.h>
#include "ngsLCA_cli.h"
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
pars *pars_init(){
  pars *p =(pars*) calloc(1,sizeof(pars));
  p->htsfile = strdup("CHL_155_12485.sort.bam");
  p->acc2taxfile=strdup("nucl_gb.accession2taxid.gz");
  p->namesfile = strdup("names.dmp.gz");
  p->nodesfile= strdup("nodes.dmp.gz");
  p->hts=NULL;
  p->header=NULL;
  p->editdistMin=0;
  p->editdistMax=10;
  p->simscoreLow=0;
  p->simscoreHigh=1;
  p->fp1=p->fp2=p->fp3=p->fp_lcadist=NULL;
  p->outnames=strdup("outnames");
  p->minmapq=0;
  p->discard=516;//discard unmapped and read fail
  p->minlength=-1;
  p->charref2taxid = NULL;
  p->lca_rank = strdup("species");
  p->norank2species = 0;
  p->skipnorank = 1;
  p->howmany = 5;
  p->usedreads_sam = NULL;
  p->fixdb = 1;
  p->nthreads = 4;
  p->weighttype = 0;
  p->tempfolder = strdup("");
  return p;
}

void pars_free(pars *p){
  fclose(p->fp1);
  //fclose(p->fp2);
  fclose(p->fp3);

  if(p->header)
    sam_hdr_destroy(p->header);
  free(p->outnames);
  free(p->lca_rank);
  free(p->htsfile);
  free(p->acc2taxfile);
  free(p->namesfile);
  free(p->nodesfile);
  free(p->tempfolder);
  free(p);
}

void *read_header_thread(void *ptr){
  time_t t=time(NULL);
  pars *p = (pars *) ptr;
  fprintf(stderr,"\t-> [thread1] Will read header\n");
  p->hts = hts_open(p->htsfile,"r");
  p->header = sam_hdr_read(p->hts);
  assert(p->header);
  fprintf(stderr,"\t-> [thread1] Done reading header: %.2f sec, header contains: %d \n",(float)(time(NULL) - t),p->header->n_targets);
  pthread_mutex_unlock(&mutex1);
  return NULL;
}


char2int * ass2bin(const char *fname,int redo){
  const char *CONSTNAME = "delmeme.bin";
  gzFile FP=Z_NULL;
  char2int *cm = new char2int;
  //redo=1;
  //load binary representation
  if(redo==0&&fexists(CONSTNAME)){
    time_t t=time(NULL);
    fprintf(stderr,"\t-> [thread2] reading binary represenation\n");
    FP = gzopen(CONSTNAME,"rb");
    int key_l;
    
    while(sizeof(int)==gzread(FP,&key_l,sizeof(int))) {
      char *key =(char*) calloc(key_l+1,sizeof(char));
      assert(key_l=gzread(FP,key,key_l));
      int val;
      assert(sizeof(int)==gzread(FP,&val,sizeof(int)));
      if(cm->find(key)!=cm->end()){
	fprintf(stderr,"\t-> Duplicate entries found \'%s\'\n",key);
      }else
	(*cm)[key] = val;
    }
    fprintf(stderr,"\t-> [thread2] Done reading binary representation: %.2f sec\n",(float)(time(NULL) - t));

  }else{
    FP=gzopen(CONSTNAME,"wb");
  
    gzFile gz= Z_NULL;
    gz=gzopen(fname,"rb");
    if(gz==Z_NULL){
      fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
      exit(0);
    }
    char buf[4096];
    int at=0;
    char buf2[4096];
    extern int SIG_COND;
    gzgets(gz,buf,4096);//skip header
    while(SIG_COND&&gzgets(gz,buf,4096)){
      if(!((at++ %100000 ) ))
	if(isatty(fileno(stderr)))
	  fprintf(stderr,"\r\t-> At linenr: %d in \'%s\'      ",at,fname);
      strcpy(buf2,buf);
      strtok(buf,"\t\n ");
      char *key =strtok(NULL,"\t\n ");
      int val = atoi(strtok(NULL,"\t\n "));
      if(FP){
	int key_l = strlen(key);
	gzwrite(FP,&key_l,sizeof(int));
	gzwrite(FP,key,key_l);
	gzwrite(FP,&val,sizeof(int));
	fprintf(stderr,"key_l: %d key:%s val:%d\n",key_l,key,val);
      }
      if(cm->find(key)!=cm->end())
	fprintf(stderr,"\t-> Duplicate entries found \'%s\'\n",key);
      (*cm)[strdup(key)]=val;
    }
  }

  if(FP)
    gzclose(FP);

  fprintf(stderr,"\n");
  fprintf(stderr,"\t-> [%s] Number of entries to use from accesion to taxid: %lu\n",fname,cm->size());
  return cm;
}



void *read_ass2taxid_thread(void *ptr){
  pars *p = (pars *) ptr;
  p->charref2taxid = ass2bin(p->acc2taxfile,0);
  pthread_mutex_unlock(&mutex2);
  return NULL;
}



pars *get_pars(int argc,char **argv){
  pars *p = pars_init();
  if(argc % 2){
    fprintf(stderr,"\t-> Must supply arguments in the form -pattern value\n");
    free(p);
    return NULL;
  }

  int make_used_reads = 0;
  
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    //    fprintf(stderr,"key: %s val: %s\n",key,val);
    if(!strcasecmp("-bam",key)){free(p->htsfile); p->htsfile=strdup(val);}
    else if(!strcasecmp("-names",key)){free(p->namesfile); p->namesfile=strdup(val);}
    else if(!strcasecmp("-nodes",key)){free(p->nodesfile); p->nodesfile=strdup(val);}
    else if(!strcasecmp("-acc2tax",key)){free(p->acc2taxfile); p->acc2taxfile=strdup(val);}
    else if(!strcasecmp("-editdistMin",key)) p->editdistMin=atoi(val);
    else if(!strcasecmp("-editdistMax",key)) p->editdistMax=atoi(val);
    else if(!strcasecmp("-minmapq",key)) p->minmapq=atoi(val);
    else if(!strcasecmp("-minlength",key)) p->minlength=atoi(val);
    else if(!strcasecmp("-simscoreLow",key)) p->simscoreLow=atof(val);
    else if(!strcasecmp("-lca_rank",key)){free(p->lca_rank); p->lca_rank=strdup(val);}
    else if(!strcasecmp("-simscoreHigh",key)) p->simscoreHigh=atof(val);
    else if(!strcasecmp("-outnames",key)){free(p->outnames); p->outnames=strdup(val);}
    else if(!strcasecmp("-out",key)){free(p->outnames); p->outnames=strdup(val);}
    else if(!strcasecmp("-fix-ncbi",key)) p->fixdb=atoi(val);
    else if(!strcasecmp("-fix_ncbi",key)) p->fixdb=atoi(val);
    else if(!strcasecmp("-discard",key)) p->discard=atoi(val);
    else if(!strcasecmp("-howmany",key)) p->howmany=atoi(val);
    else if(!strcasecmp("-usedreads",key)) make_used_reads = atoi(val);
    else if(!strcasecmp("-norank2species",key)) p->norank2species=atoi(val);
    else if(!strcasecmp("-skipnorank",key)) p->skipnorank=atoi(val);
    else if(!strcasecmp("-nthreads",key)) p->nthreads=atoi(val);
    else if(!strcasecmp("-weighttype",key)) p->weighttype=atoi(val);
    else if(!strcasecmp("-@",key)) p->nthreads=atoi(val);
    else if(!strcasecmp("-tempfolder",key)){free(p->tempfolder); p->tempfolder=strdup(val);}
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      free(p);
      return NULL;
    }
    
    ++argv;
  }

  pthread_t thread1, thread2;
  pthread_mutex_lock(&mutex1);
  //  pthread_mutex_lock(&mutex2);
  assert(pthread_create( &thread1, NULL, read_header_thread, (void*) p)==0);
  //  assert(pthread_create( &thread2, NULL, read_ass2taxid_thread, (void*) p)==0);
  
  char buf[1024];
  snprintf(buf,1024,"%s.lca",p->outnames);
  fprintf(stderr,"\t-> Will output lca results in file:\t\t\'%s\'\n",buf);
  p->fp1 = fopen(buf,"wb");
  assert(p->fp1);
  snprintf(buf,1024,"%s.lca.stat",p->outnames);
  fprintf(stderr,"\t-> Will output lca distribution in file:\t\t\'%s\'\n",buf);
  p->fp_lcadist = NULL;
  p->fp_lcadist = fopen(buf,"wb");
  assert(p->fp_lcadist);
  snprintf(buf,1024,"%s.wlca",p->outnames);
  fprintf(stderr,"\t-> Will output lca weight in file:\t\t\'%s\'\n",buf);
  //  p->fp2 = fopen(buf,"wb");

  snprintf(buf,1024,"%s.log",p->outnames);
  fprintf(stderr,"\t-> Will output log info (problems) in file:\t\'%s\'\n",buf);
  p->fp3 = fopen(buf,"wb");
  assert(p->fp3);
  if(make_used_reads){
    snprintf(buf,1024,"%s.usedreads.sam",p->outnames);
    fprintf(stderr,"\t-> Will output the reads that are used for damage file:\t\'%s\'\n",buf);
    p->usedreads_sam = strdup(buf);
  }
  
  

  pthread_mutex_lock(&mutex1);
  pthread_mutex_lock(&mutex2);
  return p;
}

void print_pars(FILE *fp,pars *p){
  fprintf(fp,"\t-> -bam  \t%s\n",p->htsfile);
  fprintf(fp,"\t-> -names\t%s\n",p->namesfile);
  fprintf(fp,"\t-> -nodes\t%s\n",p->nodesfile);
  fprintf(fp,"\t-> -acc2tax\t%s\n",p->acc2taxfile);
  fprintf(fp,"\t-> -simscoreLow\t%f\n",p->simscoreLow);
  fprintf(fp,"\t-> -simscoreHigh\t%f\n",p->simscoreHigh);
  fprintf(fp,"\t-> -editdistMin\t%d\n",p->editdistMin);
  fprintf(fp,"\t-> -editdistMax\t%d\n",p->editdistMax);
  fprintf(fp,"\t-> -outnames\t%s\n",p->outnames);
  fprintf(fp,"\t-> -minmapq\t%d\n",p->minmapq);
  fprintf(fp,"\t-> -lca_rank\t%s\n",p->lca_rank);
  fprintf(fp,"\t-> -norank2species\t%d\n",p->norank2species);
  fprintf(fp,"\t-> -howmany\t%d\n",p->howmany);
  fprintf(fp,"\t-> -fix_ncbi\t%d\n",p->fixdb);
  fprintf(fp,"\t-> -weighttype\t%d\n",p->weighttype);
  fprintf(fp,"\t-> -tempfolder\t%d\n",p->tempfolder);
}

#ifdef __WITH_MAIN__
int main(int argc,char**argv){
  pars *p=get_pars(--argc,++argv);
  assert(p);
  
  print_pars(stdout,p);
}
#endif
