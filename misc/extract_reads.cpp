#include <cstdio>
#include <zlib.h>
#include <map>
#include <cstring>
#include <cstdlib>
#include <htslib/sam.h>
#include <cassert>
#include "../shared.h"
int SIG_COND=1;
char2int getkeys(const char *key){
  char2int cmap;
  if(key==NULL)
    return cmap;
  FILE *fp = NULL;
  
  if(((fp=fopen(key,"rb")))==NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",key);
    exit(0);
  }
  char buf[4096];
  while(fgets(buf,4096,fp)){
    char *tok = strtok(buf,"\n\t ");
    if(cmap.find(tok)!=cmap.end()){
      fprintf(stderr,"\t-> key: %s already exist will skip\n",buf);
    }
    cmap[strdup(tok)] = 1;
  }
  fprintf(stderr,"\t-> Done reading keys from: \'%s\' nitems: %lu\n",key,cmap.size());
  fclose(fp);

  return cmap;
}

//this is the old that should not be used
void runextract_sam(char2int &cmap,const char *fname,int strict){
  gzFile gz = Z_NULL;
  if(((gz=gzopen(fname,"rb")))==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    return;
  }
  char buf[4096];
  char orig[4096];
  while(gzgets(gz,buf,4096)){
    if(buf[0]=='@'){//header
      fprintf(stdout,"%s",buf);
      continue;
    }
    //fprintf(stderr,"Done reading header\n");
    //noheader
    char *orig = strncpy(orig,buf,4096);
    char *tok = strtok(buf,"\t\n ");
    tok = strtok(NULL,"\t\n ");
    tok = strtok(NULL,"\t\n ");//this is now chr/contig in sam
    if(cmap.find(tok)!=cmap.end()){
      fprintf(stdout,"%s",orig);
      continue;
    }
  }
  gzclose(gz);

}

void runextract(int2int &keeplist,samFile *htsfp,bam_hdr_t *hdr,int strict,const char *outname){
  htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

  //open outputfile and write header
  samFile *outhts = NULL;
  char out_mode[5]="ws";
  if ((outhts = sam_open_format(outname,out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing: %s\n",outname);
    exit(0);
  }

  queue *myq = init_queue(5000);//very large, should be enough.
  
  if (sam_hdr_write(outhts, hdr) < 0)
      fprintf(stderr,"Problem writing headers to %s", outname);

  //now mainloop
  bam1_t *aln = bam_init1();
  char *last=NULL;
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    char *qname = bam_get_qname(aln);
    int chr = aln->core.tid ; //contig name (chromosome)
    if(last==NULL)
      last=strdup(qname);
    //change of qname
    if(strcmp(last,qname)!=0) {
      if(strict==1){
	for(int i=0;i<myq->l;i++){
	  int2int::iterator it=keeplist.find(myq->ary[i]->core.tid);
	  if(it!=keeplist.end())
	    assert(sam_write1(outhts, hdr,myq->ary[i])>=0);
	}
      }else{
	int writedata = 0;
	for(int i=0;i<myq->l;i++){
	  int2int::iterator it=keeplist.find(myq->ary[i]->core.tid);
	  if(it!=keeplist.end()){
	    writedata=1;
	    break;
	  }
	}
	if(writedata>0){
	  for(int i=0;i<myq->l;i++){
	    assert(sam_write1(outhts, hdr,myq->ary[i])>=0);
	  }
	}
      }
      myq->l =0;
      last=strdup(qname);
    }
    assert(bam_copy1(myq->ary[myq->l],aln)!=NULL);
    myq->l++;
    if(myq->l==myq->m)
      expand_queue(myq);
  }
  if(strict==1){
    for(int i=0;i<myq->l;i++){
      int2int::iterator it=keeplist.find(myq->ary[i]->core.tid);
      if(it!=keeplist.end())
	assert(sam_write1(outhts, hdr,myq->ary[i])>=0);
    }
  }else{
    int writedata = 0;
    for(int i=0;i<myq->l;i++){
      int2int::iterator it=keeplist.find(myq->ary[i]->core.tid);
      if(it!=keeplist.end()){
	writedata=1;
	break;
      }
    }
    if(writedata>0){
      for(int i=0;i<myq->l;i++){
	assert(sam_write1(outhts, hdr,myq->ary[i])>=0);
      }
    }
  }
  
  sam_hdr_destroy(hdr);
  sam_close(htsfp);
  sam_close(outhts);
}

//taxid is the root the one where will find all lower noded(further from the root)
//child are the childnode structure
//i2i contains all the taxids that is spanned by the taxid. This could be a vector but we use hash to check internal structure that we have no loops
void gettaxids_to_use(int taxid,int2intvec &child,int2int &i2i){
  int2int::iterator it = i2i.find(taxid);
  assert(it==i2i.end());
  i2i[taxid] =1;
  int2intvec::iterator it2 = child.find(taxid);
  if (it2!=child.end()){
    std::vector<int> &avec = it2->second;
    for(int i=0;i<avec.size();i++){
      //	fprintf(stderr,"%d/%d %d\n",i,avec.size(),avec[i]);
      gettaxids_to_use(avec[i],child,i2i);
      
    }
  }


}

int main(int argc,char**argv){
  
  if(argc==3){
    fprintf(stderr,"./extract_reads -hts -key -names -nodes -acc2txt -outnames -strict -type\n");
    return 0;
  }
  argv++;
  char *keyfile = NULL;
  char *hts = NULL;
  char *names = NULL;
  char *nodefile = NULL;
  char *acc2tax = NULL;
  int strict = 0;
  int type = 0;
  char *outfile = "tmp.sam";
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    //X    fprintf(stderr,"key: %s val: %s\n",key,val);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-key",key)) keyfile=strdup(val);
    else if(!strcasecmp("-names",key)) names=strdup(val);
    else if(!strcasecmp("-nodes",key)) nodefile=strdup(val);
    else if(!strcasecmp("-acc2tax",key)) acc2tax=strdup(val);
    else if(!strcasecmp("-strict",key)) strict=atoi(val);
    else if(!strcasecmp("-out",key)) outfile=strdup(val);
    else if(!strcasecmp("-type",key)) type=atoi(val);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      return 0;
    }
    ++argv;
  }
  
  fprintf(stderr,"\t-> key: %s \n\t-> hts: %s \n\t-> nodefile: %s \n\t-> acc2tax: %s \n\t-> names: %s \n\t-> strict: %d \n\t-> type: %d \n\t-> outfile: %s\n",keyfile,hts,nodefile, acc2tax, names,strict,type,outfile);

  //open inputfile and parse header
  samFile *htsfp = hts_open(hts,"r");
  bam_hdr_t *hdr = sam_hdr_read(htsfp); 


  
  char2int cmap = getkeys(keyfile);
  
  
  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  //parsenodefile
  if(nodefile!=NULL)
    parse_nodes(nodefile,rank,parent,child,1);

  //make bamrefids to taxids
  int2int *bam2tax;//pointer uhhhh
  int2int errmap;
  if(hdr)
    bam2tax=(int2int*) bamRefId2tax(hdr,acc2tax,hts,errmap);
  fprintf(stderr,"\t-> We have bam2tax.size(): %lu and errmap.size():%lu\n",bam2tax->size(),errmap.size());
  

  //make a list of taxids to use
  int2int taxlist;
  if(names!=NULL){
    gettaxids_to_use(atoi(names),child,taxlist);
    fprintf(stderr,"\t-> Number of taxids spanned by taxid: %s from nodesfile: %s is: %lu\n",names,nodefile,taxlist.size());
  }
  int2int keeplist;
  for(char2int::iterator it=cmap.begin();it!=cmap.end();it++){
    int tokeep = sam_hdr_name2tid(hdr,it->first);
    assert(tokeep>=0);
    keeplist[tokeep] =1;

  }
  
  for(int2int::iterator it=bam2tax->begin();it!=bam2tax->end();it++){
    int2int::iterator it2=taxlist.find(it->second);
    if(it2!=taxlist.end())
      keeplist[it->first] = 1;
  }
  fprintf(stderr,"\t-> number of refids to use: %lu\n",keeplist.size());
  
  runextract(keeplist,htsfp,hdr,strict,outfile);
  
  return 0;
}
