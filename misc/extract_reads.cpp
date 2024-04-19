#include <cstdio>
#include <zlib.h>
#include <map>
#include <cstring>
#include <cstdlib>
#include <htslib/sam.h>
#include <cassert>
#include <ctype.h>
#include "../shared.h"

extern int SIG_COND;
char2int getkeys(const char *key,int value){
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
    cmap[strdup(tok)] = value;
  }
  fprintf(stderr,"\t-> Done reading keys from: \'%s\' nitems: %lu\n",key,cmap.size());
  fclose(fp);

  return cmap;
}

void runextract_taxid(int2int &keeplist,samFile *htsfp,bam_hdr_t *hdr,int strict,const char *outname){
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


void runextract_readid(char2int &keeplist,samFile *htsfp,bam_hdr_t *hdr,const char *outname,char out_mode[5],int complement){
  htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

  //open outputfile and write header
  samFile *outhts = NULL;
  if ((outhts = sam_open_format(outname,out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing: %s\n",outname);
    exit(0);
  }

  if (sam_hdr_write(outhts, hdr) < 0)
      fprintf(stderr,"Problem writing headers to %s", outname);

  //now mainloop
  bam1_t *aln = bam_init1();
  
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    char *qname = bam_get_qname(aln);
    char2int::iterator it = keeplist.find(qname);
    int isthere;
    if(it==keeplist.end())
      isthere=0;
    else{
      isthere=1;
      it->second = it->second + 1;
    }
    if(complement==0&&isthere==1)
      assert(sam_write1(outhts, hdr,aln)>=0);
    else if(complement==1&&isthere==0)
      assert(sam_write1(outhts, hdr,aln)>=0);
    
  }
  for(char2int::iterator it=keeplist.begin();it!=keeplist.end();it++){
    if(it->second==0)
      fprintf(stderr,"read: %s doesnt exists in screening file\n",it->first);

  }
  
  sam_hdr_destroy(hdr);
  sam_close(htsfp);
  sam_close(outhts);
}


void runextract_refid(int2int &keeplist,samFile *htsfp,bam_hdr_t *hdr,const char *outname,char out_mode[5],int complement){
  htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

  //open outputfile and write header
  samFile *outhts = NULL;
  if ((outhts = sam_open_format(outname,out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing: %s\n",outname);
    exit(0);
  }

  if (sam_hdr_write(outhts, hdr) < 0)
      fprintf(stderr,"Problem writing headers to %s", outname);

  //now mainloop
  bam1_t *aln = bam_init1();
  
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    int2int::iterator it = keeplist.find(aln->core.tid);
    int isthere;
    if(it==keeplist.end())
      isthere=0;
    else{
      isthere=1;
      it->second = it->second + 1;
    }
    if(complement==0&&isthere==1)
      assert(sam_write1(outhts, hdr,aln)>=0);
    else if(complement==1&&isthere==0)
      assert(sam_write1(outhts, hdr,aln)>=0);
    
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

int main_byrefid(int argc,char**argv){
  
  if(argc==1){
    fprintf(stderr,"./extract_reads byrefid -hts -key -outnames -type\n");
    return 0;
  }
  argv++;
  char *keyfile = NULL;
  char *hts = NULL;
  char *type = NULL;
  char *outfile = strdup("tmp.sam");
  char out_mode[5] = "wb";
  int docomplement = 0;
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-key",key)) keyfile=strdup(val);
    else if(!strcasecmp("-type",key)) {
      type = strdup(val);
      out_mode[1]=tolower(val[0]);
    }
    else if(!strcasecmp("-docomplement",key)) docomplement=atoi(val);
    else if(!strcasecmp("-out",key)) outfile=strdup(val);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      return 0;
    }
    ++argv;
  }
  
  fprintf(stderr,"\t-> key: %s \n\t-> hts: %s \n\t-> type: %s \n\t-> outfile: %s\n",keyfile,hts,type,outfile);

  //open inputfile and parse header
  samFile *htsfp = hts_open(hts,"r");
  bam_hdr_t *hdr = sam_hdr_read(htsfp); 

  char2int cmap = getkeys(keyfile,0);
  
  int2int keeplist;
  int VERB = 4;
  size_t counter[2] = {0,0};//there, not there
  for(char2int::iterator it=cmap.begin();it!=cmap.end();it++){
    int tokeep = sam_hdr_name2tid(hdr,it->first);
    if(tokeep>=0){
      keeplist[tokeep] =1;
      counter[0] = counter[0] +1;
      //      fprintf(stderr,"%s\n",it->first);
    }else{
      if(VERB>0){
	fprintf(stderr,"\t-> This id: %s does not exist in sam/bam/cramfile: %s\n",it->first,hts);
	fprintf(stderr,"\t-> This info is only printed %d more times\n",VERB);
	VERB--;
      }
      counter[1] = counter[1] +1;
    }

  }
  fprintf(stderr,"\t-> Number of refids to use: %lu from -key \'%s\'\n\t-> Number of refids notused: %lu\n",keeplist.size(),keyfile,counter[1]);
  runextract_refid(keeplist,htsfp,hdr,outfile,out_mode,docomplement);
  
  return 0;
}

int main_bytaxid(int argc,char**argv){
  
  if(argc==1){
    fprintf(stderr,"./extract_reads -hts -key -taxid -nodes -acc2txt -outnames -strict -type\n");
    return 0;
  }
  fprintf(stderr,"\t-> 19april2024, i dont remember this functionality. Use at own risk\n");
  argv++;
  char *keyfile = NULL;
  char *hts = NULL;
  char *names = NULL;
  char *nodefile = NULL;
  char *acc2tax = NULL;
  int strict = 0;
  int type = 0;
  char *outfile = strdup("tmp.sam");
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    //X    fprintf(stderr,"key: %s val: %s\n",key,val);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-key",key)) keyfile=strdup(val);
    else if(!strcasecmp("-taxid",key)) names=strdup(val);
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

  char2int cmap = getkeys(keyfile,0);
  
  char2int taxids;
  if(!fexists(names))
    taxids[names] = 1;
  else
    taxids = getkeys(names,0);

  for(char2int::iterator it=taxids.begin();it!=taxids.end();it++)
    fprintf(stderr,"\t-> looping up taxid: %s\n",it->first);
  
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
  if(hdr||1)
    bam2tax=(int2int*) bamRefId2tax(hdr,acc2tax,hts,errmap,NULL,0);
  fprintf(stderr,"\t-> We have bam2tax.size(): %lu and errmap.size():%lu\n",bam2tax->size(),errmap.size());
  

  //make a list of taxids to use
  int2int taxlist;
  for(char2int::iterator it=taxids.begin();it!=taxids.end();it++){
    gettaxids_to_use(atoi(it->first),child,taxlist);
    fprintf(stderr,"\t-> Number of taxids spanned by taxid: %d from nodesfile: %s is: %lu\n",atoi(it->first),nodefile,taxlist.size());
  }
  fprintf(stderr,"\t-> Total number of taxids to filter from: %lu\n",taxlist.size());
  int2int keeplist;
  for(char2int::iterator it=cmap.begin();it!=cmap.end();it++){
    int tokeep = sam_hdr_name2tid(hdr,it->first);
    assert(tokeep>=0);
    keeplist[tokeep] =1;

  }
  fprintf(stderr,"\t-> Number of refids to use: %lu from -key %s\n",keeplist.size(),keyfile);
  
  for(int2int::iterator it=bam2tax->begin();it!=bam2tax->end();it++){
    int2int::iterator it2=taxlist.find(it->second);
    if(it2!=taxlist.end())
      keeplist[it->first] = 1;
  }
  fprintf(stderr,"\t-> Number of refids to use: %lu\n",keeplist.size());
  
  runextract_taxid(keeplist,htsfp,hdr,strict,outfile);
  
  return 0;
}


int main_byreadid(int argc,char**argv){
  
  if(argc==1){
    fprintf(stderr,"./extract_reads -hts -key -out -type\n");
    return 0;
  }
  argv++;
  char *keyfile = NULL;
  char *hts = NULL;
  int type = 0;
  char *outfile = strdup("tmp.sam");
  char out_mode[5] = "wb";
  int docomplement = 0;
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    //X    fprintf(stderr,"key: %s val: %s\n",key,val);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-key",key)) keyfile=strdup(val);
    else if(!strcasecmp("-out",key)) outfile=strdup(val);
    else if(!strcasecmp("-docomplement",key)) docomplement=atoi(val);
    else if(!strcasecmp("-type",key)) out_mode[1]=tolower(val[0]);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      return 0;
    }
    ++argv;
  }
  
  fprintf(stderr,"\t-> key: %s \n\t-> hts: %s \n\t-> type: %d \n\t-> outfile: %s\n",keyfile,hts,type,outfile);

  //open inputfile and parse header
  samFile *htsfp = hts_open(hts,"r");
  bam_hdr_t *hdr = sam_hdr_read(htsfp); 

  char2int cmap = getkeys(keyfile,0);
  
  fprintf(stderr,"\t-> number of refids to use: %lu\n",cmap.size());
  
  runextract_readid(cmap,htsfp,hdr,outfile,out_mode,docomplement);
  
  
  return 0;
}


int main(int argc,char**argv){
  
  if(argc==1){
    fprintf(stderr,"./extract_reads bytaxid -hts -key -taxid -nodes -acc2txt -out -strict -type [dev]\n");
    fprintf(stderr,"./extract_reads byreadid -hts -key -out -type -docomplement\n");
    fprintf(stderr,"./extract_reads byrefid -hts -key -out -type -docomplement\n");
    return 0;
  }
  if(strcasecmp(argv[1],"bytaxid")==0)
    return main_bytaxid(argc--,++argv);
  if(strcasecmp(argv[1],"byreadid")==0)
    return main_byreadid(argc--,++argv);
   if(strcasecmp(argv[1],"byrefid")==0)
    return main_byrefid(argc--,++argv);

  fprintf(stderr,"\t-> Unknown options please use bytaxid or byreadid or byrefid\n");
  
  return 0;
}
