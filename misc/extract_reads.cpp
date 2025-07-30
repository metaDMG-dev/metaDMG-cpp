#include <cstdio>
#include <zlib.h>
#include <map>
#include <cstring>
#include <cstdlib>
#include <htslib/sam.h>
#include <cassert>
#include <ctype.h>
#include "../shared.h"

//to make life easier we are making some global variables.
char out_mode[8] = "wb";//<- assume output is a bam
htsFormat fmt{};//<- stupid syntax {} ls vs {0} rs

void set_output_format_and_mode(const char *outname, const char *user_fmt_str) {
  
  // Determine format from file extension
  htsExactFormat fmt_from_ext = bam;  // Default = BAMblam
  if (strcmp(outname, "-") != 0) {
    const char *dot = strrchr(outname, '.');
    if (dot && dot[1]) {
      if (strcasecmp(dot + 1, "sam") == 0)
	fmt_from_ext = sam;
      else if (strcasecmp(dot + 1, "bam") == 0)
	fmt_from_ext = bam;
      else if (strcasecmp(dot + 1, "cram") == 0)
	fmt_from_ext = cram;
    }
  }

  // If user specifies format
  if (user_fmt_str) {
    if (strcasecmp(user_fmt_str, "sam") == 0) {
      fmt.format = sam;
      strcpy(out_mode, "w");
    } else if (strcasecmp(user_fmt_str, "bam") == 0) {
      fmt.format = bam;
      strcpy(out_mode, "wb");
    } else if (strcasecmp(user_fmt_str, "cram") == 0) {
      fmt.format = cram;
      strcpy(out_mode, "wc");
    } else {
      fprintf(stderr, "Unknown output format: %s\n", user_fmt_str);
      exit(1);
    }
    
    // Validate against extension unless writing to stdout
    if (strcmp(outname, "-") != 0 && fmt.format != fmt_from_ext) {
      fprintf(stderr,
	      "Error: user-specified format '%s' conflicts with output file extension '%s'.\n",
	      user_fmt_str, outname);
      exit(1);
    }
    return;
  }

  // No user format — use inferred or default (BAM)
  fmt.format = fmt_from_ext;
  if (fmt.format == sam)
    strcpy(out_mode, "w");
  else if (fmt.format == cram)
    strcpy(out_mode, "wc");
  else
    strcpy(out_mode, "wb");  // BAM default
}


extern int SIG_COND;
char2int getkeys(const char *key,int value,int nospace){
  char *delim=NULL;
  if(nospace==0)
    delim = strdup("\n\t ");
  else if(nospace==1)
    delim = strdup("\n");
  char2int cmap;
  if(key==NULL)
    return cmap;
  FILE *fp = NULL;
  
  if(((fp=fopen(key,"rb")))==NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",key);
    exit(-1);
  }
  char buf[4096];
  while(fgets(buf,4096,fp)){
    char *tok = strtok(buf,delim);
    if(cmap.find(tok)!=cmap.end()){
      fprintf(stderr,"\t-> key: %s already exist will skip\n",buf);
    }
    cmap[strdup(tok)] = value;
  }
  fprintf(stderr,"\t-> Done reading keys from: \'%s\' nitems: %lu\n",key,cmap.size());
  fclose(fp);
  free(delim);
  return cmap;
}

int2int getkeysint(const char *key,int value){
  int2int cmap;
  if(key==NULL)
    return cmap;
  FILE *fp = NULL;
  
  if(((fp=fopen(key,"rb")))==NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",key);
    exit(-1);
  }
  char buf[4096];
  while(fgets(buf,4096,fp)){
    char *tok = strtok(buf,"\n\t ");
    if(cmap.find(atoi(tok))!=cmap.end()){
      fprintf(stderr,"\t-> key: %s already exist will skip\n",buf);
    }
    cmap[atoi(tok)] = value;
  }
  fprintf(stderr,"\t-> Done reading keys from: \'%s\' nitems: %lu\n",key,cmap.size());
  fclose(fp);

  return cmap;
}

void doflush(queue *myq,int2int &keeplist,bam_hdr_t *hdr,samFile *outhts,int strict){
  // fprintf(stderr,"flush: %lu strictk:%d\n",myq->l,strict);
  if(strict==1){//will only print specific match
    for(int i=0;i<(int)myq->l;i++){
      int2int::iterator it=keeplist.find(myq->ary[i]->core.tid);
      if(it!=keeplist.end())
	assert(sam_write1(outhts, hdr,myq->ary[i])>=0);
    }
  }
  if(strict==0){
    //if writedata=0 then no aligments will be printed, otherwise all
    int writedata = 0;
    for(int i=0;i<(int)myq->l;i++){
      int2int::iterator it=keeplist.find(myq->ary[i]->core.tid);
      if(it!=keeplist.end()){
	writedata=1;
	break;
      }
    }
    if(writedata>0){
      for(int i=0;i<(int)myq->l;i++){
	assert(sam_write1(outhts, hdr,myq->ary[i])>=0);
      }
    }
  }
  myq->l =0;
}

//strict=0 means only refids matching
//strict=1 means all aln will be included, if there is a match for one of the alignments
void runextract_int2int(int2int &keeplist,samFile *htsfp,bam_hdr_t *hdr,const char *outname,char *type,int strict){
  fprintf(stderr,"outname: %s outformat: %s\n",outname,out_mode);
  set_output_format_and_mode(outname,type);


  //open outputfile and write header
  samFile *outhts = NULL;
  if ((outhts = sam_open_format(outname,out_mode, &fmt)) == NULL) {
    fprintf(stderr,"Error opening file for writing: %s\n",outname);
    exit(-1);
  }

  queue *myq = init_queue(5000);//very large, should be enough.
  
  if (sam_hdr_write(outhts, hdr) < 0)
      fprintf(stderr,"Problem writing headers to %s", outname);

  //now mainloop
  bam1_t *aln = bam_init1();
  char *last=NULL;
  while(sam_read1(htsfp,hdr,aln) >= 0) {
    char *qname = bam_get_qname(aln);
    if(last==NULL)
      last=strdup(qname);
    
    //change of qname
    if(strcmp(last,qname)!=0) {
      doflush(myq,keeplist,hdr,outhts,strict);
      myq->l = 0;
      last=strdup(qname);
    }
    assert(bam_copy1(myq->ary[myq->l],aln)!=NULL);
    myq->l++;
    if(myq->l==myq->m)
      expand_queue(myq);
  }
  doflush(myq,keeplist,hdr,outhts,strict);
  //  doflush(myq,strict,outhts);
  myq->l = 0;
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
    exit(-1);
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
      fprintf(stderr,"read: %s doesnt exists in inputfile file\n",it->first);

  }
  
  sam_hdr_destroy(hdr);
  sam_close(htsfp);
  sam_close(outhts);
}


//taxid is the root the one where will find all lower noded(further from the root)
//child are the childnode structure
//i2i contains all the taxids that is spanned by the taxid. This could be a vector but we use hash to check internal structure that we have no loops
void gettaxids_to_use(int taxid,int2intvec &child,int2int &i2i){
  i2i[taxid] =1;
  int2intvec::iterator it = child.find(taxid);
  if(it==child.end()){
    fprintf(stderr,"\t-> Problem finding taxid: %d from nodesfile\n",taxid);
  }else{ 
    std::vector<int> &avec = it->second;
    for(int i=0;i<(int)avec.size();i++){
      //	fprintf(stderr,"%d/%d %d\n",i,avec.size(),avec[i]);
      gettaxids_to_use(avec[i],child,i2i);
      
    }
  }


}

int main_byrefid(int argc,char**argv){
  if(argc==2){
    fprintf(stderr,"./extract_reads byrefid -hts -key -out -type -strict \n");
    return 0;
  }
  argv++;
  char *keyfile = NULL;
  char *hts = NULL;
  char *type = NULL;
  char *outfile = strdup("tmp.bam");
  int docomplement = 0;
  char out_mode[5] = "wb";
  int strict = 1;
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-key",key)) keyfile=strdup(val);
    else if(!strcasecmp("-docomplement",key)) docomplement=atoi(val);
    else if(!strcasecmp("-type",key)) {
      type = strdup(val);
      out_mode[1]=tolower(val[0]);
    }
    else if(!strcasecmp("-strict",key)) strict=atoi(val);
    else if(!strcasecmp("-out",key)) outfile=strdup(val);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      return 0;
    }
    ++argv;
  }
  
  fprintf(stderr,"\t-> key: %s \n\t-> hts: %s \n\t-> type: %s \n\t-> outfile: %s\n\t-> strict: %d\n",keyfile?keyfile:"(null)",hts?hts:"(null)",type,outfile?outfile:"(null)",strict);

  //open inputfile and parse header
  samFile *htsfp = hts_open(hts,"r");
  bam_hdr_t *hdr = sam_hdr_read(htsfp); 

  char2int cmap = getkeys(keyfile,0,0);
  
  int2int keeplist;
  int VERB = 4;
  size_t counter[2] = {0,0};//there, not there
  if(docomplement==0){
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
  }else{
    for (int i = 0; i < hdr->n_targets; i++) {
      int tokeep = -1;
      //printf("Reference %d: name = %s, length = %d\n", i, header->target_name[i], header->target_len[i]);
      char2int::iterator it = cmap.find(hdr->target_name[i]);
      if(it==cmap.end())
	tokeep = i;

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
  }
  
  fprintf(stderr,"\t-> Number of refids to use: %lu from -key \'%s\'\n\t-> Number of refids notused: %lu\n",keeplist.size(),keyfile,counter[1]);
  runextract_int2int(keeplist,htsfp,hdr,outfile,out_mode,strict);
  
  return 0;
}

int main_bytaxid(int argc,char**argv){
  if(argc==2){
    fprintf(stderr,"./extract_reads bytaxid -hts [-key file_with_refnames] -taxid -nodes -acc2tax -accout -strict -forcedump -out -type [s/b]am -taxnames\n");
    fprintf(stderr,"-------\nAlso -forcedump 1 -accout filename.txt.gz\n-strict 1 means only alignments that match\n-strict 0 (default) means all alignments associated with a read if one of the alignments match\n---------\n");
    fprintf(stderr,"examples:\n");
    fprintf(stderr,"./extract_reads bytaxid -hts yo.bam -taxid 3258 -nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp -acc2tax /projects/caeg/data/db/mikkels/combined_accession2taxid_20221112.gz -type bam -out tmp3.bam -strict 0\n");
    fprintf(stderr,"\nExtract all those reads where one of the alignments is a child to the node given by taxid 3258\n");
    return 0;
  }

  argv++;
  char *keyfile = NULL;
  char *hts = NULL;
  char *taxid = NULL;
  char *nodefile = NULL;
  //  char out_mode[5] = "wb";
  char *acc2tax = NULL;
  int strict = 0;
  char *type = NULL;
  char *outfile = strdup("tmp.bam");
  int forcedump = 0;
  char *accout = NULL;
  char *taxnames = NULL;
  char *names = NULL;
  while(*argv){
    char *key=*argv;
    char *val=*(++argv);
    //    fprintf(stderr,"key: %s val: %s\n",key,val);
    if(!strcasecmp("-hts",key)) hts=strdup(val);
    else if(!strcasecmp("-key",key)) keyfile=strdup(val);
    else if(!strcasecmp("-taxid",key)) taxid=strdup(val);
    else if(!strcasecmp("-taxnames",key)) taxnames=strdup(val);
    else if(!strcasecmp("-names",key)) names=strdup(val);
    else if(!strcasecmp("-nodes",key)) nodefile=strdup(val);
    else if(!strcasecmp("-acc2tax",key)) acc2tax=strdup(val);
    else if(!strcasecmp("-accout",key)) accout=strdup(val);
    else if(!strcasecmp("-strict",key)) strict=atoi(val);
    else if(!strcasecmp("-forcedump",key)) forcedump=atoi(val);
    else if(!strcasecmp("-out",key)) outfile=strdup(val);
    else if(!strcasecmp("-type",key))
      type=strdup(val);
    else{
      fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
      return 0;
    }
    ++argv;
  }
  
  fprintf(stderr,"\t-> key: %s \n\t-> hts: %s \n\t-> nodefile: %s \n\t-> acc2tax: %s \n\t-> taxid: %s\n\t-> taxnames: %s \n\t-> strict: %d \n\t-> type: %s \n\t-> outfile: %s\n\t-> forcedump:%d\n\t-> accout:%s\n\t-> names: %s\n",keyfile?keyfile:"(null)",hts?hts:"(null)",nodefile, acc2tax, taxid,taxnames,strict,type,outfile,forcedump,accout,names);

  if(taxid ==NULL&&taxnames==NULL){
    fprintf(stderr,"\t-> Need to supply -taxid and/or -taxnames\n");
    return 0;
  }
  if(acc2tax==NULL){
    fprintf(stderr,"\t-> Must supply -acc2tax\n");
    return 0;
  }
  int2int taxids;
  if(taxid!=NULL) {
    if(!fexists(taxid)) {
      char *tok = strtok(taxid, ",");
      while (tok != NULL) {
	taxids[atoi(tok)] = 1;
	tok = strtok(NULL, ",");
      }
    } else {
      taxids = getkeysint(taxid,0);
    }
  }
  
  //open inputfile and parse header
  samFile *htsfp = hts_open(hts,"r");
  bam_hdr_t *hdr = sam_hdr_read(htsfp); 
 
   //make bamrefids to taxids
  int2int *bam2tax;//pointer uhhhh
  int2int errmap;
  
  if(hdr||1)
    bam2tax=(int2int*) bamRefId2tax(hdr,acc2tax,hts,strdup("/tmp/"),forcedump,accout,NULL);
  fprintf(stderr,"\t-> We have bam2tax.size(): %lu and errmap.size():%lu \n",bam2tax->size(),errmap.size());

  if(taxnames!=NULL){
    if(names==NULL){
      fprintf(stderr,"\t-> if taxnames is defined then names should also be defined\n");
      exit(-1);
    }
    int2char nammap = parse_names(names);
    char2int mapnam;
    for(auto it=nammap.begin();it!=nammap.end();it++){
      auto it2=mapnam.find(it->second);
      if(it2!=mapnam.end()){
	fprintf(stderr,"\t-> Problem with duplicate name in name map: key %d val: %s\n",it->first,it->second);
	exit(-1);
      }
      mapnam[it->second] = it->first;
    }
    if(!fexists(taxnames)) {
      char *tok = strtok(taxnames, ",");
      while (tok != NULL) {
	auto it = mapnam.find(tok);
	if(it==mapnam.end()){
	  fprintf(stderr,"\t-> Problem finding taxname: %s\n",tok);
	  exit(-1);
	}
	taxids[it->second] = 1;
	tok = strtok(NULL, ",");
      }
    } else {
      char2int tmp = getkeys(taxnames,0,1);//"name"
      for(auto it=tmp.begin();it!=tmp.end();it++){
	auto it2 = mapnam.find(it->first);
	if(it2==mapnam.end()){
	  fprintf(stderr,"\t-> Problem finding taxname: %s\n",it->first);
	  exit(-1);
	}
	taxids[it2->second] = 1;
      }
    }
  }

  
  if(taxids.size()==0)
    return 0;

  for(int2int::iterator it=taxids.begin();0&&it!=taxids.end();it++)
    fprintf(stderr,"\t-> looking up taxid: %d\n",it->first);

  char2int cmap = getkeys(keyfile,0,0);
  

  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  //parsenodefile
  if(nodefile!=NULL)
    parse_nodes(nodefile,rank,parent,child,1);
 

  

  //make a list of taxids to use
  int2int taxlist;
  for(int2int::iterator it=taxids.begin();it!=taxids.end();it++){
    gettaxids_to_use(it->first,child,taxlist);
    fprintf(stderr,"\t-> Number of taxids spanned by taxid: %d from nodesfile: %s is: %lu\n",it->first,nodefile,taxlist.size());
  }
  fprintf(stderr,"\t-> Total number of taxids to filter from: %lu\n",taxlist.size());
#if 0//print out taxids to be included
  for(int2int::iterator it=taxlist.begin();it!=taxlist.end();it++)
    fprintf(stderr,"taxs: %d %d\n",it->first,it->second);
#endif
  
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
  if(keeplist.size()==0)
    fprintf(stderr,"\t-> No ids to extract\n");
  else
    runextract_int2int(keeplist,htsfp,hdr,outfile,type,strict);
  
  return 0;
}


int main_byreadid(int argc,char**argv){
  if(argc==2){
    fprintf(stderr,"./extract_reads byreadid -hts -key -out -type\n");
    fprintf(stderr,"./extract_reads byreadid -key keyfile -hts compr.bam\n");
    return 0;
  }
  
  argv++;
  char *keyfile = NULL;
  char *hts = NULL;
  int type = 0;
  char *outfile = strdup("tmp.bam");
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

  char2int cmap = getkeys(keyfile,0,0);
  
  fprintf(stderr,"\t-> number of refids to use: %lu\n",cmap.size());
  
  runextract_readid(cmap,htsfp,hdr,outfile,out_mode,docomplement);
  return 0;
}


int main(int argc,char**argv){
  
  if(argc==1){
    fprintf(stderr,"./extract_reads bytaxid -hts [-key file_with_refnames] -taxid -nodes -acc2txt -out -strict -type [dev]\n");
    fprintf(stderr,"./extract_reads byreadid -hts -key -out -type -docomplement\n");
    fprintf(stderr,"./extract_reads byrefid -hts -key -out -type -docomplement\n");
    fprintf(stderr,"-type is outputtype -type sam -type bam\n");
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
