//gpl thorfinn@binf.ku.dk
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctype.h>
#include <getopt.h>
#include <cassert>
#include <time.h>
#include <zlib.h>
#include <map>
#include "profile.h"
#include "shared.h"
#include "version.h"

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));
typedef std::map<int,char *> int2char;
int usage_getdamage(FILE *fp){
  fprintf(fp,"\nUsage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>\n");
  fprintf(fp,"\nExample: ./metadamage getdamage -l 10 -p 5 --threads 8 ../data/subs.sam\nOptions:\n");
  fprintf(fp,"  -f/--fasta\t is required with CRAM\n");
  fprintf(fp,"  -l/--minlength\t reads shorter than minlength will be discarded\n");
  fprintf(fp,"  -p/--printlength\t use this number of positions from 5' and 3'\n");
  fprintf(fp,"  -r/--runmode\t runmode 1 means that damage patterns will be calculated for each chr/scaffold contig.\n\t\t runmode 0 means one global estimate.\n");
  fprintf(fp,"  -@/--threads\t Number of threads used for reading/writing\n");
  fprintf(fp,"  -o/--outname\t Number of threads used for reading/writing\n");
  return 0;
}

double *getval(std::map<int, double *> &retmap,int2intvec &child,int taxid,int howmany){
  // fprintf(stderr,"getval\t%d\t%d\n",taxid,howmany);
  std::map<int,double *>::iterator it = retmap.find(taxid);
  if(it!=retmap.end()){
    //fprintf(stderr,"has found: %d in retmap\n",it->first);
#if 0
    fprintf(stderr,"val\t%d",taxid);
    for(int i=0;i<3*howmany+1;i++)
      fprintf(stderr,"\t%f",it->second[i]);
    fprintf(stderr,"\n");
#endif
    return it->second;
  }
  double *ret = new double [3*howmany+1];
  for(int i=0;i<3*howmany+1;i++)
    ret[i] = 0.0;
  if(child.size()>0) {// if we have supplied -nodes
    int2intvec::iterator it2 = child.find(taxid);
    if (it2!=child.end()){
      std::vector<int> &avec = it2->second;
      for(int i=0;i<avec.size();i++){
	//	fprintf(stderr,"%d/%d %d\n",i,avec.size(),avec[i]);
	double *tmp = getval(retmap,child,avec[i],howmany);
	for(int i=0;i<3*howmany+1;i++)
	  ret[i] += tmp[i];

      }
    }
  }
  
  retmap[taxid] = ret;

  return ret;
}

mydata getval_full(std::map<int, mydata> &retmap,int2intvec &child,int taxid,int howmany){
  // fprintf(stderr,"getval\t%d\t%d\n",taxid,howmany);
  std::map<int,mydata>::iterator it = retmap.find(taxid);
  if(it!=retmap.end()){
    //fprintf(stderr,"has found: %d in retmap\n",it->first);
#if 0
    fprintf(stderr,"val\t%d",taxid);
    for(int i=0;i<3*howmany+1;i++)
      fprintf(stderr,"\t%f",it->second[i]);
    fprintf(stderr,"\n");
#endif
    return it->second;
  }
  mydata ret;
  ret.nreads = 0;
  ret.fw = new size_t[16*howmany];
  ret.bw = new size_t[16*howmany];
  
  for(int i=0;i<16*howmany;i++){
    ret.fw[i] = 0;
    ret.bw[i] = 0;
  }
  if(child.size()>0) {// if we have supplied -nodes
    int2intvec::iterator it2 = child.find(taxid);
    if (it2!=child.end()){
      std::vector<int> &avec = it2->second;
      for(int i=0;i<avec.size();i++){
	//	fprintf(stderr,"%d/%d %d\n",i,avec.size(),avec[i]);
	mydata tmp = getval_full(retmap,child,avec[i],howmany);
	ret.nreads += tmp.nreads;
	for(int i=0;i<16*howmany;i++){
	  ret.fw[i] += tmp.fw[i];
	  ret.bw[i] += tmp.bw[i];
	}
      }
    }
  }
  
  retmap[taxid] = ret;

  return ret;
}

int main_getdamage(int argc,char **argv){
  if(argc==1)
    return usage_getdamage(stderr);

  //  int MAXLENGTH = 256;
  int minLength = 35;
  int printLength = 5;
  char *refName = NULL;
  char *fname = NULL;
  int runmode =0;//this means one species, runmode=1 means multi species
  htsFile *fp = NULL;
  char *onam = strdup("meta");
  int nthreads = 4;
  //fix thesepro
  static struct option lopts[] = {
    {"fasta", 1, 0, 'f'},
    {"minlength", 1, 0, 'l'},
    {"threads", 1, 0, '@'},
    {"printlength", 1, 0, 'p'},
    {"outname", 1, 0, 'o'},
    {"help", 0, 0, '?'},
    {"runmode", 1, 0, 'r'},
    {NULL, 0, NULL, 0}
  };
  
  int c;
  while ((c = getopt_long(argc, argv,
			  "f:l:p:r:o:@:",
			  lopts, NULL)) >= 0) {
    switch (c) {
    case 'f': refName = strdup(optarg); break;
    case 'l': minLength = atoi(optarg); break;
    case '@': nthreads = atoi(optarg); break;
    case 'p': printLength = atoi(optarg); break;
    case 'o': {free(onam);onam = strdup(optarg); break;}
    case 'r': runmode = atoi(optarg); break;
    case '?':
      return usage_getdamage(stdout);

    default:
      fprintf(stderr,"Never here: %s %s\n",optarg,fname);
      break;
    }
  }
  if(optind<argc)
    fname = strdup(argv[optind]);
  fprintf(stderr,"./metadamage refName: %s minLength: %d printLength: %d runmode: %d outname: %s nthreads: %d\n",refName,minLength,printLength,runmode,onam,nthreads);
  if(fname==NULL){
    usage_getdamage(stderr);
    return 0;
  }
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
  damage *dmg = new damage(printLength,nthreads,0);
  int skipper[4] = {3,3,3,3};
  while(((ret=sam_read1(fp,hdr,b)))>=0){
    if(bam_is_unmapped(b) ){
      if(skipper[0])
      	fprintf(stderr,"skipping: %s unmapped, this msg is printed: %d times more\n",bam_get_qname(b),--skipper[0]);
      continue;
    }
    if(bam_is_failed(b) ){
      if(skipper[1])
	fprintf(stderr,"skipping: %s failed: flags=%d, this msg is printed: %d times more\n",bam_get_qname(b),b->core.flag,--skipper[1]);
      continue;
    }
    if(b->core.l_qseq < minLength){
      if(skipper[2])
	fprintf(stderr,"skipping: %s too short, this msg is printed %d times more \n",bam_get_qname(b),--skipper[2]);
      continue;
    }

    dmg->damage_analysis(b,runmode!=0?b->core.tid:0);
  }


  dmg->printit(stdout,printLength);
  
  dmg->write(onam,runmode==1?hdr:NULL);
  dmg->bwrite(onam,hdr);
    
  sam_hdr_destroy(hdr);
  bam_destroy1(b);
  sam_close(fp);
  destroy_damage(dmg);
  free(fname);
  free(onam);
  return 0;
}

int main_index(int argc,char **argv){
  char *infile = argv[1];
  fprintf(stderr,"infile: %s\n",infile);
  char onam[strlen(infile)+20];
  sprintf(onam,"%s.idx",infile);
  fprintf(stderr,"outfile: %s\n",onam);
  FILE *fp = NULL;
  if(((fp=fopen(onam,"wb")))==NULL){
    fprintf(stderr,"Problem opening file\n");
    return 0;
  }

  fclose(fp);
  return 0;
}


int main_print(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"./metadamage print file.bdamage.gz [-names names.gz -bam file.bam -ctga -countout -nodes -howmany -r -doOld -nodes]\n");
    return 0;
  }
  char *infile = NULL;
  char *inbam = NULL;
  char *infile_nodes = NULL;
  char *infile_names = NULL;
  int ctga =0;//only print ctga errors
  int search = -1;
  int countout = 0;

  int howmany = 15;
  int doold = 1;
  while(*(++argv)){
    if(strcasecmp("-names",*argv)==0)
      infile_names = strdup(*(++argv));
    else if(strcasecmp("-bam",*argv)==0)
      inbam = strdup(*(++argv));
    else if(strcasecmp("-r",*argv)==0)
      search = atoi(*(++argv));
    else if(strcasecmp("-howmany",*argv)==0)
      howmany = atoi(*(++argv));
    else if(strcasecmp("-ctga",*argv)==0)
      ctga =1;
    else if(strcasecmp("-doOld",*argv)==0)
      doold = 1;
    else if(strcasecmp("-countout",*argv)==0)
      countout =1;
    else if(strcasecmp("-nodes",*argv)==0)
      infile_nodes = strdup(*(++argv));
    else
      infile = strdup(*argv);
  }


  fprintf(stderr,"infile: %s inbam: %s search: %d ctga: %d countout: %d nodes: %s names: %s howmany: %d\n",infile,inbam,search,ctga,countout,infile_nodes,infile_names,howmany);

  assert(infile);
  int2char name_map;
  if(infile_names!=NULL)
    name_map = parse_names(infile_names);

  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  if(infile_nodes!=NULL)
    parse_nodes(infile_nodes,rank,parent,child,1);
  if(search!=-1&& doold==0){
    std::map<int, double *> retmap = load_bdamage3(infile,howmany);
    double *dbl = getval(retmap,child,search,howmany);
    double dbldbl[3*howmany+1];//3 because ct,ga,other
    dbldbl[0] = dbl[0];
    for(int i=0;i<3*howmany;i++)
      dbldbl[i+1] = dbl[1+i]/dbl[0];
    
    fprintf(stdout,"%d\t%.0f",search,dbldbl[0]);
    for(int i=0;i<3*howmany;i++)
      fprintf(stdout,"\t%f",dbldbl[1+i]);
    fprintf(stdout,"\n");
    return 0;
  }
  
  BGZF *bgfp = NULL;
  samFile *samfp = NULL;
  bam_hdr_t *hdr =NULL;

  if(((bgfp = bgzf_open(infile, "r")))== NULL){
    fprintf(stderr,"Could not open input BAM file: %s\n",infile);
    return 1;
  }
  
  if(inbam!=NULL){
    if(((  samfp = sam_open_format(inbam, "r", NULL) ))== NULL){
      fprintf(stderr,"Could not open input BAM file: %s\n",inbam);
      return 1;
    }
    if(((hdr= sam_hdr_read(samfp))) == NULL){
      fprintf(stderr,"Could not read header for: %s\n",inbam);
    return 1;
    }
  }
  
  int printlength;
  assert(sizeof(int)==bgzf_read(bgfp,&printlength,sizeof(int)));
  fprintf(stderr,"\t-> printlength(howmany) from inputfile: %d\n",printlength);

  int ref_nreads[2];

  if(ctga==0){
    if(hdr!=NULL)
      fprintf(stdout,"#Reference\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
    else if(infile_names!=NULL)
      fprintf(stdout,"#FunkyName\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
    else
      fprintf(stdout,"#taxid\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
  }else{
    

  }
    
  int data[16];
  while(1){
    int nread=bgzf_read(bgfp,ref_nreads,2*sizeof(int));
    double ctgas[2*printlength];
    if(nread==0)
      break;
    assert(nread==2*sizeof(int));
    for(int at=0;at<printlength;at++){
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      if(at>howmany)
	continue;
      if(search==-1||search==ref_nreads[0]) {
	if(ctga==0)  {
	  if(hdr!=NULL)
	    fprintf(stdout,"%s\t%d\t5\'\t%d",hdr->target_name[ref_nreads[0]],ref_nreads[1],at);
	  else if(infile_names!=NULL){
	    int2char::iterator itt = name_map.find(ref_nreads[0]);
	    if(itt==name_map.end()){
	      fprintf(stderr,"\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n",ref_nreads[0],infile_names);
	      exit(0);
	    }
	    fprintf(stdout,"\"%s\"\t%d\t5\'\t%d",itt->second,ref_nreads[1],at);
	  }else
	    fprintf(stdout,"%d\t%d\t5\'\t%d",ref_nreads[0],ref_nreads[1],at);
	}else{
	  if(at==0)
	    fprintf(stdout,"%d\t%d",ref_nreads[0],ref_nreads[1]);
	}
	if(countout==1){
	  for(int i=0;i<16;i++)
	    fprintf(stdout,"\t%d",data[i]);
	  fprintf(stdout,"\n");
	}else{
	  float flt[16];
	  
	  for(int i=0;i<4;i++){
	    double tsum =0;
	    for(int j=0;j<4;j++){
	      tsum += data[i*4+j];
	      flt[i*4+j] = data[i*4+j];
	    }
	    if(tsum==0) tsum = 1;
	    for(int j=0;j<4;j++)
	      flt[i*4+j] /=tsum;
	  }
	  
	  if(ctga==0){
	    for(int j=0;j<16;j++)
	      fprintf(stdout,"\t%f",flt[j]);
	    fprintf(stdout,"\n");
	  }else
	    ctgas[at] = flt[7];
	}
      }
    }
    if(search==-1||search==ref_nreads[0]){
      if(ctga==1){
	for(int i=0;i<printlength;i++)
	  fprintf(stdout,"\t%f",ctgas[i]);
      }
    }
  
    for(int at=0;at<printlength;at++) {
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      if(at>howmany)
	continue;
      if(search==-1||search==ref_nreads[0]){
	if(ctga==0) {
	  if(hdr!=NULL)
	    fprintf(stdout,"%s\t%d\t3\'\t%d",hdr->target_name[ref_nreads[0]],ref_nreads[1],at);
	  else if(infile_names!=NULL){
	    int2char::iterator itt = name_map.find(ref_nreads[0]);
	    if(itt==name_map.end()){
	      fprintf(stderr,"\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n",ref_nreads[0],infile_names);
	      exit(0);
	    }
	    fprintf(stdout,"\"%s\"\t%d\t3\'\t%d",itt->second,ref_nreads[1],at);
	  }
	  else
	    fprintf(stdout,"%d\t%d\t3\'\t%d",ref_nreads[0],ref_nreads[1],at);
	}
	if(countout==1){
	  for(int i=0;i<16;i++)
	    fprintf(stdout,"\t%d",data[i]);
	  fprintf(stdout,"\n");
	}else{	
	  float flt[16];
	  
	  for(int i=0;i<4;i++){
	    double tsum =0;
	    for(int j=0;j<4;j++){
	      tsum += data[i*4+j];
	      flt[i*4+j] = data[i*4+j];
	    }
	    if(tsum==0) tsum = 1;
	    for(int j=0;j<4;j++)
	      flt[i*4+j] /=tsum;
	  }
	  if(ctga==0){
	    for(int j=0;j<16;j++)
	      fprintf(stdout,"\t%f",flt[j]);
	    fprintf(stdout,"\n");
	  }else
	  ctgas[at+printlength] = flt[8];
	}
      }
    }
    if(search==-1||search==ref_nreads[0]){
      if(ctga==1){
	for(int i=0;i<printlength;i++)
	  fprintf(stdout,"\t%f",ctgas[printlength+i]);
	fprintf(stdout,"\n");
      }
    }
  }

  if(bgfp)
    bgzf_close(bgfp);
  if(hdr)
    bam_hdr_destroy(hdr);
  if(samfp)
    sam_close(samfp);
  return 0;
}

int main_print2(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"./metadamage print2 file.bdamage.gz [-acc2tax file.gz -bam file.bam -ctga -countout -nodes -howmany -r -doOld -nodes]\n");
    return 0;
  }
  char *infile = NULL;
  char *inbam = NULL;
  char *acc2tax = NULL;
  int ctga =0;//only print ctga errors
  int search = -1;
  int countout = 0;
  char *infile_nodes = NULL;
  int howmany = 15;
  int doold = 0;
  while(*(++argv)){
    if(strcasecmp("-acc2tax",*argv)==0)
      acc2tax = strdup(*(++argv));
    else if(strcasecmp("-bam",*argv)==0)
      inbam = strdup(*(++argv));
    else if(strcasecmp("-r",*argv)==0)
      search = atoi(*(++argv));
    else if(strcasecmp("-howmany",*argv)==0)
      howmany = atoi(*(++argv));
    else if(strcasecmp("-ctga",*argv)==0)
      ctga =1;
    else if(strcasecmp("-doOld",*argv)==0)
      doold = 1;
    else if(strcasecmp("-countout",*argv)==0)
      countout =1;
    else if(strcasecmp("-nodes",*argv)==0)
      infile_nodes = strdup(*(++argv));
    else
      infile = strdup(*argv);
  }


  fprintf(stderr,"infile: %s inbam: %s names: %s search: %d ctga: %d countout: %d nodes: %s\n",infile,inbam,acc2tax,search,ctga,countout,infile_nodes);
  assert(infile);
  int2char name_map;
  if(acc2tax!=NULL)
    name_map = parse_names(acc2tax);

  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  if(infile_nodes!=NULL)
    parse_nodes(infile_nodes,rank,parent,child,1);
  if(search!=-1&& doold==0){
    std::map<int, double *> retmap = load_bdamage3(infile,howmany);
    double *dbl = getval(retmap,child,search,howmany);
    double dbldbl[3*howmany+1];//3 because ct,ga,other
    dbldbl[0] = dbl[0];
    for(int i=0;i<3*howmany;i++)
      dbldbl[i+1] = dbl[1+i]/dbl[0];
    
    fprintf(stdout,"%d\t%.0f",search,dbldbl[0]);
    for(int i=0;i<3*howmany;i++)
      fprintf(stdout,"\t%f",dbldbl[1+i]);
    fprintf(stdout,"\n");
    return 0;
  }
  
  BGZF *bgfp = NULL;
  samFile *samfp = NULL;
  bam_hdr_t *hdr =NULL;

  if(((bgfp = bgzf_open(infile, "r")))== NULL){
    fprintf(stderr,"Could not open input BAM file: %s\n",infile);
    return 1;
  }
  
  if(inbam!=NULL){
    if(((  samfp = sam_open_format(inbam, "r", NULL) ))== NULL){
      fprintf(stderr,"Could not open input BAM file: %s\n",inbam);
      return 1;
    }
    if(((hdr= sam_hdr_read(samfp))) == NULL){
      fprintf(stderr,"Could not read header for: %s\n",inbam);
    return 1;
    }
  }
  
  int printlength;
  assert(sizeof(int)==bgzf_read(bgfp,&printlength,sizeof(int)));
  fprintf(stderr,"\t-> printlength(howmany) from inputfile: %d\n",printlength);

  int ref_nreads[2];
  char *type_name = NULL;
  if(hdr!=NULL)
    type_name = strdup("#Reference");
  else if(acc2tax!=NULL)
    type_name = strdup("#FunkyName");
  else
    type_name = strdup("#taxid");
  
  if(ctga==0){
    fprintf(stdout,"%s\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n",type_name);
  }else{
    fprintf(stdout,"%s\tNalignment",type_name);
    for(int i=0;i<howmany;i++)
      fprintf(stdout,"\tCT_%d",i);
    for(int i=0;i<howmany;i++)
      fprintf(stdout,"\tGA_%d",i);
    fprintf(stdout,"\n");
  }
    
  int data[16];

  while(1){
    int nread=bgzf_read(bgfp,ref_nreads,2*sizeof(int));
    double ctgas[2*printlength];
    if(nread==0)
      break;
    fprintf(stderr,"ref: %d nreads: %d\n",ref_nreads[0],ref_nreads[1]);
    assert(nread==2*sizeof(int));
    for(int at=0;at<printlength;at++){
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      if((at+1)>howmany)
	continue;
      if(search==-1||search==ref_nreads[0]) {
	if(ctga==0)  {
	  if(hdr!=NULL)
	    fprintf(stdout,"%s\t%d\t5\'\t%d",hdr->target_name[ref_nreads[0]],ref_nreads[1],at);
	  else if(acc2tax!=NULL){
	    int2char::iterator itt = name_map.find(ref_nreads[0]);
	    if(itt==name_map.end()){
	      fprintf(stderr,"\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n",ref_nreads[0],acc2tax);
	      exit(0);
	    }
	    fprintf(stdout,"\"%s\"\t%d\t5\'\t%d",itt->second,ref_nreads[1],at);
	  }else
	    fprintf(stdout,"%d\t%d\t5\'\t%d",ref_nreads[0],ref_nreads[1],at);
	}else{
	  if(at==0)
	    fprintf(stdout,"%d\t%d",ref_nreads[0],ref_nreads[1]);
	}
	if(countout==1){
	  for(int i=0;i<16;i++)
	    fprintf(stdout,"\t%d",data[i]);
	  fprintf(stdout,"\n");
	}else{
	  float flt[16];
	  
	  for(int i=0;i<4;i++){
	    double tsum =0;
	    for(int j=0;j<4;j++){
	      tsum += data[i*4+j];
	      flt[i*4+j] = data[i*4+j];
	    }
	    if(tsum==0) tsum = 1;
	    for(int j=0;j<4;j++)
	      flt[i*4+j] /=tsum;
	  }
	  
	  if(ctga==0){
	    for(int j=0;j<16;j++)
	      fprintf(stdout,"\t%f",flt[j]);
	    fprintf(stdout,"\n");
	  }else
	    ctgas[at] = flt[7];
	}
      }
    }
    if(search==-1||search==ref_nreads[0]){
      if(ctga==1){
	for(int i=0;i<howmany;i++)
	  fprintf(stdout,"\t%f",ctgas[i]);
      }
    }
  
    for(int at=0;at<printlength;at++) {
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      if(at+1>howmany)
	continue;
      if(search==-1||search==ref_nreads[0]){
	if(ctga==0) {
	  if(hdr!=NULL)
	    fprintf(stdout,"%s\t%d\t3\'\t%d",hdr->target_name[ref_nreads[0]],ref_nreads[1],at);
	  else if(acc2tax!=NULL){
	    int2char::iterator itt = name_map.find(ref_nreads[0]);
	    if(itt==name_map.end()){
	      fprintf(stderr,"\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n",ref_nreads[0],acc2tax);
	      exit(0);
	    }
	    fprintf(stdout,"\"%s\"\t%d\t3\'\t%d",itt->second,ref_nreads[1],at);
	  }
	  else
	    fprintf(stdout,"%d\t%d\t3\'\t%d",ref_nreads[0],ref_nreads[1],at);
	}
	if(countout==1){
	  for(int i=0;i<16;i++)
	    fprintf(stdout,"\t%d",data[i]);
	  fprintf(stdout,"\n");
	}else{	
	  float flt[16];
	  
	  for(int i=0;i<4;i++){
	    double tsum =0;
	    for(int j=0;j<4;j++){
	      tsum += data[i*4+j];
	      flt[i*4+j] = data[i*4+j];
	    }
	    if(tsum==0) tsum = 1;
	    for(int j=0;j<4;j++)
	      flt[i*4+j] /=tsum;
	  }
	  if(ctga==0){
	    for(int j=0;j<16;j++)
	      fprintf(stdout,"\t%f",flt[j]);
	    fprintf(stdout,"\n");
	  }else
	  ctgas[at+printlength] = flt[8];
	}
      }
    }

    if(search==-1||search==ref_nreads[0]){
      if(ctga==1){
	for(int i=0;i<howmany;i++)
	  fprintf(stdout,"\t%f",ctgas[printlength+i]);
	fprintf(stdout,"\n");
      }
    }
  }

  if(bgfp)
    bgzf_close(bgfp);
  if(hdr)
    bam_hdr_destroy(hdr);
  if(samfp)
    sam_close(samfp);
  return 0;
}


int main_merge(int argc,char **argv){
  fprintf(stderr,"./metadamage merge file.lca file.bdamage.gz [-names file.gz -bam file.bam -howmany 5 -nodes trestructure.gz]\n");
  if(argc<=2)
    return 0;
  char *infile_lca = argv[1];
  char *infile_bdamage = argv[2];
  char *infile_nodes = NULL;
  
  char *inbam = NULL;
  char *acc2tax = NULL;
  int howmany = 5;
  ++argv;
  while(*(++argv)){
    if(strcasecmp("-names",*argv)==0)
      acc2tax = strdup(*(++argv));
    if(strcasecmp("-bam",*argv)==0)
      inbam = strdup(*(++argv));
    if(strcasecmp("-howmany",*argv)==0)
      howmany = atoi(*(++argv));
    if(strcasecmp("-nodes",*argv)==0)
      infile_nodes = strdup(*(++argv));
  }


  fprintf(stderr,"infile_lca: %s infile_bdamage: %s nodes: %s\n",infile_lca,infile_bdamage,infile_nodes);

  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  if(infile_nodes!=NULL)
    parse_nodes(infile_nodes,rank,parent,child,1);

  std::map<int, double *> retmap = load_bdamage3(infile_bdamage,howmany);
  //  fprintf(stderr,"retmap.size():%lu\n",retmap.size());
  int2char name_map;
  if(acc2tax!=NULL)
    name_map = parse_names(acc2tax);
  
  BGZF *bgfp = NULL;
  samFile *samfp = NULL;
  bam_hdr_t *hdr =NULL;

  if(inbam!=NULL){
    if(((bgfp = bgzf_open(inbam, "r")))== NULL){
      fprintf(stderr,"Could not open input BAM file: %s\n",inbam);
      return 1;
    }
  }
  
  if(inbam!=NULL){
    if(((  samfp = sam_open_format(inbam, "r", NULL) ))== NULL){
      fprintf(stderr,"Could not open input BAM file: %s\n",inbam);
      return 1;
    }
    if(((hdr= sam_hdr_read(samfp))) == NULL){
      fprintf(stderr,"Could not read header for: %s\n",inbam);
    return 1;
    }
  }

  gzFile fp = Z_NULL;
  fp = gzopen(infile_lca,"rb");
  assert(fp!=Z_NULL);
  char buf[4096];
  char orig[4096];
  gzgets(fp,buf,4096);//skipheader
  buf[strlen(buf)-1] = '\0';
  fprintf(stderr,"%s: VERSION:%s\n",buf,METADAMAGE_VERSION);
  float presize = retmap.size();
  while(gzgets(fp,buf,4096)){
    strncpy(orig,buf,4096);
    //    fprintf(stderr,"buf: %s\n",buf);
    char *tok=strtok(buf,"\t\n ");
    int taxid=atoi(strtok(NULL,":"));
    //    fprintf(stderr,"taxid: %d\n",taxid);
    double *dbl = getval(retmap,child,taxid,howmany);
    double dbldbl[3*howmany+1];//3 because ct,ga,other
    dbldbl[0] = dbl[0];
    for(int i=0;i<3*howmany;i++)
      dbldbl[i+1] = dbl[1+i]/dbl[0];

    orig[strlen(orig)-1] = '\0';
    fprintf(stdout,"%s\t%d:%.0f",orig,taxid,dbldbl[0]);
    for(int i=0;i<3*howmany;i++)
      fprintf(stdout,"\t%f",dbldbl[1+i]);
    fprintf(stdout,"\n");
  }
  float postsize=retmap.size();
  fprintf(stderr,"\t-> pre: %f post:%f grownbyfactor: %f\n",presize,postsize,postsize/presize);
  if(bgfp)
    bgzf_close(bgfp);
  if(hdr)
    bam_hdr_destroy(hdr);
  if(samfp)
    sam_close(samfp);
  if(fp!=Z_NULL)
    gzclose(fp);
  return 0;
}

int2int getlcadist(char *fname){
  int2int lcadist;
  int tlen =strlen(fname)+10; 
  char tmp[tlen];
  snprintf(tmp,tlen,"%sdist",fname);
  //  fprintf(stderr,"tmp: %s\n",tmp);
  FILE *fp = NULL;
  fp=fopen(tmp,"rb");
  assert(fp!=NULL);
  char buf[4096];
  while(fgets(buf,4096,fp)){
    int key = atoi(strtok(buf,"\t\n "));
    int val = atoi(strtok(NULL,"\t\n "));
    lcadist[key] = val;
  }
  fprintf(stderr,"\t-> Done reading: %lu entries from file: \'%s\'\n",lcadist.size(),tmp);
  return lcadist;
}



std::map<int,double *> getcsv(char *fname){
  std::map<int,double *> ret;
  FILE *fp = NULL;
  fp=fopen(fname,"rb");
  assert(fp!=NULL);
  char buf[4096];
  fgets(buf,4096,fp);//skip header
  while(fgets(buf,4096,fp)){
    int key = atoi(strtok(buf,"\t\n, "));
    double *valval = new double[2];
    valval[0] = atof(strtok(NULL,"\t\n, "));
    valval[1] = atof(strtok(NULL,"\t\n, "));
    ret[key] = valval;
  }
  fprintf(stderr,"\t-> Done reading: %lu entries from file: \'%s\'\n",ret.size(),fname);
  return ret;
}


int main_merge2(int argc,char **argv){
  fprintf(stderr,"./metadamage merge2 file.lca file.bdamage.gz christian.csv [-names file.gz -bam file.bam -howmany 5 -nodes trestructure.gz]\n");
  if(argc<=3)
    return 0;
  char *infile_lca = argv[1];
  char *infile_bdamage = argv[2];
  char *infile_christian = argv[3];
  char *infile_nodes = NULL;
  
  char *acc2tax = NULL;
  int howmany = 5;
  fprintf(stdout,"#./metadamage "); 
  for(int i=0;i<argc;i++)
    fprintf(stdout,"%s ",argv[i]);
  fprintf(stdout,"\n");
  ++argv;
  while(*(++argv)){
    if(strcasecmp("-names",*argv)==0)
      acc2tax = strdup(*(++argv));
    if(strcasecmp("-howmany",*argv)==0)
      howmany = atoi(*(++argv));
    if(strcasecmp("-nodes",*argv)==0)
      infile_nodes = strdup(*(++argv));
  }
 
  fprintf(stderr,"infile_lca: %s infile_bdamage: %s nodes: %s\n",infile_lca,infile_bdamage,infile_nodes);

  int2int lcadist = getlcadist(infile_lca);
  std::map<int,double *> chris = getcsv(infile_christian);
  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  if(infile_nodes!=NULL)
    parse_nodes(infile_nodes,rank,parent,child,1);

  std::map<int, double *> retmap = load_bdamage3(infile_bdamage,howmany);
  fprintf(stderr,"\t-> Number of entries in damage pattern file: %lu\n",retmap.size());
  
  int2char name_map;
  if(acc2tax!=NULL)
    name_map = parse_names(acc2tax);
  
  BGZF *bgfp = NULL;

  gzFile fp = Z_NULL;
  fp = gzopen(infile_lca,"rb");
  assert(fp!=Z_NULL);
  char buf[4096];
  char orig[4096];
  gzgets(fp,buf,4096);//skipheader
  buf[strlen(buf)-1] = '\0';
  fprintf(stdout,"%s: VERSION:%s\n",buf,METADAMAGE_VERSION);
  float presize = retmap.size();
  double *rawval = new double[2];
  rawval[0]=rawval[1]=-1.0;
  while(gzgets(fp,buf,4096)){
    strncpy(orig,buf,4096);
    //    fprintf(stderr,"buf: %s\n",buf);
    char *tok=strtok(buf,"\t\n ");
    int taxid=atoi(strtok(NULL,":"));
    //    fprintf(stderr,"taxid: %d\n",taxid);
    int nclass = -1;
    int2int::iterator itit = lcadist.find(taxid);
    if(itit!=lcadist.end())
      nclass = itit->second;
    double *valval = rawval;
    std::map<int,double *>::iterator ititit = chris.find(taxid);
    if(ititit!=chris.end())
      valval = ititit->second;
    double *dbl = getval(retmap,child,taxid,howmany);
    double dbldbl[3*howmany+1];//3 because ct,ga,other
    dbldbl[0] = dbl[0];
    for(int i=0;i<3*howmany;i++)
      dbldbl[i+1] = dbl[1+i]/dbl[0];

    //  orig[strlen(orig)-1] = '\0';
    //rollout results
    tok = strtok(orig,"\t\n ");
    fprintf(stdout,"%s\t%d:%.0f:%d:%f:%f\t",tok,taxid,dbldbl[0],nclass,valval[0],valval[1]);
    for(int i=0;i<3*howmany-1;i++)
      fprintf(stdout,"%f:",dbldbl[1+i]);
    fprintf(stdout,"%f",dbldbl[3*howmany]);
    while(((tok=strtok(NULL,"\t\n ")))){
      fprintf(stdout,"\t%s",tok);
    }
    fprintf(stdout,"\n");
    
  }
  float postsize=retmap.size();
  fprintf(stderr,"\t-> pre: %f post:%f grownbyfactor: %f\n",presize,postsize,postsize/presize);
  if(bgfp)
    bgzf_close(bgfp);
  if(fp!=Z_NULL)
    gzclose(fp);
  return 0;
}

int main_print_all(int argc,char **argv){
  fprintf(stderr,"./metadamage print_all file.bdamage.gz -names file.gz [-howmany 5] -nodes trestructure.gz -names fil.gz\n");
  if(argc<=2)
    return 0;
  char *infile_bdamage = NULL;
  char *infile_nodes = NULL;
  char *infile_names = NULL;
  int howmany = 5;
 
  while(*(++argv)){
    if(strcasecmp("-names",*argv)==0)
      infile_names = strdup(*(++argv));
    else if(strcasecmp("-howmany",*argv)==0)
      howmany = atoi(*(++argv));
    else if(strcasecmp("-nodes",*argv)==0)
      infile_nodes = strdup(*(++argv));
    else
      infile_bdamage = strdup(*argv);
  }


  fprintf(stderr,"infile_names: %s infile_bdamage: %s nodes: %s ",infile_names,infile_bdamage,infile_nodes);
  fprintf(stderr,"#VERSION:%s\n",METADAMAGE_VERSION);
  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  if(infile_nodes!=NULL)
    parse_nodes(infile_nodes,rank,parent,child,1);

  std::map<int, double *> retmap = load_bdamage3(infile_bdamage,howmany);
  fprintf(stderr,"\t-> number of entries in damage pattern file: %lu\n",retmap.size());
  int2char name = parse_names(infile_names);
  
  BGZF *bgfp = NULL;
  samFile *samfp = NULL;

  float presize = retmap.size();
  getval(retmap,child,1,howmany); //this will do everything
  float postsize=retmap.size();
  fprintf(stderr,"\t-> pre: %f post:%f grownbyfactor: %f\n",presize,postsize,postsize/presize);
  for(std::map<int, double *>::iterator it=retmap.begin();it!=retmap.end();it++){
    int taxid = it->first;
    double *dbl = it->second;
    char *myrank =NULL;
    char *myname = NULL;
    int2char::iterator itc=rank.find(taxid);
    if(itc!=rank.end())
      myrank=itc->second;
    itc=name.find(taxid);
    if(itc!=name.end())
      myname = itc->second;
    double dbldbl[3*howmany+1];//3 because ct,ga,other
    dbldbl[0] = dbl[0];
    for(int i=0;i<3*howmany;i++)
      dbldbl[i+1] = dbl[1+i]/dbl[0];
    if(dbldbl[0]>0){
      fprintf(stdout,"%d:\"%s\":\"%s\":%.0f",taxid,myname,myrank,dbldbl[0]);
      for(int i=0;i<3*howmany;i++)
	fprintf(stdout,"\t%f",dbldbl[1+i]);
      fprintf(stdout,"\n");
    };
  }
 
  
  if(bgfp)
    bgzf_close(bgfp);

  return 0;
}


int main_print_ugly(int argc,char **argv){
  fprintf(stderr,"./metadamage print_ugly file.bdamage.gz -names file.gz -nodes trestructure.gz -names fil.gz\n");
  if(argc<=2)
    return 0;
  char *infile_bdamage = NULL;
  char *infile_nodes = NULL;
  char *infile_names = NULL;
  int howmany;
 
  while(*(++argv)){
    if(strcasecmp("-names",*argv)==0)
      infile_names = strdup(*(++argv));
    else if(strcasecmp("-nodes",*argv)==0)
      infile_nodes = strdup(*(++argv));
    else
      infile_bdamage = strdup(*argv);
  }


  fprintf(stderr,"infile_names: %s infile_bdamage: %s nodes: %s ",infile_names,infile_bdamage,infile_nodes);
  fprintf(stderr,"#VERSION:%s\n",METADAMAGE_VERSION);
  //map of taxid -> taxid
  int2int parent;
  //map of taxid -> rank
  int2char rank;
  //map of parent -> child taxids
  int2intvec child;

  if(infile_nodes!=NULL)
    parse_nodes(infile_nodes,rank,parent,child,1);

  std::map<int, mydata> retmap = load_bdamage_full(infile_bdamage,howmany);
  fprintf(stderr,"\t-> number of entries in damage pattern file: %lu printlength(howmany):%d\n",retmap.size(),howmany);
  int2char name = parse_names(infile_names);
  
  BGZF *bgfp = NULL;
  samFile *samfp = NULL;

  float presize = retmap.size();
  getval_full(retmap,child,1,howmany); //this will do everything
  float postsize=retmap.size();
  fprintf(stderr,"\t-> pre: %f post:%f grownbyfactor: %f\n",presize,postsize,postsize/presize);
  for(std::map<int, mydata>::iterator it=retmap.begin();it!=retmap.end();it++){
    int taxid = it->first;
    mydata md = it->second;
    if(it->second.nreads==0)
      continue;
    char *myrank =NULL;
    char *myname = NULL;
    int2char::iterator itc=rank.find(taxid);
    if(itc!=rank.end())
      myrank=itc->second;
    itc=name.find(taxid);
    if(itc!=name.end())
      myname = itc->second;

    for(int i=0;i<howmany;i++){
      fprintf(stdout,"%d\t\"%s\"\t\"%s\"\t%d\t5'\t%d",taxid,myname,myrank,it->second.nreads,i);
      for(int ii=0;ii<16;ii++)
	fprintf(stdout,"\t%lu",it->second.fw[i*16+ii]);
      fprintf(stdout,"\n");
    }
    for(int i=0;i<howmany;i++){
      fprintf(stdout,"%d\t\"%s\"\t\"%s\"\t%d\t3'\t%d",taxid,myname,myrank,it->second.nreads,i);
      for(int ii=0;ii<16;ii++)
	fprintf(stdout,"\t%lu",it->second.bw[i*16+ii]);
      fprintf(stdout,"\n");
    }

  }
 
  
  if(bgfp)
    bgzf_close(bgfp);

  return 0;
}

//from ngsLCA.cpp
int main_lca(int argc, char **argv);
int main(int argc, char **argv){
  
  clock_t t=clock();
  time_t t2=time(NULL);
 
  if(argc==1){
    fprintf(stderr,"./metadamage getdamage file.bam\n");
    fprintf(stderr,"./metadamage mergedamage files.damage.*.gz\n");
    fprintf(stderr,"./metadamage index files.damage.gz\n");
    fprintf(stderr,"./metadamage merge files.lca files.bdamage.gz\n");
    fprintf(stderr,"./metadamage merge2 files.lca files.bdamage.gz christian.csv\n");
    fprintf(stderr,"./metadamage lca [many options]\n");
    fprintf(stderr,"./metadamage print bdamage.gz\n");
    fprintf(stderr,"./metadamage print2 [many options] bdamage.gz\n");
    fprintf(stderr,"./metadamage print_all [many options] bdamage.gz\n");
    fprintf(stderr,"./metadamage print_ugly [many options] bdamage.gz\n");
    return 0;
  }
  fprintf(stderr,"#");
  for(int i=0;i<argc;i++)
    fprintf(stderr,"%s ",argv[i]);
  fprintf(stderr,"\n");
  fflush(stderr);
  argc--;++argv;
  if(!strcmp(argv[0],"getdamage"))
    main_getdamage(argc,argv);
  if(!strcmp(argv[0],"index"))
    main_index(argc,argv);
  if(!strcmp(argv[0],"print"))
    main_print(argc,argv);
   if(!strcmp(argv[0],"print_all"))
    main_print_all(argc,argv);
   if(!strcmp(argv[0],"print_ugly"))
     main_print_ugly(argc,argv);
   if(!strcmp(argv[0],"print2"))
    main_print2(argc,argv);
  if(!strcmp(argv[0],"merge"))
    main_merge(argc,argv);
  if(!strcmp(argv[0],"merge2"))
    main_merge2(argc,argv);
  if(!strcmp(argv[0],"lca"))
    main_lca(argc,argv);

  free(dingding2);
  //  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  //fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2)); 
  return 0;
}
