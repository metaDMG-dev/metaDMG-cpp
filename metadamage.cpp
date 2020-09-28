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

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));
typedef std::map<int,char *> int2char;
int usage_getdamage(FILE *fp){
  fprintf(fp,"\nUsage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>\n");
  fprintf(fp,"\nExample: ./metadamage getdamage -l 10 -p 5 --threads 8 ../data/subs.sam\nOptions:\n");
  fprintf(fp,"  -f/--fasta\t is required with CRAM\n");
  fprintf(fp,"  -l/--minlength\t reads shorter than minlength will be discarded\n");
  fprintf(fp,"  -r/--runmode\trunmode 1 means that damage patterns will be calculated for each chr/scaffold contig.\n\t\trunmode 0 means one global estimate.\n");
  fprintf(fp,"  -@/--threads\t Number of threads used for reading/writing\n");
  return 0;
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
  char *onam = "meta";
  int nthreads = 4;
  //fix these
  static struct option lopts[] = {
    {"fasta", 1, 0, 'f'},
    {"minlength", 1, 0, 'l'},
    {"threads", 1, 0, '@'},
    {"length", 1, 0, 'p'},
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
    case 'o': onam = strdup(optarg); break;
    case 'r': runmode = atoi(optarg); break;
    case '?':
      return usage_getdamage(stdout);

    default:
      fprintf(stderr,"Never here\n",optarg,fname);
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
  while(((ret=sam_read1(fp,hdr,b)))>=0){
    if(bam_is_unmapped(b) ){
      	fprintf(stderr,"skipping: %s unmapped \n",bam_get_qname(b));
      continue;
    }
    if(bam_is_failed(b) ){
      fprintf(stderr,"skipping: %s failed: flags=%d \n",bam_get_qname(b),b->core.flag);
      continue;
    }
    if(b->core.l_qseq < minLength){
      fprintf(stderr,"skipping: %s too short \n",bam_get_qname(b));
      continue;
    }
    if(bam_is_paired(b)){
      fprintf(stderr,"skipping: %s  is paired (can be considered using the -paired flag\n",bam_get_qname(b));
      continue;
    }
    //    fprintf(stderr,"Analyzing\n");
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
}


int main_print(int argc,char **argv){
  fprintf(stderr,"./metadamage print file.bdamage.gz [-names file.gz -bam file.bam]\n");
  char *infile = argv[1];
  char *inbam = NULL;
  char *acc2tax = NULL;
  int search = -1;
  ++argv;
  while(*(++argv)){
    if(strcasecmp("-names",*argv)==0)
      acc2tax = strdup(*(++argv));
    if(strcasecmp("-bam",*argv)==0)
      inbam = strdup(*(++argv));
    if(strcasecmp("-r",*argv)==0)
      search = atoi(*(++argv));
  }


  fprintf(stderr,"infile: %s inbam: %s names: %s\n",infile,inbam,acc2tax);

  int2char name_map;
  if(acc2tax!=NULL)
    name_map = parse_names(acc2tax);

  
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
  fprintf(stderr,"printlength: %d\n",printlength);

  int ref_nreads[2];

  if(hdr!=NULL)
    fprintf(stdout,"#Reference\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
  else if(acc2tax!=NULL)
    fprintf(stdout,"#FunkyName\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
  else
    fprintf(stdout,"#taxid\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
  
  int data[16];
  while(1){
    int nread=bgzf_read(bgfp,ref_nreads,2*sizeof(int));
    if(nread==0)
      break;
    assert(nread==2*sizeof(int));
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      if(search==-1||search==ref_nreads[0]){
	if(hdr!=NULL)
	  fprintf(stdout,"%s\t%d\t5\'\t%d",hdr->target_name[ref_nreads[0]],ref_nreads[1],i);
	else if(acc2tax!=NULL){
	  int2char::iterator itt = name_map.find(ref_nreads[0]);
	  if(itt==name_map.end()){
	    fprintf(stderr,"\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n",ref_nreads[0],acc2tax);
	    exit(0);
	  }
	  fprintf(stdout,"%s\t%d\t5\'\t%d",itt->second,ref_nreads[1],i);
	}else
	  fprintf(stdout,"%d\t%d\t5\'\t%d",ref_nreads[0],ref_nreads[1],i);
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
	for(int j=0;j<16;j++)
	  fprintf(stdout,"\t%f",flt[j]);
	fprintf(stdout,"\n");
      }
    }
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      if(search==-1||search==ref_nreads[0]){
	if(hdr!=NULL)
	  fprintf(stdout,"%s\t%d\t3\'\t%d",hdr->target_name[ref_nreads[0]],ref_nreads[1],i);
	else if(acc2tax!=NULL){
	  int2char::iterator itt = name_map.find(ref_nreads[0]);
	  if(itt==name_map.end()){
	    fprintf(stderr,"\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n",ref_nreads[0],acc2tax);
	    exit(0);
	  }
	  fprintf(stdout,"%s\t%d\t3\'\t%d",itt->second,ref_nreads[1],i);
	}
	else
	  fprintf(stdout,"%d\t%d\t3\'\t%d",ref_nreads[0],ref_nreads[1],i);
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
	for(int j=0;j<16;j++)
	  fprintf(stdout,"\t%f",flt[j]);
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
}

double *getval(std::map<int, double *> &retmap,int2intvec &child,int taxid,int howmany){
  //  fprintf(stderr,"getval(taxid):%d howmany: %d\n",taxid,howmany);
  std::map<int,double *>::iterator it = retmap.find(taxid);
  if(it!=retmap.end()){
    //fprintf(stderr,"has found: %d in retmap\n",it->first);
    return it->second;
  }
  double *ret = new double [2*howmany];
  for(int i=0;i<2*howmany;i++)
    ret[i] = 0.0;
  if(child.size()>0){// if we have supplied -nodes
    int2intvec::iterator it2 = child.find(taxid);
    if (it2!=child.end()){
      std::vector<int> &avec = it2->second;
      int efsize =0;
      for(int i=0;i<avec.size();i++){
	//	fprintf(stderr,"%d/%d %d\n",i,avec.size(),avec[i]);
	double *tmp = getval(retmap,child,avec[i],howmany);
	int hasdata =0;
	for(int i=1;i<2*howmany;i++)
	  if(tmp[0]!=tmp[i]){
	    hasdata++;
	    break;
	  }
	for(int i=0;hasdata&&i<2*howmany;i++)
	  ret[i] += tmp[i];
	if(hasdata)
	  efsize++;
      }
      // fprintf(stderr,"efsize: %d aveclsize:%lu\n",efsize,avec.size());
      for(int i=0;(efsize>0)&&i<2*howmany;i++)
	ret[i] /= (1.0*efsize);
    }
  }
  
  retmap[taxid] = ret;

  return ret;
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






  std::map<int, double *> retmap = load_bdamage(infile_bdamage,5);
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
  float presize = retmap.size();
  while(gzgets(fp,buf,4096)){
    strncpy(orig,buf,4096);
    //    fprintf(stderr,"buf: %s\n",buf);
    char *tok=strtok(buf,"\t\n ");
    int taxid=atoi(strtok(NULL,":"));
    //    fprintf(stderr,"taxid: %d\n",taxid);
    double *dbl = getval(retmap,child,taxid,howmany);
    orig[strlen(orig)-1] = '\0';
    fprintf(stdout,"%s\t%d",orig,taxid);
    for(int i=0;i<2*howmany;i++)
      fprintf(stdout,"\t%f",dbl[i]);
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
    fprintf(stderr,"./metadamage lca [many options]\n");
    return 0;
  }
  argc--;++argv;
  if(!strcmp(argv[0],"getdamage"))
    main_getdamage(argc,argv);
  if(!strcmp(argv[0],"index"))
    main_index(argc,argv);
  if(!strcmp(argv[0],"print"))
    main_print(argc,argv);
  if(!strcmp(argv[0],"merge"))
    main_merge(argc,argv);
  if(!strcmp(argv[0],"lca"))
    main_lca(argc,argv);

  //  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  //fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2)); 
  return 0;
}
