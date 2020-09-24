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

//usefull little function to split
char *strpop(char **str,char split){
  char *tok=*str;
  while(**str){
    if(**str!=split)
      (*str)++;
    else{
      **str='\0'; (*str)++;
      break;
    }
  }
  return tok;
}
//usefull little function to remove tab and newlines
void strip(char *line){
  int at=0;
  //  fprintf(stderr,"%s\n",line);
  for(int i=0;i<strlen(line);i++)
    if(line[i]=='\t'||line[i]=='\n')
      continue;
    else
      line[at++]=line[i];
  //  fprintf(stderr,"at:%d line:%p\n",at,line);
  line[at]='\0';
  //fprintf(stderr,"%s\n",line);
}


 

int2char parse_names(const char *fname){

  gzFile gz= Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
    exit(0);
  }
  int2char name_map;
  char buf[4096];
  int at=0;
  char **toks = new char*[5];
  
  while(gzgets(gz,buf,4096)){
    strip(buf);//fprintf(stderr,"buf:%s\n",buf);
    char *saveptr = buf;
    toks[0]=strpop(&saveptr,'|');
    toks[1]= strpop(&saveptr,'|');
    toks[2]= strpop(&saveptr,'|');
    toks[3]= strpop(&saveptr,'|');
    for(int i=0;0&&i<4;i++)
      fprintf(stderr,"%d):\'%s\'\n",i,toks[i]);

    int key=atoi(toks[0]);
    //    fprintf(stderr,"key:%d\n",key);
    if(toks[3]&&strcmp(toks[3],"scientific name")==0){
      int2char::iterator it=name_map.find(key);
      
      if(it!=name_map.end())
	fprintf(stderr,"\t->[%s] duplicate name(column1): %s\n",fname,toks[0]);
      else
	name_map[key]=strdup(toks[1]);

    }
    if(0&&at++>10)
      break;
  }
  //  int2char::iterator it = name_map.find(61564);  assert(it!=name_map.end());
  fprintf(stderr,"\t-> [%s] Number of unique names (column1): %lu with third column 'scientific name'\n",fname,name_map.size());
  return name_map;
}


int main_print(int argc,char **argv){
  fprintf(stderr,"./metadamage print file.bdamage.gz [-names file.gz -bam file.bam]\n");
  char *infile = argv[1];
  char *inbam = NULL;
  char *acc2tax = NULL;
  ++argv;
  while(*(++argv)){
    if(strcasecmp("-names",*argv)==0)
      acc2tax = strdup(*(++argv));
    if(strcasecmp("-bam",*argv)==0)
      inbam = strdup(*(++argv));
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
    fprintf(stdout,"#Reference\tNreads\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
  else if(acc2tax!=NULL)
    fprintf(stdout,"#FunkyName\tNreads\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
  else
    fprintf(stdout,"#taxid\tNreads\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
  
  int data[16];
  while(1){
    int nread=bgzf_read(bgfp,ref_nreads,2*sizeof(int));
    if(nread==0)
      break;
    assert(nread==2*sizeof(int));
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
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
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
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

  if(bgfp)
    bgzf_close(bgfp);
  if(hdr)
    bam_hdr_destroy(hdr);
  if(samfp)
    sam_close(samfp);
}


int main_merge(int argc,char **argv){
  fprintf(stderr,"./metadamage merge file.lca file.bdamage.gz [-names file.gz -bam file.bam -howmany 5]\n");
  if(argc<=0)
    return 0;
  char *infile_lca = argv[1];
  char *infile_bdamage = argv[2];
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
  }


  fprintf(stderr,"infile_lca: %s infile_bdamage: %s\n",infile_lca,infile_bdamage);
  std::map<int, double *> retmap = load_bdamage(infile_bdamage,5);
  fprintf(stderr,"retmap.size():%lu\n",retmap.size());
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
  while(gzgets(fp,buf,4096)){
    strncpy(orig,buf,4096);
    //    fprintf(stderr,"buf: %s\n",buf);
    char *tok=strtok(buf,"\t\n ");
    int taxid=atoi(strtok(NULL,":"));
    //    fprintf(stderr,"taxid: %d\n",taxid);
    std::map<int,double *>::iterator it = retmap.find(taxid);
    orig[strlen(orig)-1] = '\0';
    fprintf(stdout,"%s\t%d",orig,it->first);
    if(it==retmap.end()){
      //      fprintf(stderr,"problem finding taxid: %d in damage analysis\n",taxid);
      for(int i=0;i<2*howmany;i++)
	fprintf(stdout,"\t0.0");
    }else{
      double *dbl=it->second;
      for(int i=0;i<2*howmany;i++)
	fprintf(stdout,"\t%f",dbl[i]);
    }
    fprintf(stdout,"\n");
  }
  
  if(bgfp)
    bgzf_close(bgfp);
  if(hdr)
    bam_hdr_destroy(hdr);
  if(samfp)
    sam_close(samfp);
  if(fp!=Z_NULL)
    gzclose(fp);
}

int main(int argc, char **argv){
  clock_t t=clock();
  time_t t2=time(NULL);
 
  if(argc==1){
    fprintf(stderr,"./metadamage getdamage file.bam\n");
    fprintf(stderr,"./metadamage mergedamage files.damage.*.gz\n");
    fprintf(stderr,"./metadamage index files.damage.gz\n");
    fprintf(stderr,"./metadamage merge files.lca files.bdamage.gz\n");
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

  //  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  //fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2)); 
  return 0;
}
