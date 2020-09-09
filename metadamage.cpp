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
#include "profile.h"

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

int usage(FILE *fp,int a){
  return 0;
}

int main_getdamage(int argc,char **argv){
  fprintf(stderr,"%s\n",__FUNCTION__);
  
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
    {"add", 1, 0, 0},
    {"append", 0, 0, 0},
    {"delete", 1, 0, 0},
    {"verbose", 0, 0, 0},
    {"create", 1, 0, 'c'},
    {"file", 1, 0, 0},
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
	  if (optopt == '?') {  // '-?' appeared on command line
	    return usage(stdout,0);
	  } else {
	    if (optopt) { // Bad short option
	      fprintf(stdout,"./metadamage invalid option -- '%c'\n", optopt);
	    } else { // Bad long option
	      // Do our best.  There is no good solution to finding
	      // out what the bad option was.
	      // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
	      if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
		fprintf(stdout,"./metadamage unrecognised option '%s'\n",argv[optind - 1]);
	      }
	    }
	    return 0;//usage(stderr, 0);
	  }
    default:
      fprintf(stderr,"Never here\n",optarg,fname);
      break;
    }
  }
  if(optind<argc)
    fname = strdup(argv[optind]);
  fprintf(stderr,"./metadamage refName: %s minLength: %d printLength: %d runmode: %d outname: %s nthreads: %d\n",refName,minLength,printLength,runmode,onam,nthreads);
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
  while(((ret=sam_read1(fp,hdr,b)))>0){
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
    dmg->damage_analysis(b,runmode!=0?b->core.tid:0);
  }


  dmg->printit(stdout,printLength);
  
  dmg->write(onam,runmode==1?hdr:NULL);
  dmg->bwrite(onam,hdr);
    
  sam_hdr_destroy(hdr);
  bam_destroy1(b);
  sam_close(fp);
  destroy_damage(dmg);

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
  fprintf(stderr,"./metadamage print file.bdamage.gz [file.bam]\n");
  char *infile = argv[1];
  char *inbam = NULL;
  if(argc>2)
    inbam = argv[2];
  fprintf(stderr,"infile: %s inbam: %s\n",infile,inbam);

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
      else
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
      else
	fprintf(stdout,"%d\t%d\t3\'\t%d",ref_nreads[0],ref_nreads[1],i);
      float flt[16];
      double tsum =0;
      for(int j=0;j<16;j++){
	tsum += data[j];
	flt[j] = data[j];
      }
      if(tsum==0) tsum = 1;
      for(int j=0;j<16;j++)
	fprintf(stdout,"\t%f",flt[j]/tsum);
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

int main(int argc, char **argv){
  clock_t t=clock();
  time_t t2=time(NULL);
 
  if(argc==1){
    fprintf(stderr,"./metadamage getdamage file.bam\n");
    fprintf(stderr,"./metadamage mergedamage files.damage.*.gz\n");
    fprintf(stderr,"./metadamage index files.damage.gz\n");
    fprintf(stderr,"./metadamage print files.damage.gz\n");
    return 0;
  }
  argc--;++argv;
  if(!strcmp(argv[0],"getdamage"))
    main_getdamage(argc,argv);
  if(!strcmp(argv[0],"index"))
    main_index(argc,argv);
  if(!strcmp(argv[0],"print"))
    main_print(argc,argv);

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2)); 
  return 0;
}
