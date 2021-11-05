//gpl thorfinn@binf.ku.dk
#include <cstring>
#include <cstdlib>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <ctime>
#include <getopt.h>
#include "profile.h"

htsFormat *dingding3 =(htsFormat*) calloc(1,sizeof(htsFormat));

int nproc = 0;//number of reads processed

void parse_sequencingdata(char *refName,char *fname,int mapped_only,int se_only,int mapq){
  char reconstructedRef[512];
  std::pair< kstring_t*, std::vector<int> >  reconstructedReference;
  kstring_t *kstr =new kstring_t;
  kstr->l=kstr->m=0;
  kstr->s=NULL;
  reconstructedReference.first = kstr;
  samFile *in=NULL;
  
  if(refName){
    char *ref =(char*) malloc(10 + strlen(refName) + 1);
    sprintf(ref, "reference=%s", refName);
    hts_opt_add((hts_opt **)&dingding3->specific,ref);
    free(ref);
  }
  
  if(strstr(fname,".cram")!=NULL && refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    exit(0);
  }
  
  if((in=sam_open_format(fname,"r",dingding3))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
		exit(0);
  }
  
  
  bam_hdr_t  *hdr = sam_hdr_read(in);
  
  bam1_t *b = bam_init1();
  
  int ret;
  int refId=-1;
  while(((ret=sam_read1(in,hdr,b)))>0){
    nproc++;

  
    //if -m was used, discard unmapped reads
    if(mapped_only!=0){
      if(b->core.flag&4)
	continue;
    }
    
    //default -a 1
    //only use single end reads
    //which is either single end or collapsed reads
    if(se_only==1){
      if(b->core.flag&1)//if paired continue
	continue;
    }
    
    
    // if mapq threshold is set
    if(mapq!=-1 && b->core.qual<mapq)
      continue;
    if(refId==-1||refId!=b->core.tid){
      refId=b->core.tid;
      fprintf(stderr,"\t-> Now at Chromosome: %s\n",hdr->target_name[refId]);
    }

    fprintf(stderr,"readid:%s len: %d\nREAD:\t\n",bam_get_qname(b),b->core.l_qseq);
    for(int i=0;i<b->core.l_qseq;i++)
      fprintf(stderr,"%c",seq_nt16_str[bam_seqi(bam_get_seq(b),i)]);
    fprintf(stderr,"\nQSCOre: \n");
    for(int i=0;i<b->core.l_qseq;i++)
      fprintf(stderr,"%c",33+bam_get_qual(b)[i]);
    
    //then we simply write it to the output
    
    memset(reconstructedRef,0,512);
    reconstructRefWithPosHTS( b,reconstructedReference,reconstructedRef);
    fprintf(stderr,"\nReference:%s\n",reconstructedRef);

  }
 
  bam_destroy1(b);
  hts_opt_free((hts_opt *)dingding3->specific);
  free(dingding3);
  sam_hdr_destroy(hdr);
  sam_close(in);
  free(fname);
   
}

int usage(FILE *fp,int val){
  fprintf(stderr,"./metadamage pmd [options]\n");
  return 0;
}


int main_pmd(int argc, char **argv){

  	int VERBOSE = 0;
	unsigned long int seed = 0;

	clock_t t=clock();
	time_t t2=time(NULL);

	char *fname,*refName;

	fname=refName=NULL;
	int c;
	int nthreads = 1;

	//thresholds for filters

	int mapq =-1;
	int mapped_only = 0;
	int se_only = 1;
	
	if(argc==1){
		usage(stdout,0);
		return 0;
	}
	static struct option lopts[] = {
		{"add", 1, 0, 0},
		{"append", 0, 0, 0},
		{"delete", 1, 0, 0},
		{"verbose", 0, 0, 0},
		{"create", 1, 0, 'c'},
		{"file", 1, 0, 0},
		{NULL, 0, NULL, 0}
	};

	// x: means x has a param; is not a switch
	while ((c = getopt_long(argc, argv,
					"T:@:a:q:m:v:?",
					lopts, NULL)) >= 0) {
		switch (c) {
		case 'T': refName = strdup(optarg); break;
		case '@': nthreads = atoi(optarg); break;
		case 'a': se_only = atoi(optarg); break;
		case 'q': mapq = atoi(optarg); break;
		case 'm': mapped_only = 1; break;
		case 'v': VERBOSE = 1; break;
		case '?':
		  if (optopt == '?') {  // '-?' appeared on command line
		    return usage(stdout,0);
		  } else {
		    if (optopt) { // Bad short option
		      fprintf(stdout,"./metadamage pmd invalid option -- '%c'\n", optopt);
		    } else { // Bad long option
		      // Do our best.  There is no good solution to finding
		      // out what the bad option was.
		      // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
		      if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
			fprintf(stdout,"./superduper unrecognised option '%s'\n",argv[optind - 1]);
		      }
		    }
		    return 0;//usage(stderr, 0);
		  }
		default:
		  fname = strdup(optarg);
		  fprintf(stderr,"assinging: %s to fname:%s\n",optarg,fname);
		  break;
		}
	}
	if(optind<argc)
	  fname = strdup(argv[optind]);
	
	if(!fname){
	  fprintf(stderr,"\t-> No input file specified\n");
	  usage(stdout,0);
	  return 0;
	}
	
	if(fname){
	  parse_sequencingdata(refName,fname,mapped_only,se_only,mapq);

	}

	fprintf(stderr,
			"\n"
			"\t[ALL done] cpu-time used =  %.2f sec\n"
			"\t[ALL done] walltime used =  %.2f sec\n"
			,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));  
	return 0;
}

