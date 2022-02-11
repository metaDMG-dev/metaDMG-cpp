//gpl thorfinn@binf.ku.dk
//lei zhao. code in bottom

#include <cstring>
#include <cstdlib>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <ctime>
#include <getopt.h>
#include <cmath>

#include "profile.h"


htsFormat *dingding3 =(htsFormat*) calloc(1,sizeof(htsFormat));

int nproc = 0;//number of reads processed

extern char refToChar[256];

char revTable[256];
char qsToProp[256];

double C0 = 0.01;
double p = 0.3;
double ppi = 0.01;

double pmd_stat(char *seq,char *ref,int len,char *qs){
  //make revcom
  static char compref1[512];
  for (int i=0; i<len; i++)
    compref1[i] = revTable[ref[len-i-1]];
  
  //do flipcount
  int count[2] = {0,0};
  for(int i=0;i<len;i++){
    if(ref[i] !=seq[i])
      count[0] = count[0] + 1;
    if(compref1[i] !=seq[i])
      count[1] = count[1] + 1;
  }

  char *myref = NULL;
  if(count[0]>count[1])
    myref = compref1;
  else
    myref = ref;

  //do stat
  double llh1 = 0;
  double llh2 = 0;
  for (int i=0; i<len; i++){
    double eps = qsToProp[qs[i]];
    if (ref[i] == 'C'){
      double Dz = C0+p*pow(1-p,i+1);
      double pm1 = (1-ppi)*(1-eps)*(1-Dz)+ppi*eps/3*(1-Dz)+ppi/3*(1-eps)*Dz+(1-ppi)*eps/3*Dz;
      double pm2 = (1-ppi)*(1-eps)+ppi*eps/3;
      if (seq[i] == 'C'){
	llh1 += log(pm1);
	llh2 += log(pm2);
      }else{
	llh1 += log(1-pm1);
	llh2 += log(1-pm2);
      }
    }else if (ref[i] == 'G'){
      double Dz = C0+p*pow(1-p,len-i);
      double pm1 = (1-ppi)*(1-eps)*(1-Dz)+ppi*eps/3*(1-Dz)+ppi/3*(1-eps)*Dz+(1-ppi)*eps/3*Dz;
      double pm2 = (1-ppi)*(1-eps)+ppi*eps/3;
      if (seq[i] == 'G'){
	llh1 += log(pm1);
	llh2 += log(pm2);
      }else{
	llh1 += log(1-pm1);
	llh2 += log(1-pm2);
      }
    }
  }
  return llh1-llh2;
}


void wrapper(const bam1_t *b,const char * reconstructedReference,const std::vector<int> &  reconstructedReferencePos,const int & minQualBase, int MAXLENGTH,float **mm5p,float **mm3p,float incval,char myread[512],char myref[512]){
  const char *alphabetHTSLIB = "NACNGNNNTNNNNNNN";
  char refeBase;
  char readBase;
  int  qualBase;
  
    int i;

    int j=0;
    for(i=0;i<int(b->core.l_qseq);i++,j++){

	refeBase=toupper(reconstructedReference[j]);
	readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);

	qualBase=int(bam_get_qual(b)[i]);  
  
	if( refeBase == 'S'){ //don't care about soft clipped or indels	    
	    j--;
	    continue;
	}
	
	if( refeBase == 'I'){ //don't care about soft clipped or indels
	  //i--;
	  continue;
	}

	if(refeBase == 'D'){//deletion
	    //j++;
	    i--;
	    continue;
	}

	if(qualBase < minQualBase)
	    continue;
	
	if(refeBase == 'M'){//match
	    refeBase =  readBase;
	}

	//	refeBase = refToChar[refeBase];
	//readBase = refToChar[readBase];
	//	fprintf(stderr,"en: %c to: %c\n",refeBase,readBase);
	myread[i] = readBase;
	myref[i] = refeBase;
	/*
	if( refeBase!=4  && readBase!=4 ){
	    int dist5p=i;
	    int dist3p=b->core.l_qseq-1-i;
	    
	    if( bam_is_rev(b) ){
		refeBase=com[refeBase];
		readBase=com[readBase];
		//dist5p=int(al.QueryBases.size())-1-i;
		dist5p=int(b->core.l_qseq)-1-i;
		dist3p=i;
	    }

	    if(dist5p<MAXLENGTH)
	      mm5p[dist5p][toIndex[refeBase][readBase]] += incval;
	    if(dist3p<MAXLENGTH)
	      mm3p[dist3p][toIndex[refeBase][readBase]] += incval;
	    
	}
	*/
    }
}


void parse_sequencingdata(char *refName,char *fname,int mapped_only,int se_only,int mapq){
  char reconstructedRef[512];
  char myread[512];
  char myrefe[512];
  std::pair< kstring_t*, std::vector<int> >  mypair;
  kstring_t *kstr =new kstring_t;
  kstr->l=kstr->m=0;
  kstr->s=NULL;
  mypair.first = kstr;
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
  char myqscore[512];
  while(((ret=sam_read1(in,hdr,b)))>0) {
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

#if 0
    fprintf(stderr,"readid:%s len: %d\nREAD:\t\n",bam_get_qname(b),b->core.l_qseq);
    for(int i=0;i<b->core.l_qseq;i++)
      fprintf(stderr,"%c",seq_nt16_str[bam_seqi(bam_get_seq(b),i)]);
    fprintf(stderr,"\nQSCOre: \n");
    for(int i=0;i<b->core.l_qseq;i++)
      fprintf(stderr,"%c",33+bam_get_qual(b)[i]);
#endif
    for(int i=0;i<b->core.l_qseq;i++)
      myqscore[i] = bam_get_qual(b)[i];
    
    //then we simply write it to the output
    memset(reconstructedRef,0,512);
    memset(myread,'N',512);
    memset(myrefe,'N',512);
    reconstructRefWithPosHTS(b,mypair,reconstructedRef);
    wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
#if 0
    fprintf(stderr,"\nmyread:\n%.*s\nmyReference:\n%.*s\n",b->core.l_qseq,myread,b->core.l_qseq,myrefe);
    fprintf(stderr,"---read[%d]----\n",nproc);
#endif

    //do stat
    double pmdstat = pmd_stat(myread,myrefe,b->core.l_qseq,myqscore);
    fprintf(stderr,"%s\tPMD: %f\n",bam_get_qname(b),pmdstat);
  }
 
  bam_destroy1(b);
  hts_opt_free((hts_opt *)dingding3->specific);
  free(dingding3);
  sam_hdr_destroy(hdr);
  sam_close(in);
  free(fname);
   
}

int usage(FILE *fp,int val){
  fprintf(stderr,"./metadamage pmd [-T ref.fa -@ threads -a se_only -q minmapQ -v VERBOSE] file.bam\n");
  fprintf(stderr,"-a is an integer zero or one, indicating if paired end reads should be discarded\n");
  return 0;
}


int main_pmd(int argc, char **argv){
  //build revtable
  memset(revTable,'n',256);
  revTable['A'] = 'T';
  revTable['C'] = 'G';
  revTable['G'] = 'C';
  revTable['T'] = 'A';

  //build qscore to prob table
  for(int i=0;i<256;i++){
    double val = i;
    qsToProp[i] = pow(10.0,-val/10.0);
  }
  
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


/*
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;
double C0 = 0.01;
double p = 0.3;
double ppi = 0.01;
int main(){
    FILE* outfile = fopen("example.txt","wt");
    ifstream infile("test1.txt");
    string line;
    while (getline(infile, line))
    {
        istringstream iss(line);
        int pos, len;
        string ref, read, bq;
        if (!(iss >> pos >> len >> ref >> read >> bq)) { break;}
        //cout << pos << "\t" << len << "\t" << read << "\t" << ref << "\t" << bq << "\n";
        cout << len <<"\n";
        cout << read.size() << "\n";
        int count1 = 0;
        int count2 = 0;
        string compref1 = "";
        for (int i=0; i<len; i++){
            if (ref[i] == 'A'){
                compref1 = compref1+"T";
            }else if (ref[i] == 'C'){
                compref1 = compref1+"G";
            }else if (ref[i] == 'G'){
                compref1 = compref1+"C";
            }else if (ref[i] == 'T'){
                compref1 = compref1+"A";
            }
        }
        //cout<<compref1<<"\n";
        string compref(compref1);
        reverse(compref.begin(),compref.end());
 //       cout<<compref<<"\n";
        if (len>read.size()){
            ref = ref.substr(0,read.size());
            compref = compref.substr(len-read.size(),read.size());
            //compref = compref.substr(0,read.size());
            len = read.size();
        }
 //       cout<<compref<<"\n";
        for (int i=0; i<read.size(); i++){
            if (ref[i] != read[i]){
                count1 = count1 + 1;
            }
            if (compref[i] != read[i]){
                count2 = count2 + 1;
            }
        }
        if (count1 > count2){ref = compref;}
        double llh1 = 0;
        double llh2 = 0;
        for (int i=0; i<len; i++){
            int Q = bq[i]-'!';
            double eps = pow(10,-Q/10);
            if (ref[i] == 'C'){
                double Dz = C0+p*pow(1-p,i+1);
                double pm1 = (1-ppi)*(1-eps)*(1-Dz)+ppi*eps/3*(1-Dz)+ppi/3*(1-eps)*Dz+(1-ppi)*eps/3*Dz;
                double pm2 = (1-ppi)*(1-eps)+ppi*eps/3;
                if (read[i] == 'C'){
                    llh1 += log(pm1);
                    llh2 += log(pm2);
                }else{
                    llh1 += log(1-pm1);
                    llh2 += log(1-pm2);
                }
            }else if (ref[i] == 'G'){
                double Dz = C0+p*pow(1-p,len-i);
                double pm1 = (1-ppi)*(1-eps)*(1-Dz)+ppi*eps/3*(1-Dz)+ppi/3*(1-eps)*Dz+(1-ppi)*eps/3*Dz;
                double pm2 = (1-ppi)*(1-eps)+ppi*eps/3;
                if (read[i] == 'G'){
                    llh1 += log(pm1);
                    llh2 += log(pm2);
                }else{
                    llh1 += log(1-pm1);
                    llh2 += log(1-pm2);
                }
            }
        }
        cout<<ref<<"\n";
        cout<<read<<"\n";
        cout<<llh1-llh2<<"\n";
        fprintf(outfile,"%s\t%s\t%f\n",ref.c_str(),read.c_str(),llh1-llh2);
//        fprintf(outfile,"%f\n",llh1-llh2);
        
    }
    fclose(outfile);
    return 0;
}

 */
