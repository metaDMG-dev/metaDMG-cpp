#include <vector>
#include <cstring>
#include <ctype.h>
#include <utility>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "profile.h"

unsigned **getmatrix(size_t x,size_t y){
  unsigned **ret = new unsigned*[x];
  for(int i=0;i<x;i++){
    ret[i] = new unsigned[y];
    for(int j=0;j<16;j++)
      ret[i][j] = 0;
  }
  return ret;
}

void destroymatrix(unsigned**d,size_t x){
  for(int i=0;i<x;i++)
    delete [] d[i];
  delete [] d;
}

void destroy_damage(damage *dmg){
  for(std::map<int,triple >::iterator it=dmg->assoc.begin();it!=dmg->assoc.end();it++){
    destroymatrix(it->second.mm5p,dmg->MAXLENGTH);
    destroymatrix(it->second.mm3p,dmg->MAXLENGTH);
  }
  free(dmg->reconstructedReference.first->s);
  
  delete dmg->reconstructedReference.first;
  delete dmg;
}

BGZF *my_bgzf_open(const char *name,int nthreads){
  BGZF *ret = NULL;
  ret = bgzf_open(name,"wb");
  assert(ret != NULL);
  if(nthreads>1){
    fprintf(stderr,"\t-> Setting threads to: %d \n",nthreads);
    bgzf_mt(ret,nthreads,64);
  }
  
  return ret;
}


//A=0,C=1,G=2,T=3
char refToChar[256] = {
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

char toIndex[4][4]={
  {0,1,2,3},
  {4,5,6,7},
  {8,9,10,11},
  {12,13,14,15}
};
//a->t,c->g,g->c,t->a
char com[4] = {3,2,1,0};

typedef struct{
    char bp;
    int offset;
} mdField;

static char DUMMYCHAR='#';

void mdString2Vector2(const uint8_t *md,std::vector<mdField> &toReturn){
  const char *mdFieldToParse =(const char*) md+1;
  toReturn.clear();
    int i=0;
    // int addToOffset=0;
    mdField toadd;

    toadd.offset=0;
    toadd.bp='N';

    while(strlen(mdFieldToParse) != i){
	if(isdigit(mdFieldToParse[i])){
	    toadd.offset=toadd.offset*10+mdFieldToParse[i]-'0';
	}else{
	    //deletions in read (insertion in reference)
	    if(mdFieldToParse[i] == '^'){
		if(toadd.offset != 0){
		    toadd.bp=DUMMYCHAR;
		    toReturn.push_back(toadd);
		    toadd.offset=0;
		    toadd.bp='N';
		}

		i++;
		mdField toadd2;
		toadd2.offset=0;
		toadd2.bp='^';
		while(isalpha(mdFieldToParse[i])){
		    i++;
		    toadd2.offset++;
		}
		toReturn.push_back(toadd2);
		i--;
	    }else{
		toadd.bp=mdFieldToParse[i];
		toReturn.push_back(toadd);

		toadd.offset=0;
		toadd.bp='N';
	    }

	}
	i++;
    }
}

void  reconstructRefWithPosHTS(const bam1_t   * b,std::pair< kstring_t *, std::vector<int> > &pp,char *reconstructedTemp){
  
  pp.first->l = 0;
  pp.second.clear();

  std::vector<mdField> parsedMD;
    //skip unmapped
    if( ((b)->core.flag&BAM_FUNMAP) != 0 ){
      fprintf(stderr,"The function reconstructRefWithPosOnReadHTS()  cannot be called for unmapped reads\n");
      exit(1);
    }
    
    uint8_t *mdptr = bam_aux_get(b, "MD");
    //    fprintf(stderr,"%s\n",bam_get_qname(b));
    if(mdptr==NULL){
      fprintf(stderr,"ReconsReferenceHTSLIB: Cannot get MD tag from:%s ",bam_get_qname(b));
      exit(1);
    }
    
    int32_t   n_cigar_op = b->core.n_cigar;
    uint32_t *cigar      = bam_get_cigar(b);

    int at =0;
    for(int32_t i = 0; i < n_cigar_op; i++){
	char opchr = bam_cigar_opchr(cigar[i]);
        int32_t oplen = bam_cigar_oplen(cigar[i]);
	memset(reconstructedTemp+at,opchr,oplen);
	at += oplen;
    }
    
    //get a vector representation of the MD field	
    mdString2Vector2(mdptr,parsedMD);
#if 0
    for(int i=0;i<parsedMD.size();i++)
      fprintf(stderr,"%d) %c %d\n",i,parsedMD[i].bp,parsedMD[i].offset);
#endif
        
    //int initialPositionControl=al->Position;
    int initialPositionControl=b->core.pos;

    //combine the CIGAR and MD into one single string
    int mdVectorIndex=0;

    for(unsigned int i=0;i<strlen(reconstructedTemp);i++){
	if(reconstructedTemp[i] == 'M' ) { //only look at matches and indels	    
		
	  if(mdVectorIndex < (int)parsedMD.size() ){ //still have mismatches

		if(parsedMD[mdVectorIndex].offset == 0){ //we have reached a mismatch				

		    if(parsedMD[mdVectorIndex].bp == DUMMYCHAR){ //no char to add, need to backtrack on the CIGAR
			i--;
		    }else{
		      kputc(parsedMD[mdVectorIndex].bp,pp.first);
		      pp.second.push_back(initialPositionControl++);
		    }
		    mdVectorIndex++;
		}else{ //wait until we reach a mismatch
		  kputc(reconstructedTemp[i],pp.first);
		  parsedMD[mdVectorIndex].offset--;
		  pp.second.push_back(initialPositionControl++);
		}

		while( (mdVectorIndex<(int)parsedMD.size()) && (parsedMD[mdVectorIndex].bp == '^' ) ){ 
		    initialPositionControl+=parsedMD[mdVectorIndex].offset;
		    mdVectorIndex++;
		}
		    
	    }else{
	      kputc(reconstructedTemp[i],pp.first);
	      pp.second.push_back(initialPositionControl++);
	    }
	}else{
	    if(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){ //soft clipped bases and indels
	      kputc(reconstructedTemp[i],pp.first);
	      pp.second.push_back(initialPositionControl);
	    }
	}
    }

    if(strlen(pp.first->s) != b->core.l_qseq){
      fprintf(stderr,"Could not recreate the sequence for read: %s pp.first->s: %s strlen():%lu\n",bam_get_qname(b),pp.first->s,strlen(pp.first->s));
      exit(1);
    }

    if(pp.second.size() != strlen(pp.first->s)){
      fprintf(stderr,"Could not determine the positions for the read: %s\n",bam_get_qname(b));
	exit(1);
    }
}

inline void increaseCounters(const bam1_t *b,const char * reconstructedReference,const std::vector<int> &  reconstructedReferencePos,const int & minQualBase, int MAXLENGTH,unsigned **mm5p,unsigned **mm3p){
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

	refeBase = refToChar[refeBase];
	readBase = refToChar[readBase];

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
	      mm5p[dist5p][toIndex[refeBase][readBase]]++;
	    if(dist3p<MAXLENGTH)
	      mm3p[dist3p][toIndex[refeBase][readBase]]++;
	    
	}
    }
}

int damage::damage_analysis(bam1_t *b,int which){
  if(assoc.find(which)==assoc.end()){
    triple val={0,getmatrix(MAXLENGTH,16),getmatrix(MAXLENGTH,16)};
    assoc[which] = val;
    mm5p = val.mm5p;
    mm3p = val.mm3p;
    //    fprintf(stderr,"has added which: %d\n",which);
  }
  std::map<int,triple >::iterator it=assoc.find(which);
  it->second.nreads++;
  if(b->core.l_qseq-10>temp_len){
    temp_len = b->core.l_qseq;
    kroundup32(temp_len);
    free(reconstructedTemp);
    reconstructedTemp =(char*) calloc(temp_len,1);
  }
  memset(reconstructedTemp,0,temp_len);
  reconstructRefWithPosHTS(b,reconstructedReference,reconstructedTemp);
  increaseCounters(b,reconstructedReference.first->s, reconstructedReference.second,minQualBase,MAXLENGTH,it->second.mm5p,it->second.mm3p);
  return 0;
}
void damage::write(char *fname,bam_hdr_t *hdr){
  //fprintf(stderr,"Dumping asso.size(): %lu\n",assoc.size());
  char *outname=strdup("metaout");
  if(fname){
    free(outname);
    outname = fname;
  }
  kstring_t kstr;
  kstr.l=kstr.m=0;
  kstr.s=NULL;
  char onam[1024];
  snprintf(onam,1024,"%s.res.gz",outname);
  fprintf(stderr,"\t-> Will dump: \'%s\' this contains damage patterns for: %lu items\n",onam,assoc.size());
  BGZF *fp= my_bgzf_open(onam,nthreads);

  for(std::map<int,triple>::iterator it=assoc.begin();it!=assoc.end();it++ ){
    if(it->second.nreads==0)//should never happen
      continue;
    if(hdr!=NULL)
      ksprintf(&kstr,"%s\t%lu",hdr->target_name[it->first],it->second.nreads);
    else
      ksprintf(&kstr,"%lu",it->second.nreads);
    for(int l=0;l<MAXLENGTH;l++){
      for(int i=0;i<16;i++)
	ksprintf(&kstr,"\t%d",it->second.mm5p[l][i]);
    }
    for(int l=0;l<MAXLENGTH;l++){
      for(int i=0;i<16;i++)
	ksprintf(&kstr,"\t%d",it->second.mm3p[l][i]);
    }
    ksprintf(&kstr,"\n");
    assert(bgzf_write(fp,kstr.s,kstr.l)==kstr.l);
    kstr.l=0;
  }
  bgzf_close(fp);
  free(kstr.s);
}
void damage::bwrite(char *fname,bam_hdr_t *hdr){
  //fprintf(stderr,"Dumping asso.size(): %lu\n",assoc.size());
  
  char onam[1024];
  snprintf(onam,1024,"%s.bdamage.gz",fname);
  fprintf(stderr,"\t-> Will dump: \'%s\' this contains damage patterns for: %lu items\n",onam,assoc.size());
  BGZF *fp= my_bgzf_open(onam,nthreads);
  assert(bgzf_write(fp,&MAXLENGTH,sizeof(int))==sizeof(int));
 
  for(std::map<int,triple>::iterator it=assoc.begin();it!=assoc.end();it++ ){
    if(it->second.nreads==0)//should never happen
      continue;
    assert(bgzf_write(fp,&it->first,sizeof(int))==sizeof(int));
    assert(bgzf_write(fp,&it->second.nreads,sizeof(int))==sizeof(int));
    for(int l=0;l<MAXLENGTH;l++)
      assert(bgzf_write(fp,it->second.mm5p[l],sizeof(int)*16)==sizeof(int)*16);

    for(int l=0;l<MAXLENGTH;l++)
      assert(bgzf_write(fp,it->second.mm3p[l],sizeof(int)*16)==sizeof(int)*16);
  }
  bgzf_close(fp);
}


int printinfo(FILE *fp){
  fprintf(fp,"./profile <options>  [in BAM file]\nThis program reads a BAM file and produces a deamination profile for the 5' and 3' ends\n");
  fprintf(fp,"Other options:\n\t-minq INT\n\t-minl INT\n\t-length INT\n\t-paired\n");
  return 0;
}


int printresults_grenaud2(FILE *fp,unsigned **mm5p,int lengthMaxToPrint){
  fprintf(fp,"pos\tA>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n");
  for(int l=0;l<lengthMaxToPrint;l++){
    fprintf(fp,"%*d\t",int(log10(lengthMaxToPrint))+1,l);
    for(int n1=0;n1<4;n1++){   
      int totalObs=0;
      for(int n2=0;n2<4;n2++)
	totalObs+=mm5p[l][4*n1+n2];
      
      for(int n2=0;n2<4;n2++){   
	if(n1==n2)
	  continue;
	fprintf(fp,"%*.*f",1+5+1,5,std::max(0.0,double(mm5p[l][4*n1+n2])/double(totalObs)));
	if(!(n1 ==3 && n2 == 2 ))
	  fprintf(fp,"\t");
      }
    }
    fprintf(fp,"\n");
  }
  return 0;
}

void damage::printit(FILE *fp,int l){
  if(mm5p)
    printresults_grenaud2(fp,mm5p,l);
  if(mm3p)
    printresults_grenaud2(fp,mm3p,l);
}

#ifdef __WITH_MAIN__
int main(int argc, char *argv[]) {
  int MAXLENGTH = 1000;
  int lengthMaxToPrint = 5;
  int minQualBase      = 0;
  int minLength        = 35;
  int quiet=0;
  
  if(argc == 1 ||(argc == 2 && (strcasecmp(argv[1],"--help")==0) )){
    printinfo(stderr);
    return 0;       
  }

  for(int i=1;i<(argc-1);i++){ //all but the last 3 args
    if(strcasecmp(argv[i],"-q")==0){
      quiet=1;
      continue;
    }
    if(strcasecmp(argv[i],"-minq")==0){
      minQualBase=atoi(argv[i+1]);
      i++;
      continue;
    }
    if(strcasecmp(argv[i],"-minl")==0){
      minLength=atoi(argv[i+1]);
      i++;
      continue;
    }
    if(strcasecmp(argv[i],"-length")==0){
      lengthMaxToPrint=atoi(argv[i+1]);
      i++;
      continue;
    }
    fprintf(stderr,"Error: unknown option: %s\n",argv[i]);
    return 1;
  }
   
  damage *dmg = init_damage(MAXLENGTH);
  char *bamfiletopen =  argv[ argc-1 ];

  ;
  bam1_t    *b;
  bam_hdr_t *h;
  
  if(((  fp = sam_open_format(bamfiletopen, "r", NULL) ))== NULL){
    fprintf(stderr,"Could not open input BAM file: %s\n",bamfiletopen);
    return 1;
  }
  
  if(((h= sam_hdr_read(fp))) == NULL){
    fprintf(stderr,"Could not read header for: %s\n",bamfiletopen);
    return 1;
  }
  b = bam_init1();
  
  while(sam_read1(fp, h, b) >= 0){
    if(bam_is_unmapped(b) ){
      if(!quiet)
	fprintf(stderr,"skipping: %s unmapped \n");
      continue;
    }
    if(bam_is_failed(b) ){
      if(!quiet)
	fprintf(stderr,"skipping: %s failed \n");
      continue;
    }
    if(b->core.l_qseq < minLength){
      if(!quiet)
	fprintf(stderr,"skipping: %s too short \n");
      continue;
    }
    if(bam_is_paired(b)){
      if(!quiet)
	fprintf(stderr,"skipping: %s  is paired (can be considered using the -paired flag\n",bam_get_qname(b));
      continue;
    }
    
    dmg->damage_analysis(b,0);
    
  }
  
  sam_hdr_destroy(h);
  bam_destroy1(b);
  sam_close(fp);
  fprintf(stderr,"nreads: %lu\n",dmg->assoc.begin()->second.nreads);
  printresults_grenaud2(stdout,dmg->mm5p,lengthMaxToPrint);
  printresults_grenaud2(stdout,dmg->mm3p,lengthMaxToPrint);
  for(int i=0;0&i<16;i++)
    fprintf(stdout,"%lu\t",dmg->mm5p[0][i]);
  destroy_damage(dmg);
  return 0;
}

#endif

std::map<int,double *> load_bdamage3(const char* fname,int howmany ){
  //  fprintf(stderr,"./metadamage print file.bdamage.gz [-names file.gz -bam file.bam]\n");
  const char *infile = fname;
  //  fprintf(stderr,"infile: %s howmany: %d \n",infile,howmany);
  
  BGZF *bgfp = NULL;

  if(((bgfp = bgzf_open(infile, "r")))== NULL){
    fprintf(stderr,"Could not open input BAM file: %s\n",infile);
    exit(0);
  }

  std::map<int,double*> retmap;
  
  int printlength;
  assert(sizeof(int)==bgzf_read(bgfp,&printlength,sizeof(int)));
  
  if(howmany>printlength){
    fprintf(stderr,"\t-> Problem binary file has data for: %d positions, but you are requesting merge with: %d positions \n",printlength,howmany);
    fprintf(stderr,"\t-> Solutions set -howmany to lower value\n");
    exit(0);
  }
  
  int ref_nreads[2];
 
  int data[16];
  while(1){
    int nread=bgzf_read(bgfp,ref_nreads,2*sizeof(int));
    if(nread==0)
      break;
    assert(nread==2*sizeof(int));
    
    double *formap = new double [1+3*howmany];
    for(int i=0;i<1+3*howmany;i++)
      formap[i] = 0.0;
    //    fprintf(stderr,"formap: %p\n",formap);
    int incer =0;
    formap[incer++] = ref_nreads[1];

    for(int at=0;at<printlength;at++) {
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      float flt[16];//this will contain the float representation of counts
      for(int i=0;i<4;i++) {//loop over A*,C*,G*,T*
	double tsum =0;
	for(int j=0;j<4;j++){//tsum is the sum of *A,*C,*G,*T
	  tsum += data[i*4+j];
	  flt[i*4+j] = data[i*4+j];
	}
	if(tsum==0)
	  tsum = 1;
	for(int j=0;j<4;j++)
	  flt[i*4+j] /=tsum; //now rescale such that *A+*C+*G+*T=1
      }
      if(at<howmany) {//carefull satan
	formap[incer++] =flt[7]*formap[0];//flt[7] := CT

	for(int n=0;n<4;n++)
	  for(int j=0;j<4;j++)
	    if(n!=j)
	      formap[1+2*howmany+at] += flt[n*4+j];//we add all mutations after format[1+ct+ga]
	formap[1+2*howmany+at] -= flt[7];//we dont want to add the ct, so we subtract
      }
    }
    for(int at=0;at<printlength;at++){
      assert(16*sizeof(int)==bgzf_read(bgfp,data,sizeof(int)*16));
      
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
      if(at<howmany) {//carefull satan
	formap[incer++] =flt[8]*formap[0];

	for(int n=0;n<4;n++)
	  for(int j=0;j<4;j++)
	    if(n!=j)
	      formap[1+2*howmany+at] += flt[n*4+j];//we add all mutations after format[1+ct+ga]
	formap[1+2*howmany+at] -= flt[8];//we dont want to add the ga, so we subtract
	formap[1+2*howmany+at] = formap[1+2*howmany+at]/10.0*formap[0];//previously we were dividing with 22?, I think it should be 10,
      }
    }
    for(int i=0;0&&i<2*howmany;i++)
      fprintf(stdout,"[%d] %f\n",i,formap[i]);
    retmap[ref_nreads[0]] = formap;
  }
  //  exit(0);
  if(bgfp)
    bgzf_close(bgfp);
  fprintf(stderr,"\t-> Done loading binary bdamage.gz file. It contains: %lu\n",retmap.size());
  for(std::map<int,double *>::iterator it = retmap.begin();0&&it!=retmap.end();it++)
    fprintf(stderr,"it->second:%p\n",it->second);

  return retmap;
}

std::map<int, mydata> load_bdamage_full(const char* fname,int &printlength){
  //  fprintf(stderr,"./metadamage print file.bdamage.gz [-names file.gz -bam file.bam]\n");
  const char *infile = fname;
  //  fprintf(stderr,"infile: %s howmany: %d \n",infile,howmany);
  
  BGZF *bgfp = NULL;

  if(((bgfp = bgzf_open(infile, "r")))== NULL){
    fprintf(stderr,"Could not open input BAM file: %s\n",infile);
    exit(0);
  }

  std::map<int,mydata> retmap;
  printlength =0;
  assert(sizeof(int)==bgzf_read(bgfp,&printlength,sizeof(int)));

  int ref_nreads[2];
 
  while(1){
    int nread=bgzf_read(bgfp,ref_nreads,2*sizeof(int));
    if(nread==0)
      break;
    assert(nread==2*sizeof(int));
    mydata md;
    
    md.fw = new size_t[16*printlength];
    md.bw = new size_t[16*printlength];
    md.nreads = ref_nreads[1];

    int tmp[16];
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(int)==bgzf_read(bgfp,tmp,sizeof(int)*16));
      for(int ii=0;ii<16;ii++)
	md.fw[i*16+ii] = tmp[ii];
    }
  
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(int)==bgzf_read(bgfp,tmp,sizeof(int)*16));
      for(int ii=0;ii<16;ii++)
	md.bw[i*16+ii] = tmp[ii];
    }
    retmap[ref_nreads[0]] = md;
  }

  if(bgfp)
    bgzf_close(bgfp);
  fprintf(stderr,"\t-> Done loading binary bdamage.gz file. It contains: %lu\n",retmap.size());
  for(std::map<int,mydata>::iterator it = retmap.begin();0&&it!=retmap.end();it++)
    fprintf(stderr,"it->second:%p\n",it->second);

  return retmap;
}
