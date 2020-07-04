#include <iostream>
#include <vector>
#include <cstring>
#include <ctype.h>
#include <string>
#include <utility>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <cstdlib>
#include <zlib.h>
#include <fstream>
#include <cmath>
#include <sstream>

#define MAXLENGTH 1000

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
char toBase[4] = {'A','C','G','T'};

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

using namespace std;
#define bam_is_sec(b)         (((b)->core.flag&BAM_FSECONDARY)        != 0)
#define bam_is_supp(b)        (((b)->core.flag&BAM_FSUPPLEMENTARY)    != 0)
#define bam_is_paired(b)      (((b)->core.flag&BAM_FPAIRED)     != 0)
#define bam_is_rmdup(b)       (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_qcfailed(b)    (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_unmapped(b)    (((b)->core.flag&BAM_FUNMAP)      != 0)
#define bam_is_read1(b)       (((b)->core.flag&BAM_FREAD1)      != 0)
#define bam_is_failed(b)      ( bam_is_qcfailed(b) || bam_is_rmdup(b) || bam_is_sec(b) || bam_is_supp(b) )

inline string printDoubleAsWhitePaddedString(double d,int digits,int precision){
    char buffer [digits+2+precision];
    sprintf(buffer,"%*.*f",digits+precision+1,precision, d);
    return string(buffer);
}

inline string printIntAsWhitePaddedString(int i,int digits){
    char buffer [digits+1];
    sprintf(buffer,"%*d",digits,i);
    return string(buffer);
}

inline string booleanAsString(bool toprint){
    if(toprint)
	return string("currently turned on/used");
    else
	return string("not on/not used");
}

template <typename T>
string stringify(const T i){
    stringstream s;
    s << i;
    return s.str();
}
	
template <typename T>
T destringify( const string& s ){
    istringstream i(s);
    T x;
    if (!(i >> x)){
	cerr<<"Utils.cpp: destringify() Unable to convert string=\""<<s<<"\""<<endl;
	exit(1);
    }
    return x;
} 

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

void  reconstructRefWithPosHTS(const bam1_t   * b,std::pair< std::string, std::vector<int> > &pp){
  pp.first = "";
  pp.second.clear();

  std::string reconstructedTemp="";
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

    for(int32_t i = 0; i < n_cigar_op; i++){
	char opchr = bam_cigar_opchr(cigar[i]);
        int32_t oplen = bam_cigar_oplen(cigar[i]);
	reconstructedTemp+=std::string(oplen,opchr);
	//cerr <<reconstructedTemp << endl;
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

    for(unsigned int i=0;i<reconstructedTemp.size();i++){
	if(reconstructedTemp[i] == 'M' ){ //only look at matches and indels	    
		
	    if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

		if(parsedMD[mdVectorIndex].offset == 0){ //we have reached a mismatch				

		    if(parsedMD[mdVectorIndex].bp == DUMMYCHAR){ //no char to add, need to backtrack on the CIGAR
			i--;
		    }else{
			pp.first += parsedMD[mdVectorIndex].bp;
			pp.second.push_back(initialPositionControl++);
		    }
		    mdVectorIndex++;
		}else{ //wait until we reach a mismatch
		    pp.first +=reconstructedTemp[i];
		    parsedMD[mdVectorIndex].offset--;
		    pp.second.push_back(initialPositionControl++);
		}

		//skipping all the positions with deletions on the read
		//if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

		while( (mdVectorIndex<int(parsedMD.size())) &&
		       (parsedMD[mdVectorIndex].bp == '^' ) ){ 
		    initialPositionControl+=parsedMD[mdVectorIndex].offset;
		    mdVectorIndex++;
		}
		    
	    }else{
		pp.first +=reconstructedTemp[i];
		pp.second.push_back(initialPositionControl++);
	    }
	}else{
	    if(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){ //soft clipped bases and indels
		pp.first +=reconstructedTemp[i];
		pp.second.push_back(initialPositionControl);
	    }
	}
    }

    if(int(pp.first.size()) != b->core.l_qseq){
	cerr << "Could not recreate the sequence for read "<<bam_get_qname(b)  << endl;
	exit(1);
    }

    if(pp.second.size() != pp.first.size()){
	cerr << "Could not determine the positions for the read "<<bam_get_qname(b) << endl;
	exit(1);
    }


    //    return pair< string, vector<int> >(reconstructed,positionsOnControl);
}


using namespace std;

const int offset=0;
int numberOfCycles;
string alphabetHTSLIB = "NACNGNNNTNNNNNNN";

#define MAXLENGTH 1000

vector< vector<unsigned int> > typesOfDimer5p; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p; //3' deam rates

//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(const   bam1_t  * b,const char * reconstructedReference,const vector<int> &  reconstructedReferencePos,const int & minQualBase, const bam_hdr_t *h,bool ispaired,bool isfirstpair){ // ,int firstCycleRead,int increment

    char refeBase;
    char readBase;
    int  qualBase;
 
    //Checking if the 5' is deaminated

    int i;
    if(ispaired){ //since we cannot evaluate the 5' ends or 3' ends
	goto iterateLoop;
    }
 iterateLoop:

    int j=0;
    for(i=0;i<int(b->core.l_qseq);i++,j++){

	refeBase=toupper(reconstructedReference[j]);
	readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);

	qualBase=int(bam_get_qual(b)[i])-offset;  
  
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


	//if( isResolvedDNA(refeBase)  && isResolvedDNA(readBase) ){
	if( refeBase!=4  && readBase!=4 ){
	  //fprintf(stderr,"BASES:%d %d\n",refeBase,readBase);
	    int dist5p=i;
	    int dist3p=b->core.l_qseq-1-i;
	    
	    if( bam_is_rev(b) ){
		refeBase=com[refeBase];
		readBase=com[readBase];
		//dist5p=int(al.QueryBases.size())-1-i;
		dist5p=int(b->core.l_qseq)-1-i;
		dist3p=i;
	    }

	    if(dist5p > MAXLENGTH ||
	       dist3p > MAXLENGTH ){
		cerr<<"Molecule found "<<bam_get_qname(b)<<" with length greater than limit"<<endl;
		exit(1);
	    }

	    if( !ispaired ||  isfirstpair){
	      //     fprintf(stderr,"increase5p: %d\n",dist5p);
	      typesOfDimer5p[dist5p][toIndex[refeBase][readBase]]++;
	    }

	    if( !ispaired || !isfirstpair){
	      //fprintf(stderr,"increase3p: %d\n",dist3p);
	      typesOfDimer3p[dist3p][toIndex[refeBase][readBase]]++;
	    }
	}
    }
}

int main(int argc, char *argv[]) {
  int lengthMaxToPrint = 5;
  int minQualBase      = 0;
  int minLength        = 35;
  
  bool paired=false;
  bool quiet=false;
  
  string usage=string(""+string(argv[0])+" <options>  [in BAM file]"+
			"\nThis program reads a BAM file and produces a deamination profile for the\n"+
			"5' and 3' ends\n"+

			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+

			"\n\n\tOther options:\n"+
			"\t\t"+"-minq\t\t\tRequire the base to have at least this quality to be considered (Default: "+stringify( minQualBase )+")\n"+
			"\t\t"+"-minl\t\t\tRequire the base to have at least this quality to be considered (Default: "+stringify( minLength )+")\n"+
			"\t\t"+"-length\t[length]\tDo not consider bases beyond this length  (Default: "+stringify(lengthMaxToPrint)+" ) \n"+
			"\t\t"+"-paired\t\t\tAllow paired reads    (Default: "+booleanAsString( paired )+" ) \n"+


			"\n\n\tYou can specify either one of the two:\n"+

			"\n\n\tOutput options:\n"+

			"\t\t"+"-q\t\t\tDo not print why reads are skipped (Default: "+booleanAsString(quiet)+")\n"+
		       
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "--help") )
    ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    
    for(int i=1;i<(argc-1);i++){ //all but the last 3 args
        if(string(argv[i]) == "-paired"  ){
            paired=true;
            continue;
        }

        if(string(argv[i]) == "-q"  ){
            quiet=true;
            continue;
        }

        if(string(argv[i]) == "-minq"  ){
            minQualBase=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-minl"  ){
            minLength=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-length"  ){
            lengthMaxToPrint=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    typesOfDimer5p       = vector< vector<unsigned int> >();
    typesOfDimer3p       = vector< vector<unsigned int> >();
        
    for(int l=0;l<MAXLENGTH;l++){
	//for(int i=0;i<16;i++){
	typesOfDimer5p.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p.push_back( vector<unsigned int> ( 16,0 ) );
    }
    
    string bamfiletopen = string( argv[ argc-1 ] );

    samFile  *fp;
    bam1_t    *b;
    bam_hdr_t *h;

    fp = sam_open_format(bamfiletopen.c_str(), "r", NULL); 
    if(fp == NULL){
	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
	return 1;
    }

    h = sam_hdr_read(fp);
    if(h == NULL){
        cerr<<"Could not read header for "<<bamfiletopen<<endl;
        return 1;
    }
    b = bam_init1();
    pair< string, vector<int> >  reconstructedReference;
    while(sam_read1(fp, h, b) >= 0){
	if(bam_is_unmapped(b) ){
	    if(!quiet)
		cerr<<"skipping "<<bam_get_qname(b)<<" unmapped"<<endl;
	    continue;
	}
	if(bam_is_failed(b) ){
	    if(!quiet)
		cerr<<"skipping "<<bam_get_qname(b)<<" failed"<<endl;
	    continue;
	}
	if(b->core.l_qseq < minLength){
	    if(!quiet)
		cerr<<"skipping "<<bam_get_qname(b)<<" too short"<<endl;
	    continue;
	}
	bool ispaired    = bam_is_paired(b);
	bool isfirstpair = bam_is_read1(b);
	if(!paired){    
	    if( ispaired   ){
		if(!quiet)
		    cerr<<"skipping "<<bam_get_qname(b)<<" is paired (can be considered using the -paired flag"<<endl;
		continue;
	    }
	}
	
	reconstructRefWithPosHTS(b,reconstructedReference);
	increaseCounters(b,reconstructedReference.first.c_str(), reconstructedReference.second,minQualBase,h,ispaired,isfirstpair); //start cycle numberOfCycles-1
    }
    
    bam_destroy1(b);
    sam_close(fp);
    
    cout <<"pos\t"<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;
  

    vector< vector<unsigned int> > * typesOfDimer5pToUse;

    typesOfDimer5pToUse     = &typesOfDimer5p;
    
    for(int l=0;l<lengthMaxToPrint;l++){
      cout<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";
	
	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer5pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		cout <<printDoubleAsWhitePaddedString( std::max(0.0,double( (*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)) ,1,5);
		
		if(!(n1 ==3 && n2 == 2 ))
		    cout <<"\t";
	    }


	}
	cout<<endl;
    }

    cout<<"pos\t"<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;


    vector< vector<unsigned int> > * typesOfDimer3pToUse;

    typesOfDimer3pToUse     = &typesOfDimer3p;
    
    for(int le=0;le<lengthMaxToPrint;le++){

	int l=le;
	
	cout <<""<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";	    
	
	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer3pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		cout<<printDoubleAsWhitePaddedString( std::max(0.0,double( (*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)) ,1,5);
		
		if(!(n1 ==3 && n2 == 2 ))
		    cout<<"\t";
	    }


	}
	cout<<endl;
    }

    return 0;
}
