#include <iostream>
#include <vector>
#include <cstring>
#include <cstring>
#include <set>
#include <ctype.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <string>
#include <utility>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <cstdlib>
#include <zlib.h>
#include <map>
#include <fstream>
#include <cmath>
#include <sstream>

#define MAXLENGTH 1000

typedef struct{
    char bp;
    int offset;
} mdField;

static int asciiOffsetZero=48;
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
	


//Returns an index for every 2mer of different
inline int twoBases2index(const char c1,const char c2){
    char _c1= toupper(c1);
    char _c2= toupper(c2);

    if(_c1     ==    'A'){

	if(_c2 ==    'A')
	    return 0;
	if(_c2 ==    'C')
	    return 1;
	if(_c2 ==    'G')
	    return 2;
	if(_c2 ==    'T')
	    return 3;

	//cerr<<"Utils.h:1 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(_c1     ==    'C'){


	if(_c2 ==    'A')
	    return 4;
	if(_c2 ==    'C')
	    return 5;
	if(_c2 ==    'G')
	    return 6;
	if(_c2 ==    'T')
	    return 7;

	//cerr<<"Utils.h:2 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(_c1     ==    'G'){

	if(_c2 ==    'A')
	    return 8;
	if(_c2 ==    'C')
	    return 9;
	if(_c2 ==    'G')
	    return 10;
	if(_c2 ==    'T')
	    return 11;


	//cerr<<"Utils.h:3 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    if(_c1     ==    'T'){

	if(_c2 ==    'A')
	    return 12;
	if(_c2 ==    'C')
	    return 13;
	if(_c2 ==    'G')
	    return 14;
	if(_c2 ==    'T')
	    return 15;

	//	cerr<<"Utils.h:4 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    //    cerr<<"Utils.h:5 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
    exit(1);
}


inline char complement(const char c){
    if(c ==    'A')
	return 'T';

    if(c ==    'C')
	return 'G';

    if(c ==    'G')
	return 'C';

    if(c ==    'T')
	return 'A';



    if(c ==    'a')
	return 't';

    if(c ==    'c')
	return 'g';

    if(c ==    'g')
	return 'c';

    if(c ==    't')
	return 'a';



    if(c ==    'N')
	return 'N';
    
    exit(1);
}
//Check if it is either A,C,G,T
inline bool isResolvedDNA(const char c){
    char _c= toupper(c);

    if(_c ==    'A')
	return true;
    if(_c ==    'C')
	return true;
    if(_c ==    'G')
	return true;
    if(_c ==    'T')
	return true;
    return false;
}

inline vector<mdField> mdString2Vector(const string & mdFieldToParse){
    vector<mdField> toReturn;
    int i=0;
    // int addToOffset=0;
    mdField toadd;

    toadd.offset=0;
    toadd.bp='N';

    while(int(mdFieldToParse.length()) != i){
	if(isdigit(mdFieldToParse[i])){
	    toadd.offset=toadd.offset*10+(int(mdFieldToParse[i])-asciiOffsetZero);
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
    return toReturn;
}

std::pair< std::string, std::vector<int> >  reconstructRefWithPosHTS(const bam1_t   * b){
  std::string mdFieldString="";
  std::string reconstructed="";
  std::string reconstructedTemp="";

    //skip unmapped
    if( ((b)->core.flag&BAM_FUNMAP) != 0 ){
      fprintf(stderr,"The function reconstructRefWithPosOnReadHTS()  cannot be called for unmapped reads\n");
      exit(1);
    }
    
    uint8_t *mdptr = bam_aux_get(b, "MD");
    
    if(mdptr){
      mdFieldString = std::string( (const char*)(mdptr+1));
    }else{
      fprintf(stderr,"ReconsReferenceHTSLIB: Cannot get MD tag from:%s ",bam_get_qname(b));
      exit(1);
    }

    int32_t   n_cigar_op = b->core.n_cigar;
    uint32_t *cigar      = bam_get_cigar(b);

    for(int32_t i = 0; i < n_cigar_op; i++){
	char opchr = bam_cigar_opchr(cigar[i]);
        int32_t oplen = bam_cigar_oplen(cigar[i]);
	reconstructedTemp+=std::string(oplen,opchr);
    }

    //get a vector representation of the MD field	

    std::vector<mdField> parsedMD=mdString2Vector(mdFieldString);

    vector<int> positionsOnControl;
    
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
			reconstructed+=parsedMD[mdVectorIndex].bp;
			positionsOnControl.push_back(initialPositionControl++);
		    }
		    mdVectorIndex++;
		}else{ //wait until we reach a mismatch
		    reconstructed+=reconstructedTemp[i];
		    parsedMD[mdVectorIndex].offset--;
		    positionsOnControl.push_back(initialPositionControl++);
		}

		//skipping all the positions with deletions on the read
		//if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

		while( (mdVectorIndex<int(parsedMD.size())) &&
		       (parsedMD[mdVectorIndex].bp == '^' ) ){ 
		    initialPositionControl+=parsedMD[mdVectorIndex].offset;
		    mdVectorIndex++;
		}
		    
	    }else{
		reconstructed+=reconstructedTemp[i];
		positionsOnControl.push_back(initialPositionControl++);
	    }
	}else{
	    if(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){ //soft clipped bases and indels
		reconstructed+=reconstructedTemp[i];
		positionsOnControl.push_back(initialPositionControl);
	    }
	}
    }

    if(int(reconstructed.size()) != b->core.l_qseq){
	cerr << "Could not recreate the sequence for read "<<bam_get_qname(b)  << endl;
	exit(1);
    }

    if(positionsOnControl.size() != reconstructed.size()){
	cerr << "Could not determine the positions for the read "<<bam_get_qname(b) << endl;
	exit(1);
    }


    return pair< string, vector<int> >(reconstructed,positionsOnControl);
}


using namespace std;
// using namespace BamTools;

const int offset=0;
int numberOfCycles;
string alphabetHTSLIB = "NACNGNNNTNNNNNNN";

#define MAXLENGTH 1000

vector< vector<unsigned int> > typesOfDimer5p; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p; //3' deam rates
vector< vector<unsigned int> > typesOfDimer5p_cpg; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p_cpg; //3' deam rates
vector< vector<unsigned int> > typesOfDimer5p_noncpg; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p_noncpg; //3' deam rates


vector< vector<unsigned int> > typesOfDimer5pDouble; //5' deam rates when the 3' is deaminated according to a double str.
vector< vector<unsigned int> > typesOfDimer3pDouble; //3' deam rates when the 5' is deaminated according to a double str.
vector< vector<unsigned int> > typesOfDimer5pSingle; //5' deam rates when the 3' is deaminated according to a single str.
vector< vector<unsigned int> > typesOfDimer3pSingle; //3' deam rates when the 5' is deaminated according to a single str.


//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(const   bam1_t  * b,string & reconstructedReference,const vector<int> &  reconstructedReferencePos,const int & minQualBase,const string & refFromFasta, const bam_hdr_t *h,bool ispaired,bool isfirstpair){ // ,int firstCycleRead,int increment

    char refeBase;
    char readBase;
    int  qualBase;
 
    //Checking if the 5' is deaminated
    bool isDeam5pS=false; //C->T 5'
    bool isDeam3pS=false; //C->T 3'
    bool isDeam5pD=false; //C->T 5'
    bool isDeam3pD=false; //G->A 3'

    int i;
    //cerr<<"read  "<<bam_get_qname(b)<<endl;
    if(ispaired){ //since we cannot evaluate the 5' ends or 3' ends
	goto iterateLoop;
    }


    i=0; //5p for forward str, 3p for reverse

    refeBase=toupper(reconstructedReference[i]);
    
    //readBase=toupper(         al.QueryBases[i]);
    readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);
    //qualBase=int(             al.Qualities[i])-offset;
    qualBase=int(             bam_get_qual(b)[i])-offset;  

    if(qualBase < minQualBase)
	goto eval3pdeam;

    if(refeBase == 'S' ||refeBase == 'I'){ //don't care about soft clipped or indels
	goto eval3pdeam;
    }
    
    if(refeBase == 'M'){//match
	refeBase =  readBase;
    }

    if( isResolvedDNA(refeBase)  && 
	isResolvedDNA(readBase) ){
	// if(al.IsReverseStrand()){ //need to take the complement
	if( bam_is_rev(b) ){
	    refeBase=complement(refeBase);
	    readBase=complement(readBase);
	}
	
	if(refeBase == 'C' &&
	   readBase == 'T' ){ //C->T

	    //if(al.IsReverseStrand()){ //3'
	    if( bam_is_rev(b) ){
		isDeam3pS=true;		
	    }else{                    //5'
		isDeam5pS=true;
		isDeam5pD=true;
	    }
	}


	if(refeBase == 'G' &&
	   readBase == 'A' ){ //G->A

	    //if(al.IsReverseStrand()){ //3'
	    if( bam_is_rev(b) ){
		isDeam3pD=true;		
	    }else{                    //5'
	    }
	}
	   
    }


 eval3pdeam:
    //i=int(al.QueryBases.size())-1; //3p for forward str, 5p for reverse
    i=b->core.l_qseq-1;
    
    refeBase=toupper(reconstructedReference[i]);
    // readBase=toupper(         al.QueryBases[i]);
    // qualBase=int(              al.Qualities[i])-offset;
    readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);
    qualBase=int(             bam_get_qual(b)[i])-offset;  
    
    if(qualBase < minQualBase)
	goto iterateLoop;
    
    if(refeBase == 'S' || refeBase == 'I'){ //don't care about soft clipped or indels
	goto iterateLoop;
    }
    
    if(refeBase == 'M'){//match
	refeBase =  readBase;
    }

    if( isResolvedDNA(refeBase)  && 
	isResolvedDNA(readBase) ){
	//if(al.IsReverseStrand()){ //need to take the complement
	if( bam_is_rev(b) ){
	    refeBase=complement(refeBase);
	    readBase=complement(readBase);
	}
	
	if(refeBase == 'C' &&
	   readBase == 'T' ){ //C->T

	    // if(al.IsReverseStrand()){ //5'
	    if( bam_is_rev(b) ){

		isDeam5pS=true;
		isDeam5pD=true;		
	    }else{                    //3'
		isDeam3pS=true;
	    }
	}

	if(refeBase == 'G' &&
	   readBase == 'A' ){ //G->A

	    // if(al.IsReverseStrand()){ //5'
	    if( bam_is_rev(b) ){
    
	    }else{                    //3'
		isDeam3pD=true;
	    }
	}
	   


    }

 iterateLoop:

    
    char refBaseFromFasta      = 'N';
    char refBaseFromFastaPrev  = 'N';
    char refBaseFromFastaNext  = 'N';
    int j=0;
    for(i=0;i<int(b->core.l_qseq);i++,j++){
	// cout<<i<<endl;

	refeBase=toupper(reconstructedReference[j]);

	// readBase=toupper(          al.QueryBases[i]);
	// qualBase=int(              al.Qualities[i])-offset;

	readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);
	qualBase=int(             bam_get_qual(b)[i])-offset;  
  
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

	    if(!refFromFasta.empty()){
		refBaseFromFasta         = refFromFasta[j+1];
		refBaseFromFastaPrev     = refFromFasta[j  ];
		refBaseFromFastaNext     = refFromFasta[j+2];		
		if(refeBase != refBaseFromFasta){
		    cerr<<"Discrepency#1 for "<<bam_get_qname(b)<<" where the reference base at position "<<i<<" "<<refeBase<<" "<<refBaseFromFasta<<endl;
		    exit(1);
		}

	    }
	    
	
	    
	}else{
	    if(!refFromFasta.empty()){
		refBaseFromFasta         = refFromFasta[j+1];
		refBaseFromFastaPrev     = refFromFasta[j  ];
		refBaseFromFastaNext     = refFromFasta[j+2];		
		if(refeBase != refBaseFromFasta){
		    cerr<<"Discrepency#2 for "<<bam_get_qname(b)<<" where the reference base at position "<<i<<" "<<refeBase<<" "<<refBaseFromFasta<<endl;
		    exit(1);
		}

	    }

	}

	// cout<<refBaseFromFastaPrev<<" "<<refBaseFromFasta<<" "<<refBaseFromFastaNext<<endl;
	
	if( isResolvedDNA(refeBase)  && 
	    isResolvedDNA(readBase) ){
	    int dist5p=i;
	    //int dist3p=int(al.QueryBases.size())-1-i;
	    int dist3p=b->core.l_qseq-1-i;
	    
	    //  if(al.IsReverseStrand()){ //need to take the complement
	    if( bam_is_rev(b) ){
		refeBase=complement(refeBase);
		readBase=complement(readBase);
		//dist5p=int(al.QueryBases.size())-1-i;
		dist5p=int(b->core.l_qseq)-1-i;
		dist3p=i;
	    }

	    if(dist5p > MAXLENGTH ||
	       dist3p > MAXLENGTH ){
		cerr<<"Molecule found "<<bam_get_qname(b)<<" with length greater than limit"<<endl;
		exit(1);
	    }
	       

	    //mismatches[cycleToUse]++;
	    if( !ispaired ||  isfirstpair){
		//cerr<<"increase 5p"<<endl;
		typesOfDimer5p[dist5p][twoBases2index(refeBase,readBase)]++;
	    }

	    if( !ispaired || !isfirstpair){
		//cerr<<"increase 3p"<<endl;
		typesOfDimer3p[dist3p][twoBases2index(refeBase,readBase)]++;
	    }
	    
	    if(!refFromFasta.empty()){
		if(
		    ( (refBaseFromFasta     == 'C' && refBaseFromFastaNext == 'G') && !bam_is_rev(b) ) //!al.IsReverseStrand() )
		    ||
		    ( (refBaseFromFastaPrev == 'C' && refBaseFromFasta     == 'G') &&  bam_is_rev(b) ) //al.IsReverseStrand() )
		){
		    //cout<<"   CPG: "<<refBaseFromFastaPrev<<" "<<refBaseFromFasta<<" "<<refBaseFromFastaNext<<" ref:"<<refeBase<<" read:"<<readBase<<" "<<al.IsReverseStrand()<<" same="<<(refeBase==readBase)<<endl;
		    if( !ispaired ||  isfirstpair)
			typesOfDimer5p_cpg[dist5p][twoBases2index(refeBase,readBase)]++;
		    if( !ispaired || !isfirstpair)
			typesOfDimer3p_cpg[dist3p][twoBases2index(refeBase,readBase)]++;

		}else{
		    if( isResolvedDNA(refBaseFromFasta)                               &&
			isResolvedDNA(refBaseFromFastaPrev)                           &&
			isResolvedDNA(refBaseFromFastaNext)                           &&
			!(refBaseFromFasta     == 'C' && refBaseFromFastaNext == 'G') &&
			!(refBaseFromFastaPrev == 'C' && refBaseFromFasta     == 'G')
		    ){

		       //cout<<"nonCPG: "<<refBaseFromFastaPrev<<" "<<refBaseFromFasta<<" "<<refBaseFromFastaNext<<" ref:"<<refeBase<<" read:"<<readBase<<" "<<al.IsReverseStrand()<<" same="<<(refeBase==readBase)<<endl;
			if( !ispaired ||  isfirstpair)
			    typesOfDimer5p_noncpg[dist5p][twoBases2index(refeBase,readBase)]++;
			if( !ispaired || !isfirstpair)
			    typesOfDimer3p_noncpg[dist3p][twoBases2index(refeBase,readBase)]++;
			
		    }
		}
	    }

	    if(isDeam5pS){
		if( !ispaired || !isfirstpair)
		    typesOfDimer3pSingle[dist3p][twoBases2index(refeBase,readBase)]++;
	    }

	    if(isDeam3pS){
		if( !ispaired ||  isfirstpair)
		    typesOfDimer5pSingle[dist5p][twoBases2index(refeBase,readBase)]++;
	    }


	    if(isDeam5pD){
		if( !ispaired || !isfirstpair)
		    typesOfDimer3pDouble[dist3p][twoBases2index(refeBase,readBase)]++;
	    }

	    if(isDeam3pD){
		if( !ispaired ||  isfirstpair)
		    typesOfDimer5pDouble[dist5p][twoBases2index(refeBase,readBase)]++;
	    }


	}
    }
}


double dbl2log(const double d,bool phred){
    double t= -10.0*(log(d)/log(10.0));
    // if(d == 0){
    // 	t = 
    // }
    if(phred)
	return t;
    else 
	return d;
}


int main (int argc, char *argv[]) {

    string file5p="/dev/stdout";
    string file3p="/dev/stdout";

    bool allStr   =true;
    bool singleStr=false;
    bool doubleStr=false;
    bool singAnddoubleStr=false;

    int lengthMaxToPrint = 5;
    int minQualBase      = 0;
    int minLength        = 35;
    bool dpFormat=false;
    bool hFormat=false;
    double errorToRemove=0.0;
    bool phred=false;
    bool cpg=false;
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
			"\t\t"+"-err\t[error rate]\tSubstract [error rate] from the rates to account for sequencing errors  (Default: "+stringify(errorToRemove)+" ) \n"+
			"\t\t"+"-log\t\t\tPrint substitutions on a PHRED logarithmic scale  (Default: "+stringify(phred)+" ) \n"+
			"\t\t"+"-paired\t\t\tAllow paired reads    (Default: "+booleanAsString( paired )+" ) \n"+


			"\n\n\tYou can specify either one of the two:\n"+
			"\t\t"+"-single\t\t\tUse the deamination profile of a single strand library  (Default: "+booleanAsString( singleStr )+")\n"+
			"\t\t"+"-double\t\t\tUse the deamination profile of a double strand library  (Default: "+booleanAsString( doubleStr )+")\n"+
			"\n\tor specify this option:\n"+
			"\t\t"+"-both\t\t\tReport both C->T and G->A regardless of stand  (Default: "+booleanAsString( singAnddoubleStr )+")\n"+

			"\n\n\tOutput options:\n"+
			"\t\t"+"-5p\t[output file]\tOutput profile for the 5' end (Default: "+stringify(file5p)+")\n"+
			"\t\t"+"-3p\t[output file]\tOutput profile for the 3' end (Default: "+stringify(file3p)+")\n"+
			"\t\t"+"-dp\t\t\tOutput in damage-patterns format (Default: "+booleanAsString(dpFormat)+")\n"+
			"\t\t"+"-h\t\t\tMore human readible output (Default: "+booleanAsString(hFormat)+")\n"+
			"\t\t"+"-q\t\t\tDo not print why reads are skipped (Default: "+booleanAsString(quiet)+")\n"+
		       
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "--help") )
    ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    
    for(int i=1;i<(argc-1);i++){ //all but the last 3 args


        if(string(argv[i]) == "-dp"  ){
            dpFormat=true;
            continue;
        }

        if(string(argv[i]) == "-log"  ){
            phred=true;
            continue;
        }


        if(string(argv[i]) == "-paired"  ){
            paired=true;
            continue;
        }

        if(string(argv[i]) == "-h"  ){
            hFormat=true;
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

	if(string(argv[i]) == "-cpg"  ){
	    cpg=true;
            continue;
        }

        if(string(argv[i]) == "-length"  ){
            lengthMaxToPrint=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-err"  ){
            errorToRemove=destringify<double>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-5p" ){
	    file5p = string(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-3p" ){
	    file3p = string(argv[i+1]);
	    i++;
            continue;
        }

	if(string(argv[i]) == "-both" ){
	    //doubleStr=true;

	    allStr           = false;
	    singleStr        = false;
	    doubleStr        = false;
	    singAnddoubleStr = true;
            continue;
        }


        if(string(argv[i]) == "-single" ){

	    allStr    = false;
	    singleStr = true;
	    doubleStr = false;

            continue;
        }

        if(string(argv[i]) == "-double" ){
	    //doubleStr=true;

	    allStr    = false;
	    singleStr = false;
	    doubleStr = true;

            continue;
        }


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    if(phred && hFormat){
	cerr<<"Error: cannot specify both -log and -h"<<endl;
	return 1;
    }

    if(dpFormat && hFormat){
	cerr<<"Error: cannot specify both -dp and -h"<<endl;
	return 1;
    }

    typesOfDimer5p       = vector< vector<unsigned int> >();
    typesOfDimer3p       = vector< vector<unsigned int> >();
    typesOfDimer5p_cpg   = vector< vector<unsigned int> >();
    typesOfDimer3p_cpg   = vector< vector<unsigned int> >();
    typesOfDimer5p_noncpg= vector< vector<unsigned int> >();
    typesOfDimer3p_noncpg= vector< vector<unsigned int> >();
    
    typesOfDimer5pDouble = vector< vector<unsigned int> >();
    typesOfDimer3pDouble = vector< vector<unsigned int> >();
    typesOfDimer5pSingle = vector< vector<unsigned int> >();
    typesOfDimer3pSingle = vector< vector<unsigned int> >();

    for(int l=0;l<MAXLENGTH;l++){
	//for(int i=0;i<16;i++){
	typesOfDimer5p.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5p_cpg.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p_cpg.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5p_noncpg.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p_noncpg.push_back( vector<unsigned int> ( 16,0 ) );


	typesOfDimer5pDouble.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3pDouble.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5pSingle.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3pSingle.push_back( vector<unsigned int> ( 16,0 ) );

	//}
    }
    
    string bamfiletopen = string( argv[ argc-1 ] );
    // string deambam      = string( argv[ argc-2 ] );
    // string nondeambam   = string( argv[ argc-1 ] );

    samFile  *fp;
    bam1_t    *b;
    bam_hdr_t *h;

    string refFromFasta_;
    string refFromFasta;
    
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
	
	pair< string, vector<int> >  reconstructedReference = reconstructRefWithPosHTS(b);
	increaseCounters(b,reconstructedReference.first, reconstructedReference.second,minQualBase,refFromFasta,h,ispaired,isfirstpair); //start cycle numberOfCycles-1
    }
    
    bam_destroy1(b);
    sam_close(fp);
    
    
  
    ofstream file5pFP;
    file5pFP.open(file5p.c_str());

    if (!file5pFP.is_open()){
	cerr << "Unable to write to 5p file "<<file5p<<endl;
	exit(1);
    }

    
    if(dpFormat)
	file5pFP<<"\t";
    if(hFormat)
	file5pFP<<"pos\t";
    file5pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;
  

    vector< vector<unsigned int> > * typesOfDimer5pToUse;

    typesOfDimer5pToUse     = &typesOfDimer5p;
    
    for(int l=0;l<lengthMaxToPrint;l++){
	if(dpFormat)
	    file5pFP<<l<<"\t";

	if(hFormat)
	    file5pFP<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";
	
	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer5pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		if(allStr){
		    if(dpFormat)
		      file5pFP<<dbl2log(std::max(0.0,double( (*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
		    else
			if(hFormat)
			  file5pFP<<printDoubleAsWhitePaddedString( std::max(0.0,double( (*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
			else
			  file5pFP<<dbl2log( std::max(0.0,double( (*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred);  
		}else{ 
		    if(singAnddoubleStr){
			if(         (n1==1 && n2==3) || (n1==2 && n2==0 )  ) { 
			    if(dpFormat)
			      file5pFP<<dbl2log( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
			    else
				if(hFormat)
				  file5pFP<<printDoubleAsWhitePaddedString( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),1,5);
				else
				  file5pFP<<dbl2log( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			} else { 
			    if(dpFormat)
				file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
			    else
				if(hFormat)
				    file5pFP<<printDoubleAsWhitePaddedString( 0.0,1,5);
				else
				    file5pFP<<(phred?"-Inf":"0.0")<<"";
			}

		    }else{
			if(doubleStr){
			    //          C        T
			    if(         n1==1 && n2==3  ) { 
				if(dpFormat)
				  file5pFP<<dbl2log( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
				else
				    if(hFormat)
				      file5pFP<<printDoubleAsWhitePaddedString( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5); 
				    else
					file5pFP<<dbl2log( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			    } else { 
				if(dpFormat)
				    file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
				else
				    if(hFormat)
					file5pFP<<printDoubleAsWhitePaddedString( 0.0,1,5);
				    else					
					file5pFP<<(phred?"-Inf":"0.0"); 
			    }
			}else{ 
			    if(singleStr){
				//      C        T
				if(     n1==1 && n2==3  ) { 
				    if(dpFormat)
					file5pFP<<dbl2log( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
				    else
					if(hFormat)
					    file5pFP<<printDoubleAsWhitePaddedString(std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),1,5);
					else
					    file5pFP<<dbl2log( std::max(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
					    
				} else { 
				    if(dpFormat)
					file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
				    else
					if(hFormat)
					    file5pFP<<printDoubleAsWhitePaddedString(0.0,1,5);
					else										    
					    file5pFP<<(phred?"-Inf":"0.0"); 
				}
			    }
			}
		    }
		}

		
		if(!(n1 ==3 && n2 == 2 ))
		    file5pFP<<"\t";
	    }


	}
	file5pFP<<endl;
    }


    file5pFP.close();

    ofstream file3pFP;
    if(file3p == "/dev/stdout"){
	file3pFP.open(file3p.c_str(), ofstream::out | ofstream::app);
    }else{
	file3pFP.open(file3p.c_str());
    }

    if (!file3pFP.is_open()){
	cerr << "Unable to write to 3p file "<<file3p<<endl;
	exit(1);
    }

    if(dpFormat)
	file3pFP<<"\t";
    if(hFormat)
	file3pFP<<"pos\t";

    file3pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;


    vector< vector<unsigned int> > * typesOfDimer3pToUse;

    typesOfDimer3pToUse     = &typesOfDimer3p;
    
    for(int le=0;le<lengthMaxToPrint;le++){

	int l=le;
	if(dpFormat){
	    l=lengthMaxToPrint-1-le;
	    if(l==0)
		file3pFP<<""<<l<<"\t";	    
	    else
		file3pFP<<"-"<<l<<"\t";	    
	}

	if(hFormat){
	    //l=lengthMaxToPrint-1-le;
	    // if(l==0)
	    // 	file3pFP<<""<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint)))<<"\t";	    
	    // else
	    // 	file3pFP<<"-"<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint)))<<"\t";	    
	    file3pFP<<""<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";	    
	}

	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer3pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		if(allStr){
		    if(dpFormat)
			file3pFP<<dbl2log( std::max(0.0,double( (*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
		    else
			if(hFormat)
			    file3pFP<<printDoubleAsWhitePaddedString( std::max(0.0,double( (*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
			else
			    file3pFP<<dbl2log( std::max(0.0,double( (*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred);
		}else{ 
		    if(singAnddoubleStr){			
			if(   (n1==1 && n2==3) || (n1==2 && n2==0 )  ) { 
			    if(dpFormat)
				file3pFP<<dbl2log( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
			    else
				if(hFormat)
				    file3pFP<<printDoubleAsWhitePaddedString( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
				else
				    file3pFP<<dbl2log( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			} else { 
			    if(dpFormat)
				file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
			    else
				if(hFormat)
				    file3pFP<<printDoubleAsWhitePaddedString( 0.0 ,1,5);
				else
				    file3pFP<<(phred?"-Inf":"0.0"); 
			}

		    }else{
			if(doubleStr){
			    //          G        A
			    if(         n1==2 && n2==0  ) { 
				if(dpFormat)
				    file3pFP<<dbl2log( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
				else
				    if(hFormat)
					file3pFP<<printDoubleAsWhitePaddedString( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
				    else					
					file3pFP<<dbl2log( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			    } else { 				
				if(dpFormat)
				    file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
				else
				    if(hFormat)
					file3pFP<<printDoubleAsWhitePaddedString(  0.0 ,1,5);
				    else	
					file3pFP<<(phred?"-Inf":"0.0"); 
			    }
			}else{ 
			    if(singleStr){
				//      C        T
				if(     n1==1 && n2==3  ) { 
				    if(dpFormat)
					file3pFP<<dbl2log( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
				    else
					if(hFormat)
					    file3pFP<<printDoubleAsWhitePaddedString( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
					else					
					    file3pFP<<dbl2log( std::max(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
				} else { 				    
				    if(dpFormat)
					file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
				    else
					if(hFormat)
					    file3pFP<<printDoubleAsWhitePaddedString(  0.0 , 1,5);
					else	
					    file3pFP<<(phred?"-Inf":"0.0"); 

				}
				
			    }
			}
		    }
		}

		
		if(!(n1 ==3 && n2 == 2 ))
		    file3pFP<<"\t";
	    }


	}
	file3pFP<<endl;
    }

    file3pFP.close();
    return 0;
}
