#include <cassert>
#include <cstdio>
#include <htslib/bgzf.h>
#include <climits>
#include "profile.h"


int main_mergedamage(int argc, char **argv){
  fprintf(stderr,"information");
  fprintf(stderr,"argc: %d\n",argc);
  argc--;argv++;
  std::vector<std::map<int, mydataD> > myvec;
  char *outname = NULL;
  int globhowmany = -1;
  for(int i=0;i<argc;i++){
    // fprintf(stderr,"argc: %d val: %s\n",i,argv[i]);
    if(strcasecmp("-out",argv[i])==0||strcasecmp("-outnames",argv[i])==0){
      fprintf(stderr,"\t-> Masterhit match: argc: %d i: %d\n",argc,i);
      outname = strdup(argv[i+1]);
      i++;
      continue;
    }
    int howmany;
    char *fname = argv[i];
    std::map<int, mydataD> retmap = load_bdamage_full(fname, howmany);
    auto stupid = retmap.begin();
    fprintf(stderr, "\t-> Number of entries in damage pattern file: %lu printlength(howmany):%d howmany2: %d\n", retmap.size(), howmany,stupid->second.howmany);
    assert(howmany==stupid->second.howmany);
    if(globhowmany==-1)
      globhowmany = stupid->second.howmany;
    if(globhowmany!=stupid->second.howmany){
      fprintf(stderr,"\t-> We need the same number of mismatch cycles for merging bdamage files \n");
    }
    myvec.push_back(retmap);
  }
  fprintf(stderr,"\t-> outnames: %s size of bdamages: %lu\n",outname,myvec.size());
  if(outname==NULL){
    fprintf(stderr,"\t-> Please supply -out or -outnames as output prefix\n");
    exit(1);
  }
  if(myvec.size()<2){
    fprintf(stderr,"\t> Not really meaningfull to only merge one bdamage\n");
  }

  // now merge all but the first one, into the first one. Either adding or making new ids
  fprintf(stderr,"\t-> Doing merging old size: %lu\n",myvec[0].size());
  for(int i=1;i<myvec.size();i++){
    fprintf(stderr,"\t-> Merging %d into 0\n",i);
    std::map<int, mydataD> &small= myvec[1];
    
    for(auto slave =small.begin();slave!=small.end();slave++){
      std::map<int,mydataD>::iterator master = myvec[0].find(slave->first);//master will contain the results
      if(master!=myvec[0].end()){
	fprintf(stderr,"\t-> ID exists: %d\n",slave->first);// so id 
	//SO taxid in myvec[i] exists in myvec[0]. Let us just update the values;
	assert(slave->first==master->first);
	for(int cyc=0;cyc<globhowmany;cyc++){
	  master->second.fwD[cyc] += slave->second.fwD[cyc];
	  master->second.bwD[cyc] += slave->second.bwD[cyc];
	}
	master->second.nal += slave->second.nal;//update the nreads
      }else{
	fprintf(stderr,"\t-> ID does not exist will do full merge of entry: %d\n",slave->first);//do we need deep copy or is this ok? We havent dealloced so it should be fine...
	myvec[0][slave->first] = slave->second;
      }
    }
  }
  fprintf(stderr,"\t-> Done merging new size: %lu\n",myvec[0].size());
  char onam[1024];
  snprintf(onam, 1024, "%s.bdamage.gz", outname);
  fprintf(stderr, "\t-> Will dump: \'%s\' this contains damage patterns for: %lu items\n", onam, myvec.size());
  BGZF *fp = my_bgzf_open(onam, 2);//two threads, doesnt matter
  assert(bgzf_write(fp, &globhowmany, sizeof(int)) == sizeof(int));
  for(std::map<int,mydataD>::iterator it=myvec[0].begin();it!=myvec[0].end();it++){
    int taxid_nreads[2];
    taxid_nreads[0] = it->first;
    if(it->second.nal>INT_MAX){
      fprintf(stderr,"\t-> Problem merging entries sum of alignments are larger than maximum integer\n");
      assert(1!=0);
    }
    taxid_nreads[1] = it->second.nal;
    assert(bgzf_write(fp, taxid_nreads, 2*sizeof(int)) == 2*sizeof(int));
    for (int i = 0; i < globhowmany; i++) {
      float tmp[16];
      for (int ii = 0; ii < 16; ii++)
	tmp[ii] = it->second.fwD[i * 16 + ii];
      assert(16 * sizeof(float) == bgzf_write(fp, tmp, sizeof(float) * 16));
    }
    for (int i = 0; i < globhowmany; i++) {
      float tmp[16];
      for (int ii = 0; ii < 16; ii++)
	tmp[ii] = it->second.bwD[i * 16 + ii];
      assert(16 * sizeof(float) == bgzf_write(fp, tmp, sizeof(float) * 16));
    }
  }
  bgzf_close(fp);
  return 0;
}
