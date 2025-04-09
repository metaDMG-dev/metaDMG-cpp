#include <cassert>
#include <cstdio>
#include <htslib/bgzf.h>
#include <climits>
#include "profile.h"
#include <string>     
#include <vector>  
#include <array>     
#include <map>    
#include <cstdlib>
#include <htslib/kstring.h>


void merge_bdamage(const std::vector<std::string> &bdamage_files, const char* outname) {


  std::vector<std::map<int, mydataD> > myvec;
  int globhowmany = -1;
  for (const std::string &fname : bdamage_files) {

    int howmany;
    std::map<int, mydataD> retmap = load_bdamage_full(fname.c_str(), howmany);
    auto stupid = retmap.begin();

    fprintf(stderr, "\t-> Number of entries in damage pattern file: %lu printlength(howmany):%d howmany2: %d\n", retmap.size(), howmany,stupid->second.howmany);
    assert(howmany==stupid->second.howmany);
    if(globhowmany==-1)
      globhowmany = stupid->second.howmany;
    if(globhowmany!=stupid->second.howmany){
      fprintf(stderr,"\t-> We need the same number of mismatch cycles for merging bdamage files \n");
    }
    myvec.push_back(std::move(retmap));
  }

  // now merge all but the first one, into the first one. Either adding or making new ids
  fprintf(stderr,"\t-> Doing merging old size: %lu\n",myvec[0].size());
  for(int i=1;i<myvec.size();i++){
    fprintf(stderr,"\t-> Merging %d into 0\n",i);
    std::map<int, mydataD> &small= myvec[i];
    
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
}

void merge_rlens(const std::vector<std::string> &rlens_files, const char* outname) {
  
  // results 
  std::map<int, std::array<size_t, 200>> rlens_merged;

  // each file 
  for (const std::string &fname : rlens_files) {

    BGZF *fp = bgzf_open(fname.c_str(), "r");
    if (!fp) {
      fprintf(stderr, "could not open rlens file: %s\n", fname.c_str());
      exit(1);
    }

    kstring_t line = {0, 0, NULL};
    int ret = bgzf_getline(fp, '\n', &line); 
    if (ret < 0) {
      fprintf(stderr, "empty or invalid rlens file: %s\n", fname.c_str());
      exit(1);
    }

    while ((ret = bgzf_getline(fp, '\n', &line)) >= 0) {
      int taxid;
      size_t tmp[200];
      char *ptr = line.s;
      char *endptr;

      // taxid is first column
      taxid = strtol(ptr, &endptr, 10); 
      ptr = endptr;

      // each new cell (200)
      for (int j = 0; j < 200; j++) {
        tmp[j] = strtoull(ptr, &endptr, 10);
        ptr = endptr;
      }

      // sum and merge
      auto &entry = rlens_merged[taxid];
      for (int j = 0; j < 200; j++) {
        entry[j] += tmp[j];
      }
    }

    bgzf_close(fp);
    free(line.s);
  }


  // stolen from bwrite
  char onam[1024];
  snprintf(onam, sizeof(onam), "%s.rlens.gz", outname);
  BGZF *fp = my_bgzf_open(onam, 2);
  assert(fp);

  kstring_t kstr = {0, 0, NULL};
  ksprintf(&kstr, "id");
  for (int i = 0; i < 200; i++)
    ksprintf(&kstr, "\trlen%d", i);
  ksprintf(&kstr, "\n");

  for (const auto &it : rlens_merged) {
    ksprintf(&kstr, "%d", it.first);
    for (int i = 0; i < 199; i++)
      ksprintf(&kstr, "\t%lu", it.second[i]);
    ksprintf(&kstr, "\t%lu\n", it.second[199]);

    if (kstr.l > 1000000) {
      assert(bgzf_write(fp, kstr.s, kstr.l) == (ssize_t)kstr.l);
      kstr.l = 0;
    }
  }

  if (kstr.l > 0) {
    assert(bgzf_write(fp, kstr.s, kstr.l) == (ssize_t)kstr.l);
  }

  bgzf_close(fp);
  free(kstr.s);

  fprintf(stderr, "\t-> Done writing merged rlens: %s\n", onam);
}


int main_mergedamage(int argc, char **argv) {
  fprintf(stderr, "[mergedamage] Starting mergedamage\n");

  std::vector<std::string> bdamage_files;
  std::vector<std::string> rlens_files;
  const char* outname = nullptr;

  argc--; argv++;

  enum ParseMode { NONE, BFILES, RFILES } mode = NONE;

  // get args 
  for (int i = 0; i < argc; i++) {
    if (strcasecmp(argv[i], "-b") == 0) {
      mode = BFILES;
      continue;
    }
    if (strcasecmp(argv[i], "-r") == 0) {
      mode = RFILES;
      continue;
    }
    if (strcasecmp(argv[i], "-out") == 0 || strcasecmp(argv[i], "-outnames") == 0) {
    if (i + 1 >= argc) {
      fprintf(stderr, "Error: -out flag requires a value\n");
      exit(1);  
    }
    outname = argv[++i];
    continue;
    }

  if (mode == BFILES) {
      bdamage_files.push_back(argv[i]);
    } else if (mode == RFILES) {
      rlens_files.push_back(argv[i]);
    } 
  }

  // checks 
  if (bdamage_files.empty() && rlens_files.empty()) {
    fprintf(stderr, "no input given\n");
    return 1;
  }

  if (!outname) {
  fprintf(stderr, "no output given\n");
  return 1;
  }

  if (bdamage_files.size() == 1) {
    fprintf(stderr, "\t> merging one bdamage file is not meaningful\n");
  }

  if (rlens_files.size() == 1) {
    fprintf(stderr, "\t> merging one rlens file is not meaningful\n");
  }

  // do the thing
  if (!bdamage_files.empty()) {
    fprintf(stderr, "[mergedamage] merging %lu bdamage files\n", bdamage_files.size());
    merge_bdamage(bdamage_files, outname);
  }

  if (!rlens_files.empty()) {
    fprintf(stderr, "[mergedamage] merging %lu rlens.gz files\n", rlens_files.size());
    merge_rlens(rlens_files, outname);
  }

  return 0;
}

