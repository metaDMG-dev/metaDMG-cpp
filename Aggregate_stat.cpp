/*
  16 july, this file contains functionality from the print_ugly function and codeparts in the misc.
  The function will implement the ML estimate of the dfit from metaDMG-cpp bioxarhive paper.
  
 */

#include <random>
#include <cstdio>
#include <htslib/hts.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/sam.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/bgzf.h>
#include <ctime>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <algorithm> // for std::sort
#include <cassert>

#include <iostream>
#include "profile.h"
#include "shared.h"
#include "ngsLCA.h" //<- print_chain
#include "types.h"       // for int2intvec, int2int
#include "version.h"     // for METADAMAGE_VERSION
#include "dfit_optim.h"
#include "pval.h"
#include "dfit_helppage.h"

extern htsFormat *dingding2;

mydataD getval_full(std::map<int, mydataD> &retmap, int2intvec &child, int taxid, int howmany);
mydata2 getval_stats(std::map<int, mydata2> &retmap, int2intvec &child, int taxid) ;

int HelpPageAggregate(FILE *fp){
  fprintf(fp,"Aggregation of lca produced statistics (mean length, variance length, mean GC, variance GC) when transversing up the nodes of the tree structure\n");

  fprintf(stderr, "\t\t./metaDMG-cpp aggregate file.bdamage.gz --names file.gz --nodes trestructure.gz --lcastat file.stat --out filename\n");
  fprintf(fp,"\n--help \t\t\t\t Print extended help page to see all options.\n\n");
  fprintf(fp,"\n--names \t\t\t\t names.dmp.gz\n\n");
  fprintf(fp,"\n--nodes \t\t\t\t nodes.dmp.gz\n\n");
  fprintf(fp,"\n--lca \t\t\t\t lcaout.stat lca produced statistics\n\n");

  exit(1);
  return 0;
}

void to_root(int from,int to,std::map<int,mydata2> &stats,int2int &parent,int nreads){
  //  fprintf(stderr,"from: %d to: %d nreads:%d\n",from,to,nreads);
  mydata2 &md1 = stats.find(from)->second;
  
  std::map<int,mydata2>::iterator it=stats.find(to);
  if(it==stats.end()){
    mydata2 mdmis;
    mdmis.nreads = md1.nreads;

    mdmis.data = new double[4];
    for(int iii=0;iii<4;iii++)
      mdmis.data[iii] = md1.data[iii];

    stats[to] = mdmis;
  }
  else{
    mydata2 &md2 = stats.find(to)->second;
    md2.data[0] = ((double) md1.data[0]*md1.nreads+md2.data[0]*md2.nreads)/((double) md1.nreads+md2.nreads); //weighted mean of read length
    md2.data[2] = ((double) md1.data[2]*md1.nreads+md2.data[2]*md2.nreads)/((double) md1.nreads+md2.nreads); //weighted mean of gc
    
    double variance1;
    int nreads1;
    double variance2;
    int nreads2;
    if(((double) md1.nreads+md2.nreads)>2){//pooled variance of length and GC
      /*variance1 = ((double) md1.nreads-1)*md1.data[1];
      nreads1 += md1.nreads;
      variance2 = ((double) md2.nreads-1)*md2.data[1];
      nreads2 += md2.nreads;
      md2.data[1] = (variance1 + variance2) / (nreads1 + nreads2 - 2);
      md2.data[3] = (variance1 + variance2) / (nreads1 + nreads2 - 2);*/
      md2.data[1] = (((double) md1.nreads-1)*md1.data[1]+((double) md2.nreads-1)*md2.data[1])/((double)md1.nreads+md2.nreads-2);
      md2.data[3] = (((double) md1.nreads-1)*md1.data[3]+((double) md2.nreads-1)*md2.data[3])/((double)md1.nreads+md2.nreads-2);;
    }
    else{
      md2.data[1] = md2.data[1];
      md2.data[3] = md2.data[3];
    }

    md2.nreads += nreads;
  }

  int newto = parent.find(to)->second;
  
  if(newto!=to)
    to_root(to,newto,stats,parent,nreads);
  
}

void aggr_stat2000(std::map<int, mydata2> &stats,int2int &parent){
  std::map<int,int> dasmap;
  for(std::map<int,mydata2>::iterator it = stats.begin();it!=stats.end();it++)
    dasmap[it->first] = it->second.nreads;
  //fprintf(stderr,"dasmap.size(): %lu stats.size():%lu\n",dasmap.size(), stats.size());

  for(std::map<int,int>::iterator itt=dasmap.begin();itt!=dasmap.end();itt++){
  //  for(int i=0;i<dasvector.size();i++){
    int focal_taxid = itt->first;
    int2int::iterator it = parent.find(focal_taxid);
    assert(it!=parent.end());
    int target = it->second;
    to_root(focal_taxid,target,stats,parent,itt->second);
  }
  
}

int main_aggregate(int argc, char **argv) {
    if (argc <= 1){
      fprintf(stderr,"help\n");
      return 0;
    }
    char *infile_bdamage = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;
    char *infile_lcastat = NULL;
    char *outfile_name = NULL;
    int howmany;//this is the cycle

    while (*(++argv)) {
        if (strcasecmp("-h", *argv) == 0)
          fprintf(stderr,"help\n");
        else if (strcasecmp("--names", *argv) == 0 || strcasecmp("-names", *argv) == 0)
            infile_names = strdup(*(++argv));
        else if (strcasecmp("--nodes", *argv) == 0 || strcasecmp("-nodes", *argv) == 0)
            infile_nodes = strdup(*(++argv));
        else if (strcasecmp("-lca", *argv) == 0|| strcasecmp("--lcastat", *argv) == 0|| strcasecmp("-lcastat", *argv) == 0)
            infile_lcastat = strdup(*(++argv));
        else if (strcasecmp("-o", *argv) == 0 || strcasecmp("--out", *argv) == 0 || strcasecmp("--out_prefix", *argv) == 0)
            outfile_name = strdup(*(++argv));
        else
          infile_bdamage = strdup(*argv);
    }
    if(infile_nodes&&!infile_names){
      fprintf(stderr,"\t-> --names file.txt.gz must be defined with --nodes is defined\n");
      exit(1);
    }
    
    if(outfile_name==NULL)
      outfile_name = strdup(infile_bdamage);
    char buf[1024];
    snprintf(buf, 1024, "%s.aggregate.stat.txt.gz", outfile_name);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    BGZF *fpfpfp = bgzf_open(buf, "wb");
    kstring_t *kstr = new kstring_t;
    kstr->s = NULL; kstr->l = kstr->m = 0;
    ksprintf(kstr, "#taxid\tname\trank\tnalign\tnreads\tmean_rlen\tvar_rlen\tmean_gc\tvar_gc\tlca\ttaxa_path\n");

    // map of taxid -> taxid
    int2int parent;
    // map of taxid -> rank
    int2char rank;
    // map of parent -> child taxids
    int2intvec child;

    if (infile_nodes != NULL)
        parse_nodes(infile_nodes, rank, parent, child, 1);

    std::map<int, mydataD> retmap = load_bdamage_full(infile_bdamage, howmany);
    fprintf(stderr, "\t-> Number of entries in damage pattern file: %lu printlength(howmany):%d\n", retmap.size(), howmany);
    
    int2char name_map;

    if (infile_names)
        name_map = parse_names(infile_names);

    float presize = retmap.size();
    if(child.size()>0)
      getval_full(retmap, child, 1, howmany);  // this will do everything
    float postsize = retmap.size();
    fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);

    for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        mydataD md = it->second;
        if (it->second.nreads == 0)
	  continue;
	//	fprintf(stderr,"retmap taxid:%d nreads: %d\n",it->first,it->second.nreads);
    }

    std::map<int, mydata2> stats;
    if (infile_lcastat){
      stats = load_lcastat(infile_lcastat,1);
      presize = stats.size();
    }
    
#if 0
    if(child.size()>0)//this will not work if we have data at internal nodes 
      getval_stats(stats, child, 1);  // this will do everything
#endif
    aggr_stat2000(stats,parent);
    //    exit(0);
    if(stats.size()>0){
      postsize = stats.size();
      fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);
    }
    for (std::map<int, mydata2>::iterator it = stats.begin(); it != stats.end(); it++) {
      std::map<int, mydataD>::iterator itold = retmap.find(it->first);
        int nalign = -1;
        if (itold == retmap.end()) {
            fprintf(stderr, "\t-> Problem finding taxid: %d\n", it->first);
        } else
            nalign = itold->second.nreads;
        char *myrank = NULL;
        char *myname = NULL;
        if (it->second.nreads > 0) {
            int2char::iterator itc = rank.find(it->first);
            if (itc != rank.end())
                myrank = itc->second;
            itc = name_map.find(it->first);
            if (itc != name_map.end())
                myname = itc->second;
            ksprintf(kstr, "%d\t\"%s\"\t\"%s\"\t%d\t%d\t%f\t%f\t%f\t%f", it->first, myname, myrank, nalign, it->second.nreads, it->second.data[0], it->second.data[1], it->second.data[2], it->second.data[3]);
	    if(child.size()>0)
	      print_chain(kstr, it->first, parent, rank, name_map);
	    else
	      ksprintf(kstr,"NA\tNA\n");
            //      fprintf(stderr,"%d->(%d,%f,%f,%f,%f)\n",it->first,it->second.nreads,it->second.data[0],it->second.data[1],it->second.data[2],it->second.data[3]);
        }
    }
    //cleanup
    
    if(fpfpfp){
      if(bgzf_write(fpfpfp,kstr->s,kstr->l) == 0){
        fprintf(stderr, "\t-> Cannot write to output BGZ file\n");
        exit(1);
      }
      bgzf_close(fpfpfp);
    }

    for(int2char::iterator it=name_map.begin();it!=name_map.end();it++)
      free(it->second);
    for(int2char::iterator it=rank.begin();it!=rank.end();it++)
      free(it->second);

    for( std::map<int, mydataD>::iterator it = retmap.begin();it!=retmap.end();it++){
      mydataD md = it->second;
      delete [] md.fwD;
      delete [] md.bwD;
    }

    for( std::map<int, mydata2>::iterator it = stats.begin();it!=stats.end();it++){
      mydata2 md = it->second;
      delete [] md.data;
    }
    free(kstr->s);
    delete kstr;
    
    if(infile_bdamage)
      free(infile_bdamage);
    if(infile_nodes)
      free(infile_nodes);
    if(infile_names)
      free(infile_names);
    if(infile_lcastat)
      free(infile_lcastat);

  return 0;
}