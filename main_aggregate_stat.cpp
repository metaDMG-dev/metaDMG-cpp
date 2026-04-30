/*
  16 july, this file contains functionality from the print_ugly function and codeparts in the misc.
  The function will implement the ML estimate of the dfit from metaDMG-cpp bioxarhive paper.
  //tsk,rh 12march 2024
 */

#include <random>
#include <cstdio>
#include <htslib/hts.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/sam.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/bgzf.h>
#include <ctime>
#include <sys/time.h>
#include <math.h>
#include <algorithm> // for std::sort
#include <cstdint>
#include <vector>

#include <iostream>
#include "profile.h"
#include "shared.h"
#include "ngsLCA.h" //<- print_chain
#include "types.h"       // for int2intvec, int2int
#include "version.h"     // for METADAMAGE_VERSION
#include "dfit_optim.h"
#include "pval.h"
#include "dfit_helppage.h"
#include "main_print.h"

extern htsFormat *dingding2;

int helppage_aggregate(FILE *fp){
  fprintf(fp,"Aggregation of lca produced statistics (mean length, variance length, mean GC, variance GC) when transversing up the nodes of the tree structure\n");

  fprintf(stderr, "\t\t./metaDMG-cpp aggregate file.bdamage.gz --names file.gz --nodes trestructure.gz --lcastat file.stat --out filename\n");
  fprintf(fp,"--help \t\t Print extended help page to see all options.\n");
  fprintf(fp,"--names \t names.dmp.gz\n");
  fprintf(fp,"--nodes \t nodes.dmp.gz\n");
  fprintf(fp,"--lcastat \t\t lcaout.stat lca produced statistics\n");
  fprintf(fp,"[--dfit] \t\t output from dfit function. Optional\n");
  fprintf(fp,"[--rlens] \t\t input rlens.gz to aggregate up taxonomy tree. Optional\n");
  fprintf(fp,"--out \t\t Suffix of outputname with the predetermined prefix (.stat.gz)\n");

  return 0;
}
/*
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

#if 0
    double variance1,variance2;
    int nreads1, nreads2;
#endif
    if(((double) md1.nreads+md2.nreads)>2){//pooled variance of length and GC
$if 0
    variance1 = ((double) md1.nreads-1)*md1.data[1];
      nreads1 += md1.nreads;
      variance2 = ((double) md2.nreads-1)*md2.data[1];
      nreads2 += md2.nreads;
      md2.data[1] = (variance1 + variance2) / (nreads1 + nreads2 - 2);
      md2.data[3] = (variance1 + variance2) / (nreads1 + nreads2 - 2);
      #endif
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
*/
static void merge_mean_var_pairs(double *dst,int old_nreads,const double *src,int src_nreads,int nstats){
  int new_nreads = old_nreads + src_nreads;
  for(int i=0;i+1<nstats;i+=2){
    double old_mean = dst[i];
    double old_var  = dst[i+1];
    double src_mean = src[i];
    double src_var  = src[i+1];

    dst[i] = ((double)src_mean*src_nreads + old_mean*old_nreads)/((double)src_nreads + old_nreads);

    if(new_nreads>1){
      if(src_nreads>1 && old_nreads>1){
        double ss =
          ((double)src_nreads-1)*src_var +
          ((double)old_nreads-1)*old_var +
          ((double)src_nreads*old_nreads/(double)new_nreads)*
          (src_mean-old_mean)*(src_mean-old_mean);
        dst[i+1] = ss/((double)new_nreads-1);
      }
      else if(src_nreads>1 && old_nreads<=1){
        dst[i+1] = src_var;
      }
      else if(src_nreads<=1 && old_nreads>1){
        dst[i+1] = old_var;
      }
    }
  }
}

void to_root(int from,int to,std::map<int,mydata2> &stats,int2int &parent,int src_nreads,double *src_data,int nstats){
  // fprintf(stderr,"from: %d to: %d src_nreads:%d\n",from,to,src_nreads);

  std::map<int,mydata2>::iterator it=stats.find(to);
  if(it==stats.end()){
    mydata2 mdmis;
    mdmis.nreads = src_nreads;

    mdmis.data = new double[nstats];
    for(int iii=0;iii<nstats;iii++)
      mdmis.data[iii] = src_data[iii];

    stats[to] = mdmis;
  }
  else{
    mydata2 &md2 = stats.find(to)->second;

    int old_nreads = md2.nreads;
    merge_mean_var_pairs(md2.data, old_nreads, src_data, src_nreads, nstats);

    md2.nreads += src_nreads;
  }

  int2int::iterator parent_it = parent.find(to);
  if(parent_it==parent.end()){
    fprintf(stderr,"\t-> Error: parent for taxid:%d not found while traversing tree, will exit\n",to);
    exit(1);
  }
  int newto = parent_it->second;

  if(newto!=to)
    to_root(from,newto,stats,parent,src_nreads,src_data,nstats);

}

std::map<int,char *> read_dfit(char *fname){
  //  fprintf(stderr,"fname: %s\n",fname);
  std::map<int, char *> ret;
  
  BGZF *fpfpfp = NULL;
  fpfpfp = bgzf_open(fname, "rb");
  if(fpfpfp==NULL){
    fprintf(stderr,"\t-> Problem opening file: %s will exit\n",fname);
    exit(1);
  }
  
  kstring_t *kstr = new kstring_t;
  kstr->s=NULL;kstr->l=kstr->m = 0;
  while(bgzf_getline(fpfpfp,'\n',kstr)){
    if(kstr->l==0)
      break;
    //    fprintf(stderr,"%s len:%d",kstr->s,kstr->l);
    char *taxid = kstr->s;
    char *firsttab = strchr(kstr->s,'\t');
    if(firsttab==NULL){
      fprintf(stderr,"\t-> Warning: malformed dfit line without tab, skipping line: %s\n",kstr->s);
      kstr->l = 0;
      continue;
    }
    char *sectab = firsttab+1;
    firsttab[0] = '\0';
    //    fprintf(stderr,"taxid: %d\nstr: %s\n",atoi(taxid),kstr->s);
    if(ret.size()==0)
      ret[-1] = strdup(sectab);
    else
      ret[atoi(taxid)] = strdup(sectab);
    //    if(ret.size()>3)       break;
    kstr->l = 0;
    
  }
  fprintf(stderr,"\t-> Done reading file: \"%s\", contains: %lu taxid\n",fname,ret.size());
  free(kstr->s);
  delete kstr;
  bgzf_close(fpfpfp);
  return ret;
}

typedef std::map<int, std::vector<uint64_t> > rlensmap_t;

static void rlens_add_vec(std::vector<uint64_t> &dst, const std::vector<uint64_t> &src){
  if(dst.size()<src.size())
    dst.resize(src.size(),0);
  for(size_t i=0;i<src.size();i++)
    dst[i] += src[i];
}

static rlensmap_t read_rlens(const char *fname){
  rlensmap_t ret;
  BGZF *fp = bgzf_open(fname, "rb");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem opening rlens file: %s will exit\n",fname);
    exit(1);
  }

  kstring_t *kstr = new kstring_t;
  kstr->s=NULL;kstr->l=kstr->m=0;
  while(bgzf_getline(fp,'\n',kstr)>=0){
    if(kstr->l==0){
      kstr->l = 0;
      continue;
    }
    if(!strncmp(kstr->s,"id\t",3) || !strncmp(kstr->s,"ID\t",3)){
      kstr->l = 0;
      continue;
    }

    char *saveptr = NULL;
    char *tok = strtok_r(kstr->s,"\t\n\r ",&saveptr);
    if(tok==NULL){
      kstr->l = 0;
      continue;
    }

    int taxid = atoi(tok);
    std::vector<uint64_t> &hist = ret[taxid];

    // Sparse format: id \t rlen:count \t rlen:count ...
    char *next = strtok_r(NULL,"\t\n\r ",&saveptr);
    if(next && strchr(next,':')){
      do{
        char *colon = strchr(next,':');
        if(colon==NULL)
          continue;
        *colon = '\0';
        long rlen = atol(next);
        unsigned long long count = strtoull(colon+1,NULL,10);
        if(rlen<0)
          continue;
        if(hist.size() <= (size_t)rlen)
          hist.resize((size_t)rlen+1,0);
        hist[(size_t)rlen] += (uint64_t)count;
      }while((next = strtok_r(NULL,"\t\n\r ",&saveptr))!=NULL);
    }
    // Flat format: ID \t rlen \t count
    else if(next){
      char *next2 = strtok_r(NULL,"\t\n\r ",&saveptr);
      if(next2){
        long rlen = atol(next);
        unsigned long long count = strtoull(next2,NULL,10);
        if(rlen>=0){
          if(hist.size() <= (size_t)rlen)
            hist.resize((size_t)rlen+1,0);
          hist[(size_t)rlen] += (uint64_t)count;
        }
      }
    }
    kstr->l = 0;
  }

  fprintf(stderr,"\t-> Done reading rlens file: \"%s\", contains: %lu ids\n",fname,ret.size());
  free(kstr->s);
  delete kstr;
  bgzf_close(fp);
  return ret;
}

static void aggregate_rlens_to_root(rlensmap_t &rlens, int2int &parent){
  rlensmap_t source = rlens;
  for(rlensmap_t::iterator it=source.begin();it!=source.end();it++){
    int current = it->first;
    while(1){
      int2int::iterator pit = parent.find(current);
      if(pit==parent.end()){
        fprintf(stderr,"\t-> Error: parent for taxid:%d missing while aggregating rlens, will exit\n",current);
        exit(1);
      }
      int up = pit->second;
      if(up==current)
        break;
      rlens_add_vec(rlens[up], it->second);
      current = up;
    }
  }
}

static void write_rlens(const char *fname, const rlensmap_t &rlens){
  BGZF *fp = bgzf_open(fname, "wb");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem opening output rlens file: %s will exit\n",fname);
    exit(1);
  }
  kstring_t *kstr = new kstring_t;
  kstr->s=NULL;kstr->l=kstr->m=0;

  ksprintf(kstr,"id\trlen:count\n");
  for(rlensmap_t::const_iterator it=rlens.begin();it!=rlens.end();it++){
    ksprintf(kstr,"%d",it->first);
    for(size_t i=0;i<it->second.size();i++){
      if(it->second[i]>0)
        ksprintf(kstr,"\t%zu:%llu",i,(unsigned long long)it->second[i]);
    }
    ksprintf(kstr,"\n");
    if(kstr->l>1000000){
      if(bgzf_write(fp,kstr->s,kstr->l)!=(ssize_t)kstr->l){
        fprintf(stderr,"\t-> Problem writing output rlens file: %s\n",fname);
        exit(1);
      }
      kstr->l = 0;
    }
  }
  if(kstr->l>0){
    if(bgzf_write(fp,kstr->s,kstr->l)!=(ssize_t)kstr->l){
      fprintf(stderr,"\t-> Problem writing output rlens file: %s\n",fname);
      exit(1);
    }
  }
  free(kstr->s);
  delete kstr;
  bgzf_close(fp);
}

void aggr_stat3000(std::map<int, mydata2> &stats,int2int &parent,int nstats){
  std::map<int,int> dasmap;
  std::map<int,double *> datamap;

  for(std::map<int,mydata2>::iterator it = stats.begin();it!=stats.end();it++){
    dasmap[it->first] = it->second.nreads;

    datamap[it->first] = new double[nstats];
    for(int i=0;i<nstats;i++)
      datamap[it->first][i] = it->second.data[i];
  }

  for(std::map<int,int>::iterator itt=dasmap.begin();itt!=dasmap.end();itt++){
    int focal_taxid = itt->first;
    int2int::iterator it = parent.find(focal_taxid);
    if (it == parent.end()) {
      fprintf(stderr, "\t-> Error: iterator reached end (key:%d not found),will exit\n",focal_taxid);
      exit(1);
    }
    int target = it->second;

    if(target!=focal_taxid)
      to_root(focal_taxid,target,stats,parent,itt->second,datamap[focal_taxid],nstats);
  }

  for(std::map<int,double *>::iterator it=datamap.begin();it!=datamap.end();it++)
    delete [] it->second;
}
int main_aggregate(int argc, char **argv) {
    if (argc <= 1){
      helppage_aggregate(stderr);
      return 0;
    }
    char *infile_bdamage = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;
    char *infile_lcastat = NULL;
    char *outfile_name = NULL;
    char *infile_dfit = NULL;
    char *infile_rlens = NULL;
    int howmany;//this is the cycle
    while (*(++argv)) {
      if (strcasecmp("-h", *argv) == 0 || strcasecmp("--help", *argv) == 0)
	return helppage_aggregate(stderr);
      else if (strcasecmp("--names", *argv) == 0 || strcasecmp("-names", *argv) == 0)
	infile_names = strdup(*(++argv));
      else if (strcasecmp("--nodes", *argv) == 0 || strcasecmp("-nodes", *argv) == 0)
	infile_nodes = strdup(*(++argv));
      else if (strcasecmp("--dfit", *argv) == 0 || strcasecmp("-dfit", *argv) == 0)
	infile_dfit = strdup(*(++argv));
      else if (strcasecmp("--rlens", *argv) == 0 || strcasecmp("-rlens", *argv) == 0)
	infile_rlens = strdup(*(++argv));
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
    if(infile_bdamage==NULL){
      fprintf(stderr,"\t-> bdamage input file is required\n");
      return helppage_aggregate(stderr);
    }
    fprintf(stderr,
	    "aggregate infile_bdamage: %s infile_names: %s infile_nodes: %s infile_lcastat: %s infile_dfit: %s outfile_name: %s\n",
	    infile_bdamage ? infile_bdamage : "NULL",
	    infile_names   ? infile_names   : "NULL",
	    infile_nodes   ? infile_nodes   : "NULL",
	    infile_lcastat ? infile_lcastat : "NULL",
	    infile_dfit    ? infile_dfit    : "NULL",
	    outfile_name   ? outfile_name   : "NULL"
	    );
    if(infile_rlens && infile_nodes==NULL){
      fprintf(stderr,"\t-> Warning: --rlens was supplied without --nodes, rlens output will not be taxonomically expanded\n");
    }
    if(outfile_name==NULL)
      outfile_name = strdup(infile_bdamage);
    fprintf(stderr,
	    "aggregate infile_bdamage: %s infile_names: %s infile_nodes: %s infile_lcastat: %s infile_dfit: %s outfile_name: %s\n",
	    infile_bdamage ? infile_bdamage : "NULL",
	    infile_names   ? infile_names   : "NULL",
	    infile_nodes   ? infile_nodes   : "NULL",
	    infile_lcastat ? infile_lcastat : "NULL",
	    infile_dfit    ? infile_dfit    : "NULL",
	    outfile_name   ? outfile_name   : "NULL"
	    );
    char buf[1024];
    snprintf(buf, 1024, "%s.stat.gz", outfile_name);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    BGZF *fpfpfp = bgzf_open(buf, "wb");
    if(fpfpfp==NULL){
      fprintf(stderr,"\t-> Cannot open output BGZ file: %s\n",buf);
      exit(1);
    }
    kstring_t *kstr = new kstring_t;
    kstr->s = NULL; kstr->l = kstr->m = 0;
  
    
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
    if(child.size()>0){
	std::map<int,mydataD> results = getval_full_norec(retmap,parent,howmany);//lizard king 2000.
	retmap = results;
    }
    float postsize = retmap.size();
    fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);

    int stat_dims = 4;
    std::map<int, mydata2> stats;
    if (infile_lcastat){
      stats = load_lcastat(infile_lcastat,1,&stat_dims);
      presize = stats.size();
    }
    ksprintf(kstr, "taxid\tname\trank\tnalign\tnreads\tmean_rlen\tvar_rlen\tmean_gc\tvar_gc");
    if(stat_dims>=8)
      ksprintf(kstr, "\tmean_dust\tvar_dust\tmean_nspec\tvar_nspec");
    ksprintf(kstr, "\tlca\ttaxa_path");
    
    std::map<int, char *> dfit_int_char;
    if(infile_dfit!=NULL){
      dfit_int_char = read_dfit(infile_dfit);
      //      fprintf(stderr,"Donedone\n");
    }
    if(dfit_int_char.size()>0){
      std::map<int, char *>::iterator it = dfit_int_char.find(-1);
      if (it == dfit_int_char.end()) {
	fprintf(stderr, "\t-> Error: iterator reached end in dfit_int_char, will exit\n");
	exit(1);
      }
      ksprintf(kstr,"\t%s",it->second);
    }
    ksprintf(kstr,"\n");

#if 0
    if(child.size()>0)//this will not work if we have data at internal nodes 
      getval_stats(stats, child, 1);  // this will do everything
#endif
    aggr_stat3000(stats,parent,stat_dims);
    //    exit(0);
    if(stats.size()>0){
      postsize = stats.size();
      fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);
    }

    if(infile_rlens){
      rlensmap_t rlens = read_rlens(infile_rlens);
      if(child.size()>0)
        aggregate_rlens_to_root(rlens,parent);
      char rbuf[1024];
      snprintf(rbuf, 1024, "%s.rlens.gz", outfile_name);
      write_rlens(rbuf, rlens);
      fprintf(stderr,"\t-> Dumped aggregated rlens file: '%s' with %lu ids\n",rbuf,rlens.size());
    }

    for (std::map<int, mydata2>::iterator it = stats.begin(); it != stats.end(); it++) {
      std::map<int, mydataD>::iterator itold = retmap.find(it->first);
      std::map<int, char*>::iterator itchar = dfit_int_char.find(it->first);
        int nalign = -1;
        if (itold == retmap.end()) {
            fprintf(stderr, "\t-> Problem finding taxid: %d\n", it->first);
        } else
            nalign = itold->second.nal;
        char *myrank = NULL;
        char *myname = NULL;
        if (it->second.nreads > 0) {
            int2char::iterator itc = rank.find(it->first);
            if (itc != rank.end())
                myrank = itc->second;
            itc = name_map.find(it->first);
            if (itc != name_map.end())
                myname = itc->second;
            const char *safe_name = myname ? myname : "NA";
            const char *safe_rank = myrank ? myrank : "NA";
            ksprintf(kstr, "%d\t\"%s\"\t\"%s\"\t%d\t%d\t%f\t%f\t%f\t%f", it->first, safe_name, safe_rank, nalign, it->second.nreads, it->second.data[0], it->second.data[1], it->second.data[2], it->second.data[3]);
            if(stat_dims>=8)
              ksprintf(kstr,"\t%f\t%f\t%f\t%f",it->second.data[4],it->second.data[5],it->second.data[6],it->second.data[7]);
	    if(child.size()>0)
	      print_chain(kstr, it->first, parent, rank, name_map,0);
	    else
	      ksprintf(kstr,"NA\tNA");
            //      fprintf(stderr,"%d->(%d,%f,%f,%f,%f)\n",it->first,it->second.nreads,it->second.data[0],it->second.data[1],it->second.data[2],it->second.data[3]);
	    if(itchar!=dfit_int_char.end())
	      ksprintf(kstr,"\t%s",itchar->second);
	    ksprintf(kstr,"\n");
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
    for(std::map<int,char *>::iterator it=dfit_int_char.begin();it!=dfit_int_char.end();it++)
      free(it->second);
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
    if(infile_dfit)
      free(infile_dfit);
    if(infile_rlens)
      free(infile_rlens);
    if(outfile_name)
      free(outfile_name);

  return 0;
}
