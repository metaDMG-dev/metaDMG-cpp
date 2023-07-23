/*
  16 july, this file contains functionality from the print_ugly function and codeparts in the misc.
  The function will implement the ML estimate of the dfit from metaDMG-cpp bioxarhive paper.
  
 */

#include <cstdio>
#include <htslib/hts.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/sam.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/bgzf.h>

#include "profile.h"
#include "shared.h"
#include "ngsLCA.h" //<- print_chain
#include "types.h"       // for int2intvec, int2int
#include "version.h"     // for METADAMAGE_VERSION
#include "dfit_optim.h"
extern htsFormat *dingding2;

//aa,ac,ag,at,ca,cc,cg,ct,ga
//ct and ga has index 7,8 when zero indexed
void make_dfit_format(mydataD &md,double **dat,int howmany){
  dat[0][0] = 2*howmany;
  dat[0][1] = 0;
  for(int i=0;i<howmany;i++){
    dat[0][i+2] =i;

    dat[2][i] = 0;//initialize kcol to zero, a few lines down we will sum over all C*
    dat[1][i] = md.fwD[i*16+7];//plugin ct at kcol
    
    for(int at=0;at<4;at++)
      dat[2][i] = dat[2][i]+md.fwD[i*16+4+at];
  }

  for(int i=0;i<howmany;i++){
    dat[0][howmany+i+2] =i;

    dat[2][howmany+i] = 0;//initialize kcol to zero, a few lines dows we will sum over all G*
    dat[1][howmany+i] = md.bwD[i*16+8];//plugin ga at kcol
    
    for(int at=0;at<4;at++)
      dat[2][howmany+i] = dat[2][howmany+i]+md.bwD[i*16+8+at];
  }

}

mydataD getval_full(std::map<int, mydataD> &retmap, int2intvec &child, int taxid, int howmany);
mydata2 getval_stats(std::map<int, mydata2> &retmap, int2intvec &child, int taxid) ;

int main_dfit(int argc, char **argv) {
    fprintf(stderr, "./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz -bam file.bam -showfits int -nopt int -out file\n");
    if (argc <= 1)
        return 0;
    char *infile_bdamage = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;
    char *infile_lcastat = NULL;
    char *infile_bam = NULL;
    char *outfile_name = NULL;
    int howmany;//this is the cycle
    int showfits=0;
    int nopt = 5;
    while (*(++argv)) {
        if (strcasecmp("-names", *argv) == 0)
            infile_names = strdup(*(++argv));
        else if (strcasecmp("-nodes", *argv) == 0)
            infile_nodes = strdup(*(++argv));
        else if (strcasecmp("-lcastat", *argv) == 0)
            infile_lcastat = strdup(*(++argv));
        else if (strcasecmp("-bam", *argv) == 0)
            infile_bam = strdup(*(++argv));
	 else if (strcasecmp("-out", *argv) == 0)
            outfile_name = strdup(*(++argv));
	else if (strcasecmp("-nopt", *argv) == 0)
	  nopt = atoi(*(++argv));
	else if (strcasecmp("-showfits", *argv) == 0)
	  showfits = atoi(*(++argv));
        else
            infile_bdamage = strdup(*argv);
    }
    if(infile_nodes&&!infile_names){
      fprintf(stderr,"\t-> -names file.txt.gz must be defined with -nodes is defined\n");
      exit(1);
    }
    htsFile *samfp = NULL;
    sam_hdr_t *hdr = NULL;
    if (infile_bam) {
        if ((samfp = sam_open_format(infile_bam, "r", dingding2)) == NULL) {
            fprintf(stderr, "[%s] nonexistant file: %s\n", __FUNCTION__, infile_bam);
            exit(1);
        }
        hdr = sam_hdr_read(samfp);
    }
    if(infile_bam!=NULL&&infile_nodes!=NULL){
      fprintf(stderr,"\t-> Potential issue, supplying both -bam and -nodes might not be meaningfull\n");
      exit(1);
    }
    
    fprintf(stderr, "infile_names: %s infile_bdamage: %s nodes: %s lca_stat: %s infile_bam: %s showfits: %d nopt: %d outname: %s ", infile_names, infile_bdamage, infile_nodes, infile_lcastat, infile_bam,showfits,nopt,outfile_name);
    fprintf(stderr, "#VERSION:%s\n", METADAMAGE_VERSION);
    if(outfile_name==NULL)
      outfile_name = strdup(infile_bdamage);
    char buf[1024];
    snprintf(buf, 1024, "%s.dfit.txt.gz", outfile_name);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    BGZF *fpfpfp = bgzf_open(buf, "wb");
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
    if(child.size()>0)
      getval_full(retmap, child, 1, howmany);  // this will do everything
    float postsize = retmap.size();
    fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);

    //prepare matrix to be passed to optimization
    double **dat = new double*[3];
    dat[0] = new double[2+2*howmany];
    dat[1] = new double [2*howmany];
    dat[2] = new double [2*howmany];
    fprintf(stderr,"\t-> Will do optimization of %lu different taxids/chromosomes/scaffolds\n",retmap.size());
    if(showfits==0){
      ksprintf(kstr,"id\tA\tq\tc\tphi\tllh\tncall\tZ\tsign\n");
    }else{
      ksprintf(kstr,"id\tA\tq\tc\tphi\tllh\tncall\tZ\tsign");
      for(int i=0;i<howmany;i++)
	ksprintf(kstr,"\tfwK%d\tfwN%d\tfwdx%d\tfwdxConf%d",i,i,i,i);
      for(int i=0;i<howmany;i++)
	ksprintf(kstr,"\tbwK%d\tbwN%d\tbwdx%d\tbwtdxConf%d",i,i,i,i);
      ksprintf(kstr,"\n");
    }

    for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        mydataD md = it->second;
        if (it->second.nreads == 0)
            continue;

	if(hdr!=NULL){
	  ksprintf(kstr, "%s\t", sam_hdr_tid2name(hdr, it->first));
	}else{
	  if(name_map.size()==0)
	    ksprintf(kstr,"%d\t",it->first);
	  else{
	    int2char::iterator nit = name_map.find(it->first);
	    if(nit==name_map.end()){
	      fprintf(stderr,"\t-> Problem finding taxid: %d in names database: %s\n",it->first,infile_names);
	      exit(1);
	    }
	    ksprintf(kstr,"%d:%s\t",it->first,nit->second);
	  }
	    
	}

	make_dfit_format(md,dat,howmany);
	double pars[6] = {0.1,0.1,0.01,1000};//last one will contain the llh,and the ncall for the objective function
	optimoptim(pars,dat,5);
	double stats[2+2*(int)dat[0][0]];
	getstat(dat,pars,stats);
	
	//printit

	for(int i=0;i<6;i++)
	  ksprintf(kstr,"%f\t",pars[i]);
	ksprintf(kstr,"%f\t%f",stats[0],stats[1]);
	if(showfits){
	  int nrows = (int) dat[0][0];
	  int ncycle = nrows/2;
	  double *dx = stats+2;
	  double *dx_conf = stats+2+nrows;
	  for(int i=0;i<ncycle;i++)
	    ksprintf(kstr,"\t%.0f\t%0.f\t%f\t%f",dat[1][i],dat[2][i],dx[i],dx_conf[i]);
	  dx = stats+2+ncycle;
	  dx_conf = stats+2+nrows+ncycle;
	  for(int i=0;i<ncycle;i++)
	    ksprintf(kstr,"\t%.0f\t%0.f\t%f\t%f",dat[1][i+ncycle],dat[2][i+ncycle],dx[i],dx_conf[i]);
	}
	ksprintf(kstr,"\n");
    }
    bgzf_write(fpfpfp,kstr->s,kstr->l);
    kstr->l = 0;
    bgzf_close(fpfpfp);
    fpfpfp = NULL;
    if(infile_lcastat){
      snprintf(buf, 1024, "%s.dfit.stat.txt.gz", outfile_name);
      fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
      fpfpfp = bgzf_open(buf, "wb");
      ksprintf(kstr, "#taxid\tname\trank\tnalign\tnreads\tmean_rlen\tvar_rlen\tmean_gc\tvar_gc\n");
    }
    std::map<int, mydata2> stats;
    if (infile_lcastat){
        stats = load_lcastat(infile_lcastat);
	presize = stats.size();
    }
    if(child.size()>0)
      getval_stats(stats, child, 1);  // this will do everything
    if(stats.size()>0){
      postsize = stats.size();
      fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);
    }
    for (std::map<int, mydata2>::iterator it = stats.begin(); it != stats.end(); it++) {
      std::map<int, mydataD>::iterator itold = retmap.find(it->first);
        int nalign = -1;
        if (itold == retmap.end()) {
            fprintf(stderr, "\t-> Problem finding taxid: %d\n", it->first);
	    exit(1);
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
            ksprintf(kstr, "%d\t\"%s\"\t\"%s\"\t%d\t%d\t%f\t%f\t%f\t%f\t", it->first, myname, myrank, nalign, it->second.nreads, it->second.data[0], it->second.data[1], it->second.data[2], it->second.data[3]);
	    if(child.size()>0)
	      print_chain(kstr, it->first, parent, rank, name_map);
	    else
	      ksprintf(kstr,"\n");
            //      fprintf(stderr,"%d->(%d,%f,%f,%f,%f)\n",it->first,it->second.nreads,it->second.data[0],it->second.data[1],it->second.data[2],it->second.data[3]);
        }
    }
    //cleanup
    if(fpfpfp){
      bgzf_write(fpfpfp,kstr->s,kstr->l);
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
    
    if (hdr)
        bam_hdr_destroy(hdr);
    if (samfp)
        sam_close(samfp);
    if(infile_bdamage)
      free(infile_bdamage);
    if(infile_nodes)
      free(infile_nodes);
    if(infile_names)
      free(infile_names);
    if(infile_bam)
      free(infile_bam);
    if(infile_lcastat)
      free(infile_lcastat);
    return 0;
}
