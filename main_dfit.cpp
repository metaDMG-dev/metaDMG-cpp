/*
  16 july, this file contains functionality from the print_ugly function and codeparts in the misc.
  The function will implement the ML estimate of the dfit from metaDMG-cpp bioxarhive paper.
  
 */

#include <cstdio>
#include <zlib.h>
#include <htslib/hts.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/sam.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt

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
    fprintf(stderr, "./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz\n");
    if (argc <= 1)
        return 0;
    char *infile_bdamage = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;
    char *infile_lcastat = NULL;
    char *infile_bam = NULL;
    int howmany;

    while (*(++argv)) {
        if (strcasecmp("-names", *argv) == 0)
            infile_names = strdup(*(++argv));
        else if (strcasecmp("-nodes", *argv) == 0)
            infile_nodes = strdup(*(++argv));
        else if (strcasecmp("-lcastat", *argv) == 0)
            infile_lcastat = strdup(*(++argv));
        else if (strcasecmp("-bam", *argv) == 0)
            infile_bam = strdup(*(++argv));
        else
            infile_bdamage = strdup(*argv);
    }

    htsFile *samfp = NULL;
    sam_hdr_t *hdr = NULL;
    if (infile_bam) {
        if ((samfp = sam_open_format(infile_bam, "r", dingding2)) == NULL) {
            fprintf(stderr, "[%s] nonexistant file: %s\n", __FUNCTION__, infile_bam);
            exit(0);
        }
        hdr = sam_hdr_read(samfp);
    }

    fprintf(stderr, "infile_names: %s infile_bdamage: %s nodes: %s lca_stat: %s infile_bam: %s", infile_names, infile_bdamage, infile_nodes, infile_lcastat, infile_bam);
    fprintf(stderr, "#VERSION:%s\n", METADAMAGE_VERSION);
    char buf[1024];
    snprintf(buf, 1024, "%s.dfit.gz", infile_bdamage);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    gzFile fpfpfp = gzopen(buf, "wb");

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
    getval_full(retmap, child, 1, howmany);  // this will do everything
    float postsize = retmap.size();
    fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);

    //prepare matrix to be passed to optimization
    double **dat = new double*[3];
    dat[0] = new double[2+2*howmany];
    dat[1] = new double [2*howmany];
    dat[2] = new double [2*howmany];

    for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        mydataD md = it->second;
        if (it->second.nreads == 0)
            continue;

	if(hdr==NULL){
	  int taxid_id_refid =1;
	  //gzprintf(fpfpfp, "%s\t5'\t%d", sam_hdr_tid2name(hdr, taxid), i);
	}

	make_dfit_format(md,dat,howmany);
	double pars[6] = {0.1,0.1,0.01,1000};//last one will contain the llh,and the ncall for the objective function
	optimoptim(pars,dat,20);
	double stats[2+2*(int)dat[0][0]];
	getstat(dat,pars,stats);
	
	//printit
	fprintf(stderr,"(A,q,c,phi,llh,ncall,Z,sign)\n");
	for(int i=0;i<6;i++)
	  fprintf(stderr,"%f\t",pars[i]);
	fprintf(stderr,"%f\t%f\n",stats[0],stats[1]);
	
	fprintf(stderr,"direction,cycle,k,n,dx,dx_conf\n");
	int nrows = (int) dat[0][0];
	int ncycle = nrows/2;
	double *dx = stats+2;
	double *dx_conf = stats+2+nrows;
	for(int i=0;i<ncycle;i++)
	  fprintf(stderr,"5\t%d\t%.0f\t%0.f\t%f\t%f\n",i,dat[1][i],dat[2][i],dx[i],dx_conf[i]);
	dx = stats+2+ncycle;
	dx_conf = stats+2+nrows+ncycle;
	for(int i=0;i<ncycle;i++)
	  fprintf(stderr,"3\t%d\t%.0f\t%0.f\t%f\t%f\n",i,dat[1][i+ncycle],dat[2][i+ncycle],dx[i],dx_conf[i]);
	
    }
    gzclose(fpfpfp);
    snprintf(buf, 1024, "%s.uglyprint.stat.txt.gz", infile_bdamage);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    fpfpfp = gzopen(buf, "wb");
    gzprintf(fpfpfp, "#taxid\tname\trank\tnalign\tnreads\tmean_rlen\tvar_rlen\tmean_gc\tvar_gc\n");
    std::map<int, mydata2> stats;
    if (infile_lcastat)
        stats = load_lcastat(infile_lcastat);
    getval_stats(stats, child, 1);  // this will do everything
    for (std::map<int, mydata2>::iterator it = stats.begin(); 1 && it != stats.end(); it++) {
        std::map<int, mydataD>::iterator itold = retmap.find(it->first);
        int nalign = -1;
        if (itold == retmap.end()) {
            fprintf(stderr, "\t-> Problem finding taxid: %d\n", it->first);
            //      exit(0);
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
            gzprintf(fpfpfp, "%d\t\"%s\"\t\"%s\"\t%d\t%d\t%f\t%f\t%f\t%f\t", it->first, myname, myrank, nalign, it->second.nreads, it->second.data[0], it->second.data[1], it->second.data[2], it->second.data[3]);
            print_chain(fpfpfp, it->first, parent, rank, name_map);
            //      fprintf(stderr,"%d->(%d,%f,%f,%f,%f)\n",it->first,it->second.nreads,it->second.data[0],it->second.data[1],it->second.data[2],it->second.data[3]);
        }
    }
    //cleanup
    gzclose(fpfpfp);
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
