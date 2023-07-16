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

extern htsFormat *dingding2;

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
    snprintf(buf, 1024, "%s.uglyprint.mismatch.txt.gz", infile_bdamage);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    gzFile fpfpfp = gzopen(buf, "wb");
    gzprintf(fpfpfp, "#taxidStr\tdirection\tposition\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
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

    for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        mydataD md = it->second;
        if (it->second.nreads == 0)
            continue;
        /*
        char *myrank =NULL;
        char *myname = NULL;
        int2char::iterator itc=rank.find(taxid);
        if(itc!=rank.end())
          myrank=itc->second;
          itc=name.find(taxid);
        if(itc!=name.end())
          myname = itc->second;
        */
        for (int i = 0; i < howmany; i++) {
            //      fprintf(stdout,"%d\t\"%s\"\t\"%s\"\t%d\t5'\t%d",taxid,myname,myrank,it->second.nreads,i);
            //      fprintf(fpfpfp,"%d\t%d\t5'\t%d",taxid,it->second.nreads,i);
            if (hdr == NULL)
                gzprintf(fpfpfp, "%d\t5'\t%d", taxid, i);
            else
                gzprintf(fpfpfp, "%s\t5'\t%d", sam_hdr_tid2name(hdr, taxid), i);

            for (int ii = 0; ii < 16; ii++)
                gzprintf(fpfpfp, "\t%.0f", it->second.fwD[i * 16 + ii]);
            gzprintf(fpfpfp, "\n");
        }
        for (int i = 0; i < howmany; i++) {
            //      fprintf(stdout,"%d\t\"%s\"\t\"%s\"\t%d\t3'\t%d",taxid,myname,myrank,it->second.nreads,i);
            // fprintf(fpfpfp,"%d\t%d\t3'\t%d",taxid,it->second.nreads,i);
            if (hdr == NULL)
                gzprintf(fpfpfp, "%d\t3'\t%d", taxid, i);
            else
                gzprintf(fpfpfp, "%s\t3'\t%d", sam_hdr_tid2name(hdr, taxid), i);
            for (int ii = 0; ii < 16; ii++)
                gzprintf(fpfpfp, "\t%.0f", it->second.bwD[i * 16 + ii]);
            gzprintf(fpfpfp, "\n");
        }
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
