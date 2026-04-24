// gpl thorfinn@binf.ku.dk
#include <getopt.h>     // for optarg, getopt_long, optind, option
#include <htslib/hts.h> // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/sam.h> // for sam_hdr_read, bam_hdr_destroy, bam_hdr_t
#include <zlib.h>       // for gzprintf, gzclose, gzopen

#include <cstdio>   // for fprintf, stderr, stdout
#include <cstdlib>  // for atoi, exit, free, calloc, malloc
#include <cstring>  // for strdup, strlen
#include <map>      // for map
#include <vector>   // for vector

#include "main_getdamage.h"
#include "ngsLCA.h"  // for gccontent, mean, var
#include "profile.h" // for damage, destroy_damage, triple

static int usage_getdamage(FILE *fp) {
    fprintf(fp, "\nUsage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>\n");
    fprintf(fp, "\nExample: ./metaDMG-cpp getdamage -l 10 -p 5 --threads 8 ../data/subs.sam\nOptions:\n");
    fprintf(fp, "  -n/--threads\t\t number of threads used for reading/writing (default: 4)\n");
    fprintf(fp, "  -f/--fasta\t\t reference genome (required with CRAM)\n");
    fprintf(fp, "  -l/--min_length\t minimum read length (default: 35)\n");
    fprintf(fp, "  -p/--print_length\t number of base pairs from read termini to estimate damage (default: 5)\n");
    fprintf(fp, "  -r/--run_mode\t\t 0: global estimate (default)\n\t\t\t 1: damage patterns will be calculated for each chr/scaffold contig\n");
    fprintf(fp, "  -i/--ignore_errors\t continue analyses even if there are errors.\n");
    fprintf(fp, "  -o/--out_prefix\t output prefix (default: meta)\n");
    fprintf(fp, "  -z/--rlens_flat_out\t make flat output of bins. Nice for computers\n");
    return 1;
}

static int write_getdamage_stats(char *onam,
                                 int runmode,
                                 bam_hdr_t *hdr,
                                 damage *dmg,
                                 std::map<int, std::vector<float> > &gcconts,
                                 std::map<int, std::vector<float> > &seqlens) {
    char buf[1024];
    snprintf(buf, 1024, "%s.stat.gz", onam);
    fprintf(stderr, "\t-> Outputting overall statistic in file: \"%s\"\n", buf);

    gzFile fpstat = NULL;
    if ((fpstat = gzopen(buf, "wb")) == NULL) {
      fprintf(stderr, "\t-> Error problem opening file %s will exit\n", buf);
      return 1;
    }
    gzprintf(fpstat, "taxid\tnreads\tmean_len\tvar_len\tmean_gc\tvar_gc\tlca\trank\n");
    for (std::map<int, std::vector<float> >::iterator it = gcconts.begin(); it != gcconts.end(); it++) {
        std::map<int, triple>::iterator it2 = dmg->assoc.find(it->first);
        if (it2 == dmg->assoc.end()) {
          fprintf(stderr, "\t-> Error: iterator reached end in dmg->assoc, will exit\n");
          gzclose(fpstat);
          return 1;
        }
        std::map<int, std::vector<float> >::iterator it3 = seqlens.find(it->first);
        if (it3 == seqlens.end()) {
          fprintf(stderr, "\t-> Error: iterator reached end in seqlens (it3), will exit\n");
          gzclose(fpstat);
          return 1;
        }
        if (runmode == 1)
          gzprintf(fpstat, "%s\t%lu\t%f\t%f\t%f\t%f\tNA\tNA\n", sam_hdr_tid2name(hdr, it->first), it2->second.nreads, mean(it3->second), var(it3->second), mean(it->second), var(it->second));
        else
          gzprintf(fpstat, "global\t%lu\t%f\t%f\t%f\t%f\tNA\tNA\n", it2->second.nreads, mean(it3->second), var(it3->second), mean(it->second), var(it->second));
    }

    gzclose(fpstat);
    return 0;
}

int main_getdamage(int argc, char **argv) {
    if (argc == 1)
        return usage_getdamage(stderr);

    //  int MAXLENGTH = 256;
    int rc = 0;
    int minLength = 35;
    int printLength = 5;
    char *refName = NULL;
    char *fname = NULL;
    int runmode = 0;  // this means one species, runmode=1 means multi species
    htsFile *fp = NULL;
    bam1_t *b = NULL;
    bam_hdr_t *hdr = NULL;
    damage *dmg = NULL;
    int skipper[4] = {3, 3, 3, 3};
    std::map<int, std::vector<float> > gcconts;
    std::map<int, std::vector<float> > seqlens;
    char *onam = strdup("meta");
    int nthreads = 4;
    int ignore_errors = 0;
    int rlens_flat_out = 0;
    htsFormat *dingding2 = (htsFormat *)calloc(1, sizeof(htsFormat));
    // fix thesepro
    static struct option lopts[] = {
        {"threads", required_argument, 0, 'n'},
        {"rlens_flat_out", required_argument, 0, 'z'},
        {"fasta", required_argument, 0, 'f'},
        {"min_length", required_argument, 0, 'l'},
        {"print_length", required_argument, 0, 'p'},
        {"run_mode", required_argument, 0, 'r'},
        {"ignore_errors", no_argument, 0, 'i'},
        {"out_prefix", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}};

    int c;
    while ((c = getopt_long(argc, argv,
                            "iz:n:f:l:p:r:o:h",
                            lopts, NULL)) >= 0) {
      switch (c) {
      case 'i':
        ignore_errors = 1;
        break;
      case 'n':
        nthreads = atoi(optarg);
        break;
      case 'f':
        refName = strdup(optarg);
        break;
      case 'z':
        rlens_flat_out = atoi(optarg);
        break;
      case 'l':
        minLength = atoi(optarg);
        break;
      case 'p':
        printLength = atoi(optarg);
        break;
      case 'r':
        runmode = atoi(optarg);
        break;
      case 'o':
        free(onam);
        onam = strdup(optarg);
        break;
      case 'h':
        rc = usage_getdamage(stdout);
        goto cleanup;
      default:
        fprintf(stderr, "Never here: %s\n", optarg);
        break;
      }
    }
    if (optind < argc)
        fname = strdup(argv[optind]);
    fprintf(stderr, "\t-> ./metaDMG-cpp refName: %s min_length: %d print_length: %d run_mode: %d out_prefix: %s nthreads: %d ignore_errors: %d rlens_flat_out: %d\n", refName ? refName : "NULL", minLength, printLength, runmode, onam, nthreads, ignore_errors, rlens_flat_out);
    if (fname == NULL) {
        rc = usage_getdamage(stderr);
        goto cleanup;
    }
    if (refName) {
      int lenlen = 10 + strlen(refName) + 1;
      char *ref = (char *)malloc(lenlen);
      snprintf(ref, lenlen, "reference=%s", refName);
      hts_opt_add((hts_opt **)&dingding2->specific, ref);
      free(ref);
    }

    if ((fp = sam_open_format(fname, "r", dingding2)) == NULL) {
        fprintf(stderr, "[%s] nonexistant file: %s\n", __FUNCTION__, fname);
        rc = 1;
        goto cleanup;
    }

    b = bam_init1();
    hdr = sam_hdr_read(fp);
    if (hdr == NULL) {
      fprintf(stderr, "\t-> Hello Doctor! im afraid there is an error reading the header\n");
      rc = 1;
      goto cleanup;
    }
    int checkIfSorted(char *str);
    if (checkIfSorted(hdr->text)) {
      fprintf(stderr, "Input alignment file is not sorted.");
      if (!ignore_errors) {
        rc = 1;
        goto cleanup;
      }
    }
    int ret;
    dmg = new damage(printLength, nthreads, 0);
    while (((ret = sam_read1(fp, hdr, b))) >= 0) {
        if (bam_is_unmapped(b)) {
            if (skipper[0])
                fprintf(stderr, "skipping: %s unmapped, this msg is printed: %d times more\n", bam_get_qname(b), --skipper[0]);
            continue;
        }
        if (bam_is_failed(b)) {
            if (skipper[1])
                fprintf(stderr, "skipping: %s failed: flags=%d, this msg is printed: %d times more\n", bam_get_qname(b), b->core.flag, --skipper[1]);
            continue;
        }
        if (b->core.l_qseq < minLength) {
            if (skipper[2])
                fprintf(stderr, "skipping: %s too short, this msg is printed %d times more \n", bam_get_qname(b), --skipper[2]);
            continue;
        }

        dmg->damage_analysis(b, runmode != 0 ? b->core.tid : 0, 1);

        float mygc = gccontent(b);
        float mylen = b->core.l_qseq;

        int whichref = 0;
        if (runmode == 1)//<- runmode 0 means local one means global
            whichref = b->core.tid;
        std::map<int, std::vector<float> >::iterator it = gcconts.find(whichref);
        if (it == gcconts.end()) {
          std::vector<float> tmp1;
          tmp1.push_back(mygc);
          gcconts[whichref] = tmp1;
          std::vector<float> tmp2;
          tmp2.push_back(mylen);
          seqlens[whichref] = tmp2;
        } else {
          if (it->second.size() < 1000000)//<- we dont need more than a million values do we?
            it->second.push_back(mygc);
          it = seqlens.find(whichref);
          if (it == seqlens.end()) {
            fprintf(stderr, "\t-> Error: iterator reached end in seqlens, will exit\n");
            rc = 1;
            goto cleanup;
          }
          if (it->second.size() < 1000000)
            it->second.push_back(mylen);
        }
    }

    if (ret < -1) {
      fprintf(stderr, "\t-> Error: sam_read1 failed for input file: %s (ret=%d)\n", fname, ret);
      rc = 1;
      goto cleanup;
    }

    dmg->printit(stdout, printLength);
    dmg->write(onam, runmode == 1 ? hdr : NULL);
    dmg->bwrite(onam, rlens_flat_out);

    rc = write_getdamage_stats(onam, runmode, hdr, dmg, gcconts, seqlens);
    if (rc != 0)
      goto cleanup;
cleanup:
    if (hdr)
      sam_hdr_destroy(hdr);
    if (b)
      bam_destroy1(b);
    if (fp)
      sam_close(fp);
    if (dmg)
      destroy_damage(dmg);
    if (fname)
      free(fname);
    if (onam)
      free(onam);
    if (refName)
      free(refName);
    if (dingding2)
      free(dingding2);
    return rc;
}
