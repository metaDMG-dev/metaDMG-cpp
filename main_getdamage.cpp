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

struct getdamage_args {
    int minLength;
    int printLength;
    int editdistMin;
    int editdistMax;
    double minAni;
    double maxAni;
    char *refName;
    char *fname;
    int runmode;
    char *onam;
    int nthreads;
    int ignore_errors;
    int rlens_flat_out;
};

static int usage_getdamage(FILE *fp) {
    fprintf(fp, "\nUsage: metadamage getdamage [options] <in.bam>|<in.sam>|<in.cram>\n");
    fprintf(fp, "\nExample: ./metaDMG-cpp getdamage -l 10 -p 5 --threads 8 ../data/subs.sam\nOptions:\n");
    fprintf(fp, "  -n/--threads\t\t number of threads used for reading/writing (default: 4)\n");
    fprintf(fp, "  -f/--fasta\t\t reference genome (required with CRAM)\n");
    fprintf(fp, "  -l/--min_length\t minimum read length (default: 35)\n");
    fprintf(fp, "  --edit_dist_min\t minimum read edit distance (default: -1, disabled)\n");
    fprintf(fp, "  --edit_dist_max\t maximum read edit distance (default: -1, disabled)\n");
    fprintf(fp, "  --minAni/--min_ani\t minimum ANI from NM/read_length (default: -1, disabled)\n");
    fprintf(fp, "  --maxAni/--max_ani\t maximum ANI from NM/read_length (default: -1, disabled)\n");
    fprintf(fp, "  -p/--print_length\t number of base pairs from read termini to estimate damage (default: 5)\n");
    fprintf(fp, "  -r/--run_mode\t\t 0: global estimate (default)\n\t\t\t 1: damage patterns will be calculated for each chr/scaffold contig\n");
    fprintf(fp, "  -i/--ignore_errors\t continue analyses even if there are errors.\n");
    fprintf(fp, "  -o/--out_prefix\t output prefix (default: meta)\n");
    fprintf(fp, "  -z/--rlens_flat_out\t make flat output of bins. Nice for computers\n");
    return 0;
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

static int process_getdamage_reads(htsFile *fp,
                                   bam_hdr_t *hdr,
                                   bam1_t *b,
                                   damage *dmg,
                                   int minLength,
                                   int editMin,
                                   int editMax,
                                   double minAni,
                                   double maxAni,
                                   int runmode,
                                   std::map<int, std::vector<float> > &gcconts,
                                   std::map<int, std::vector<float> > &seqlens,
                                   const char *fname) {
    int ret;
    int skipper[4] = {3, 3, 3, 3};

    while ((ret = sam_read1(fp, hdr, b)) >= 0) {
        if (bam_is_unmapped(b)) {
            if (skipper[0])
                fprintf(stderr, "skipping: %s unmapped, this msg is printed: %d times more\n",
                        bam_get_qname(b), --skipper[0]);
            continue;
        }
        if (bam_is_failed(b)) {
            if (skipper[1])
                fprintf(stderr, "skipping: %s failed: flags=%d, this msg is printed: %d times more\n",
                        bam_get_qname(b), b->core.flag, --skipper[1]);
            continue;
        }
        if (b->core.l_qseq < minLength) {
            if (skipper[2])
                fprintf(stderr, "skipping: %s too short, this msg is printed %d times more \n",
                        bam_get_qname(b), --skipper[2]);
            continue;
        }
        uint8_t *nm = bam_aux_get(b, "NM");
        if (nm != NULL) {
            const int thiseditdist = (int)bam_aux2i(nm);
            if (editMin != -1 && thiseditdist < editMin)
                continue;
            if (editMax != -1 && thiseditdist > editMax)
                continue;

            if (minAni != -1.0 || maxAni != -1.0) {
                if (b->core.l_qseq <= 0)
                    continue;
                const double ani = 1.0 - (((double)thiseditdist) / ((double)b->core.l_qseq));
                if (minAni != -1.0 && ani < minAni)
                    continue;
                if (maxAni != -1.0 && ani > maxAni)
                    continue;
            }
        }
	dmg->damage_analysis(b, runmode != 0 ? b->core.tid : 0, 1);

        float mygc = gccontent(b);
        float mylen = b->core.l_qseq;

        int whichref = 0;
        if (runmode == 1)
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
            if (it->second.size() < 1000000)
                it->second.push_back(mygc);
            it = seqlens.find(whichref);
            if (it == seqlens.end()) {
                fprintf(stderr, "\t-> Error: iterator reached end in seqlens, will exit\n");
                return 1;
            }
            if (it->second.size() < 1000000)
                it->second.push_back(mylen);
        }
    }

    if (ret < -1) {
        fprintf(stderr, "\t-> Error: sam_read1 failed for input file: %s (ret=%d)\n", fname, ret);
        return 1;
    }

    return 0;
}

static void write_getdamage_outputs(const getdamage_args &args,
                                    bam_hdr_t *hdr,
                                    damage *dmg) {
    dmg->printit(stdout, args.printLength);
    (void)hdr;
    dmg->bwrite(args.onam, args.rlens_flat_out);
}

static getdamage_args init_getdamage_args() {
    getdamage_args args;
    args.minLength = 35;
    args.printLength = 5;
    args.editdistMin = -1;
    args.editdistMax = -1;
    args.minAni = -1.0;
    args.maxAni = -1.0;
    args.refName = NULL;
    args.fname = NULL;
    args.runmode = 0;
    args.onam = strdup("meta");
    args.nthreads = 4;
    args.ignore_errors = 0;
    args.rlens_flat_out = 0;
    return args;
}

static int parse_getdamage_args(int argc, char **argv, getdamage_args &args) {
    static struct option lopts[] = {
        {"threads", required_argument, 0, 'n'},
        {"rlens_flat_out", required_argument, 0, 'z'},
        {"fasta", required_argument, 0, 'f'},
        {"min_length", required_argument, 0, 'l'},
        {"edit_dist_min", required_argument, 0, 1000},
        {"edit_dist_max", required_argument, 0, 1001},
        {"min_ani", required_argument, 0, 1002},
        {"max_ani", required_argument, 0, 1003},
        {"minAni", required_argument, 0, 1002},
        {"maxAni", required_argument, 0, 1003},
        {"print_length", required_argument, 0, 'p'},
        {"run_mode", required_argument, 0, 'r'},
        {"ignore_errors", no_argument, 0, 'i'},
        {"out_prefix", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}};

    optind = 1;

    int c;
    while ((c = getopt_long(argc, argv,
                            "iz:n:f:l:p:r:o:h",
                            lopts, NULL)) >= 0) {
      switch (c) {
      case 'i':
        args.ignore_errors = 1;
        break;
      case 'n':
        args.nthreads = atoi(optarg);
        break;
      case 'f':
        free(args.refName);
        args.refName = strdup(optarg);
        break;
      case 'z':
        args.rlens_flat_out = atoi(optarg);
        break;
      case 'l':
        args.minLength = atoi(optarg);
        break;
      case 'p':
        args.printLength = atoi(optarg);
        break;
      case 1000:
        args.editdistMin = atoi(optarg);
        break;
      case 1001:
        args.editdistMax = atoi(optarg);
        break;
      case 1002:
        args.minAni = atof(optarg);
        break;
      case 1003:
        args.maxAni = atof(optarg);
        break;
      case 'r':
        args.runmode = atoi(optarg);
        break;
      case 'o':
        free(args.onam);
        args.onam = strdup(optarg);
        break;
      case 'h':
        usage_getdamage(stdout);
        return 2;
      default:
        fprintf(stderr, "Never here: %s\n", optarg);
        return 1;
      }
    }

    if ((argc - optind) > 1) {
        fprintf(stderr, "\t-> Error: expected exactly one input alignment file\n");
        return 1;
    }

    if (optind < argc) {
        free(args.fname);
        args.fname = strdup(argv[optind]);
    }

    return 0;
}

static int open_getdamage_input(const getdamage_args &args,
                                htsFormat *dingding2,
                                htsFile *&fp,
                                bam1_t *&b,
                                bam_hdr_t *&hdr) {
    if (args.refName) {
      int lenlen = 10 + strlen(args.refName) + 1;
      char *ref = (char *)malloc(lenlen);
      snprintf(ref, lenlen, "reference=%s", args.refName);
      hts_opt_add((hts_opt **)&dingding2->specific, ref);
      free(ref);
    }

    fp = sam_open_format(args.fname, "r", dingding2);
    if (fp == NULL) {
        fprintf(stderr, "[%s] nonexistant file: %s\n", __FUNCTION__, args.fname);
        return 1;
    }

    b = bam_init1();
    if (b == NULL) {
      fprintf(stderr, "\t-> Error: failed to allocate BAM record\n");
      return 1;
    }
    hdr = sam_hdr_read(fp);
    if (hdr == NULL) {
      fprintf(stderr, "\t-> Hello Doctor! im afraid there is an error reading the header\n");
      return 1;
    }

    int checkIfSorted(char *str);
    if (checkIfSorted(hdr->text)) {
      fprintf(stderr, "Input alignment file is not sorted.");
      if (!args.ignore_errors)
        return 1;
    }

    return 0;
}

int main_getdamage(int argc, char **argv) {
    if (argc == 1)
        return usage_getdamage(stderr);

    int rc = 0;
    getdamage_args args = init_getdamage_args();
    htsFile *fp = NULL;
    bam1_t *b = NULL;
    bam_hdr_t *hdr = NULL;
    damage *dmg = NULL;
    std::map<int, std::vector<float> > gcconts;
    std::map<int, std::vector<float> > seqlens;
    htsFormat *dingding2 = (htsFormat *)calloc(1, sizeof(htsFormat));

    rc = parse_getdamage_args(argc, argv, args);
    if (rc == 2) {
        rc = 0;
        goto cleanup;
    }
    if (rc != 0)
        goto cleanup;

    fprintf(stderr, "\t-> ./metaDMG-cpp refName: %s min_length: %d edit_dist_min: %d edit_dist_max: %d minAni: %f maxAni: %f print_length: %d run_mode: %d out_prefix: %s nthreads: %d ignore_errors: %d rlens_flat_out: %d\n", args.refName ? args.refName : "NULL", args.minLength, args.editdistMin, args.editdistMax, args.minAni, args.maxAni, args.printLength, args.runmode, args.onam, args.nthreads, args.ignore_errors, args.rlens_flat_out);
    if (args.fname == NULL) {
        rc = usage_getdamage(stderr);
        goto cleanup;
    }
    if (args.nthreads <= 0) {
        fprintf(stderr, "\t-> Error: threads must be greater than 0\n");
        rc = 1;
        goto cleanup;
    }
    if (args.printLength <= 0) {
        fprintf(stderr, "\t-> Error: print_length must be greater than 0\n");
        rc = 1;
        goto cleanup;
    }
    if (args.minLength < 0) {
        fprintf(stderr, "\t-> Error: min_length must be 0 or greater\n");
        rc = 1;
        goto cleanup;
    }
    if (args.editdistMin < -1) {
        fprintf(stderr, "\t-> Error: edit_dist_min must be -1 or greater\n");
        rc = 1;
        goto cleanup;
    }
    if (args.editdistMax < -1) {
        fprintf(stderr, "\t-> Error: edit_dist_max must be -1 or greater\n");
        rc = 1;
        goto cleanup;
    }
    if (args.editdistMin != -1 && args.editdistMax != -1 && args.editdistMin > args.editdistMax) {
        fprintf(stderr, "\t-> Error: edit_dist_min cannot be greater than edit_dist_max\n");
        rc = 1;
        goto cleanup;
    }
    if (args.minAni < -1.0 || args.minAni > 1.0) {
        fprintf(stderr, "\t-> Error: minAni/min_ani must be -1 (disabled) or between 0 and 1\n");
        rc = 1;
        goto cleanup;
    }
    if (args.maxAni < -1.0 || args.maxAni > 1.0) {
        fprintf(stderr, "\t-> Error: maxAni/max_ani must be -1 (disabled) or between 0 and 1\n");
        rc = 1;
        goto cleanup;
    }
    if (args.minAni != -1.0 && args.maxAni != -1.0 && args.minAni > args.maxAni) {
        fprintf(stderr, "\t-> Error: minAni/min_ani cannot be greater than maxAni/max_ani\n");
        rc = 1;
        goto cleanup;
    }
    if (args.runmode != 0 && args.runmode != 1) {
        fprintf(stderr, "\t-> Error: run_mode must be 0 or 1\n");
        rc = 1;
        goto cleanup;
    }

    rc = open_getdamage_input(args, dingding2, fp, b, hdr);
    if (rc != 0)
      goto cleanup;

    dmg = new damage(args.printLength, args.nthreads, 0);
    rc = process_getdamage_reads(fp, hdr, b, dmg, args.minLength, args.editdistMin, args.editdistMax, args.minAni, args.maxAni, args.runmode, gcconts, seqlens, args.fname);
    if (rc != 0)
      goto cleanup;

    write_getdamage_outputs(args, hdr, dmg);
    rc = write_getdamage_stats(args.onam, args.runmode, hdr, dmg, gcconts, seqlens);
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
    if (args.fname)
      free(args.fname);
    if (args.onam)
      free(args.onam);
    if (args.refName)
      free(args.refName);
    if (dingding2)
      free(dingding2);
    return rc;
}
