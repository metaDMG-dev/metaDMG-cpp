// gpl thorfinn@binf.ku.dk
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <strings.h>
#include <zlib.h>
#include <climits>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>

#include "main_print.h"
#include "ngsLCA.h"
#include "shared.h"
#include "version.h"

typedef std::map<int, char *> int2char;

double *getval(std::map<int, double *> &retmap, int2intvec &child, int taxid, int howmany) {
    // fprintf(stderr,"getval\t%d\t%d\n",taxid,howmany);
    std::map<int, double *>::iterator it = retmap.find(taxid);
    if (it != retmap.end()) {
        // fprintf(stderr,"has found: %d in retmap\n",it->first);
#if 0
    fprintf(stderr,"val\t%d",taxid);
    for(int i=0;i<3*howmany+1;i++)
      fprintf(stderr,"\t%f",it->second[i]);
    fprintf(stderr,"\n");
#endif
        return it->second;
    }
    double *ret = new double[3 * howmany + 1];
    for (int i = 0; i < 3 * howmany + 1; i++)
        ret[i] = 0.0;
    if (child.size() > 0) {  // if we have supplied -nodes
        int2intvec::iterator it2 = child.find(taxid);
        if (it2 != child.end()) {
            std::vector<int> &avec = it2->second;
            for (size_t i = 0; i < avec.size(); i++) {
                //	fprintf(stderr,"%d/%d %d\n",i,avec.size(),avec[i]);
                double *tmp = getval(retmap, child, avec[i], howmany);
                for (int i = 0; i < 3 * howmany + 1; i++)
                    ret[i] += tmp[i];
            }
        }
    }

    retmap[taxid] = ret;

    return ret;
}
int INT_WARN = 1;


//apparantly there is an issue when data is not only leaf.
std::map<int,mydataD> getval_full_norec(std::map<int, mydataD> &retmap, int2int &parent, int howmany) {
    // fprintf(stderr,"getval\t%d\t%d\n",taxid,howmany);

  std::map<int, mydataD> results;
  
  //funky modern syntax below
  //loop over all entries in retmap. lizard king 
  //  for (const auto &[taxid, data] : retmap) {

  for(  std::map<int, mydataD>::iterator it=retmap.begin();it!=retmap.end();it++){
    int taxid = it->first;
    mydataD data = it->second;
    
    int current = taxid;
    //   fprintf(stderr,"taxid: %d\n",taxid);
     while (true) {
       //   fprintf(stderr,"current: %d\n",current);
       mydataD &target = results[current]; //<- will call constructor if doesnt exists. magic magic
       if (!target.fwD) {//if it doesnt exists then allocate. Just like fw.D!=NULL
	 target.fwD = new double[16 * howmany]();
	 target.bwD = new double[16 * howmany]();
	 target.nal = 0;
       }
       
       target.nal += data.nal;
       
       if (target.nal > INT_MAX && INT_WARN) {
	 fprintf(stderr, "\t-> Potential issue, sum of alignment counts exceeds INT_MAX\n");
	 INT_WARN = 0;
       }
       
       for (int i = 0; i < 16 * howmany; ++i) {
	 target.fwD[i] += data.fwD[i];
	 target.bwD[i] += data.bwD[i];
       }
       
       auto it = parent.find(current);
       //       fprintf(stderr,"second: %d\n",it->second);
       if (it == parent.end()) break;
       //break if up is same as current. That only happens with root that has tqxid=1
       if(current == it->second)
	 break;
       current = it->second;
     }
     
  }

  return results;
}

int mywarn =1;

mydata2 getval_stats(std::map<int, mydata2> &retmap, int2intvec &child, int taxid) {
  if(mywarn>0){
    fprintf(stderr,"PAS PAA SATAN\n");
    mywarn--;
  }
    // fprintf(stderr,"getval\t%d\t%d\n",taxid,howmany);
    std::map<int, mydata2>::iterator it = retmap.find(taxid);
    if (it != retmap.end()) {
        // fprintf(stderr,"has found: %d in retmap\n",it->first);
#if 0
    fprintf(stderr,"val\t%d",taxid);
    for(int i=0;i<3*howmany+1;i++)
      fprintf(stderr,"\t%f",it->second[i]);
    fprintf(stderr,"\n");
#endif
        return it->second;
    }
    mydata2 ret;
    ret.nreads = 0;
    ret.data = new double[4];
    ret.data[0] = ret.data[1] = ret.data[2] = ret.data[3] = 0;

    if (child.size() > 0) {  // if we have supplied -nodes
        int2intvec::iterator it2 = child.find(taxid);
        if (it2 != child.end()) {
            std::vector<int> &avec = it2->second;
            for (size_t i = 0; i < avec.size(); i++) {
                mydata2 tmp = getval_stats(retmap, child, avec[i]);
                ret.nreads += tmp.nreads;
            }

            for (size_t i = 0; i < avec.size(); i++) {
                //	fprintf(stderr,"%d/%d %d\n",i,avec.size(),avec[i]);
                mydata2 tmp = getval_stats(retmap, child, avec[i]);
                if (tmp.nreads == 0)
                    continue;
                double a = tmp.nreads;
                double b = ret.nreads;
                double scal = a / b;
                //	fprintf(stderr,"a: %f b: %f scal: %f\n",a,b,scal);
                for (int i = 0; i < 4; i++)
                    ret.data[i] += tmp.data[i] * scal;
            }
        }
    }

    retmap[taxid] = ret;

    return ret;
}

static int usage_print(FILE *fp) {
    fprintf(fp, "Usage: ./metaDMG-cpp print file.bdamage.gz [options]\n");
    fprintf(fp, "Options: -names FILE -bam FILE -nodes FILE -r INT -howmany INT -ctga -countout -doOld\n");
    return 0;
}

static int usage_print2(FILE *fp) {
    fprintf(fp, "Usage: ./metaDMG-cpp print2 file.bdamage.gz [options]\n");
    fprintf(fp, "Options: -acc2tax FILE -bam FILE -nodes FILE -r INT -howmany INT -ctga -countout -doOld\n");
    return 0;
}

static int usage_print_all(FILE *fp) {
    fprintf(fp, "Usage: ./metaDMG-cpp print_all file.bdamage.gz -names FILE [-nodes FILE]\n");
    return 0;
}

static int usage_print_ugly(FILE *fp) {
    fprintf(fp, "Usage: ./metaDMG-cpp print_ugly file.bdamage.gz --names FILE [options]\n");
    fprintf(fp, "Options: --nodes FILE --lcastat FILE --bam FILE --out_prefix STR\n");
    return 0;
}

int main_print(int argc, char **argv) {
  //  fprintf(stderr,"\nYOYOYOYOYYOOYOOY\n");
    if (argc == 1 || (argc == 2 && (!strcasecmp(argv[1], "-h") || !strcasecmp(argv[1], "--help"))))
        return usage_print(stderr);
    char *infile = NULL;
    char *inbam = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;
    int ctga = 0;  // only print ctga errors
    int search = -1;
    int countout = 0;

    int howmany = 15;
    int doold = 1;
    while (*(++argv)) {
        if (!strcasecmp("-h", *argv) || !strcasecmp("--help", *argv))
            return usage_print(stderr);
        else if (strcasecmp("-names", *argv) == 0)
            infile_names = strdup(*(++argv));
        else if (strcasecmp("-bam", *argv) == 0)
            inbam = strdup(*(++argv));
        else if (strcasecmp("-r", *argv) == 0)
            search = atoi(*(++argv));
        else if (strcasecmp("-howmany", *argv) == 0)
            howmany = atoi(*(++argv));
        else if (strcasecmp("-ctga", *argv) == 0)
            ctga = 1;
        else if (strcasecmp("-doOld", *argv) == 0)
            doold = 1;
        else if (strcasecmp("-countout", *argv) == 0)
            countout = 1;
        else if (strcasecmp("-nodes", *argv) == 0)
            infile_nodes = strdup(*(++argv));
        else
            infile = strdup(*argv);
    }

    fprintf(stderr,
	    "infile: %s inbam: %s search: %d ctga: %d countout: %d nodes: %s names: %s howmany: %d\n",
	    infile ? infile : "NULL",
	    inbam ? inbam : "NULL",
	    search,
	    ctga,
	    countout,
	    infile_nodes ? infile_nodes : "NULL",
	    infile_names ? infile_names : "NULL",
	    howmany
	    );
    if (!infile) {
      fprintf(stderr, "\t-> Error: infile is NULL or could not be opened, will exit\n");
      exit(1);
    }

    int2char name_map;
    if (infile_names != NULL)
        name_map = parse_names(infile_names);

    // map of taxid -> taxid
    int2int parent;
    // map of taxid -> rank
    int2char rank;
    // map of parent -> child taxids
    int2intvec child;

    if (infile_nodes != NULL)
        parse_nodes(infile_nodes, rank, parent, child, 1);
    if (search != -1 && doold == 0) {
        std::map<int, double *> retmap = load_bdamage3(infile, howmany);
        double *dbl = getval(retmap, child, search, howmany);
        double *dbldbl = new double[3 * howmany + 1];  // 3 because ct,ga,other
        dbldbl[0] = dbl[0];
        for (int i = 0; i < 3 * howmany; i++)
            dbldbl[i + 1] = dbl[1 + i] / dbl[0];

        fprintf(stdout, "%d\t%.0f", search, dbldbl[0]);
        for (int i = 0; i < 3 * howmany; i++)
            fprintf(stdout, "\t%f", dbldbl[1 + i]);
        fprintf(stdout, "\n");
	delete [] dbldbl;
        return 0;
    }

    BGZF *bgfp = NULL;
    samFile *samfp = NULL;
    bam_hdr_t *hdr = NULL;

    if (((bgfp = bgzf_open(infile, "r"))) == NULL) {
        fprintf(stderr, "Could not open input bdamage.gz file: %s\n", infile);
        return 1;
    }

    if (inbam != NULL) {
        if (((samfp = sam_open_format(inbam, "r", NULL))) == NULL) {
            fprintf(stderr, "Could not open input BAM file: %s\n", inbam);
            return 1;
        }
        if (((hdr = sam_hdr_read(samfp))) == NULL) {
            fprintf(stderr, "Could not read header for: %s\n", inbam);
            return 1;
        }
    }

    int printlength;
    if (bgzf_read(bgfp, &printlength, sizeof(int)) != sizeof(int)) {
      fprintf(stderr, "\t-> Error: failed to read expected number of bytes from bgzf file, will exit\n");
      exit(1);
    }
    fprintf(stderr, "\t-> printlength(howmany) from inputfile: %d\n", printlength);

    int ref_nreads[2];

    if (ctga == 0) {
        if (hdr != NULL)
            fprintf(stdout, "Reference\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
        else if (infile_names != NULL)
            fprintf(stdout, "FunkyName\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
        else
            fprintf(stdout, "taxid\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
    } else {
    }

    float data[16];
    while (1) {
        int nread = bgzf_read(bgfp, ref_nreads, 2 * sizeof(int));
	double *ctgas = (double *)malloc(2 * printlength * sizeof(double));
	if (ctgas == NULL) {
	  fprintf(stderr, "\t-> Error: failed to allocate memory for ctgas, will exit\n");
	  exit(1);
	}
        if (nread == 0)
            break;
	if (nread != 2 * sizeof(int)) {
	  fprintf(stderr, "\t-> Error: unexpected number of bytes read (nread != 2*sizeof(int)), will exit\n");
	  exit(1);
	}
        for (int at = 0; at < printlength; at++) {
	  if (bgzf_read(bgfp, data, sizeof(float) * 16) != 16 * sizeof(float)) {
	    fprintf(stderr, "\t-> Error: failed to read expected number of bytes (16 floats) from bgzf file, will exit\n");
	    exit(1);
	  }
	    
            if (at >= howmany)
                continue;
            if (search == -1 || search == ref_nreads[0]) {
                if (ctga == 0) {
                    if (hdr != NULL)
                        fprintf(stdout, "%s\t%d\t5\'\t%d", hdr->target_name[ref_nreads[0]], ref_nreads[1], at);
                    else if (infile_names != NULL) {
                        int2char::iterator itt = name_map.find(ref_nreads[0]);
                        if (itt == name_map.end()) {
                            fprintf(stderr, "\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n", ref_nreads[0], infile_names);
                            exit(0);
                        }
                        fprintf(stdout, "\"%s\"\t%d\t5\'\t%d", itt->second, ref_nreads[1], at);
                    } else
                        fprintf(stdout, "%d\t%d\t5\'\t%d", ref_nreads[0], ref_nreads[1], at);
                } else {
                    if (at == 0)
                        fprintf(stdout, "%d\t%d", ref_nreads[0], ref_nreads[1]);
                }
                if (countout == 1) {
                    for (int i = 0; i < 16; i++)
                        fprintf(stdout, "\t%f", data[i]);
                    fprintf(stdout, "\n");
                } else {
                    float flt[16];

                    for (int i = 0; i < 4; i++) {
                        double tsum = 0;
                        for (int j = 0; j < 4; j++) {
                            tsum += data[i * 4 + j];
                            flt[i * 4 + j] = data[i * 4 + j];
                        }
                        if (tsum == 0) tsum = 1;
                        for (int j = 0; j < 4; j++)
                            flt[i * 4 + j] /= tsum;
                    }

                    if (ctga == 0) {
                        for (int j = 0; j < 16; j++)
                            fprintf(stdout, "\t%f", flt[j]);
                        fprintf(stdout, "\n");
                    } else
                        ctgas[at] = flt[7];
                }
            }
        }
        if (search == -1 || search == ref_nreads[0]) {
            if (ctga == 1) {
                for (int i = 0; i < printlength; i++)
                    fprintf(stdout, "\t%f", ctgas[i]);
            }
        }

        for (int at = 0; at < printlength; at++) {
	  if (bgzf_read(bgfp, data, sizeof(float) * 16) != 16 * sizeof(float)) {
	    fprintf(stderr, "\t-> Error: failed to read expected number of bytes (16 floats) from bgzf file, will exit\n");
	    exit(1);
	  }
            if (at >= howmany)
                continue;
            if (search == -1 || search == ref_nreads[0]) {
                if (ctga == 0) {
                    if (hdr != NULL)
                        fprintf(stdout, "%s\t%d\t3\'\t%d", hdr->target_name[ref_nreads[0]], ref_nreads[1], at);
                    else if (infile_names != NULL) {
                        int2char::iterator itt = name_map.find(ref_nreads[0]);
                        if (itt == name_map.end()) {
                            fprintf(stderr, "\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n", ref_nreads[0], infile_names);
                            exit(0);
                        }
                        fprintf(stdout, "\"%s\"\t%d\t3\'\t%d", itt->second, ref_nreads[1], at);
                    } else
                        fprintf(stdout, "%d\t%d\t3\'\t%d", ref_nreads[0], ref_nreads[1], at);
                }
                if (countout == 1) {
                    for (int i = 0; i < 16; i++)
                        fprintf(stdout, "\t%f",data[i]);
                    fprintf(stdout, "\n");
                } else {
                    float flt[16];

                    for (int i = 0; i < 4; i++) {
                        double tsum = 0;
                        for (int j = 0; j < 4; j++) {
                            tsum += data[i * 4 + j];
                            flt[i * 4 + j] = data[i * 4 + j];
                        }
                        if (tsum == 0) tsum = 1;
                        for (int j = 0; j < 4; j++)
                            flt[i * 4 + j] /= tsum;
                    }
                    if (ctga == 0) {
                        for (int j = 0; j < 16; j++)
                            fprintf(stdout, "\t%f", flt[j]);
                        fprintf(stdout, "\n");
                    } else
                        ctgas[at + printlength] = flt[8];
                }
            }
        }
        if (search == -1 || search == ref_nreads[0]) {
            if (ctga == 1) {
                for (int i = 0; i < printlength; i++)
                    fprintf(stdout, "\t%f", ctgas[printlength + i]);
                fprintf(stdout, "\n");
            }
        }
	free(ctgas);
    }
    //cleanup
    for(int2char::iterator it=name_map.begin();it!=name_map.end();it++)
      free(it->second);
    for(int2char::iterator it=rank.begin();it!=rank.end();it++)
      free(it->second);

    if (bgfp)
        bgzf_close(bgfp);
    if (hdr)
        bam_hdr_destroy(hdr);
    if (samfp)
        sam_close(samfp);
    if(infile)
      free(infile);
    if(infile_nodes)
      free(infile_nodes);
    if(infile_names)
      free(infile_names);
    if(inbam)
      free(inbam);
    return 0;
}

int main_print2(int argc, char **argv) {
    if (argc == 1 || (argc == 2 && (!strcasecmp(argv[1], "-h") || !strcasecmp(argv[1], "--help"))))
        return usage_print2(stderr);
    char *infile = NULL;
    char *inbam = NULL;
    char *acc2tax = NULL;
    int ctga = 0;  // only print ctga errors
    int search = -1;
    int countout = 0;
    char *infile_nodes = NULL;
    int howmany = 15;
    int doold = 0;
    while (*(++argv)) {
        if (!strcasecmp("-h", *argv) || !strcasecmp("--help", *argv))
            return usage_print2(stderr);
        else if (strcasecmp("-acc2tax", *argv) == 0)
            acc2tax = strdup(*(++argv));
        else if (strcasecmp("-bam", *argv) == 0)
            inbam = strdup(*(++argv));
        else if (strcasecmp("-r", *argv) == 0)
            search = atoi(*(++argv));
        else if (strcasecmp("-howmany", *argv) == 0)
            howmany = atoi(*(++argv));
        else if (strcasecmp("-ctga", *argv) == 0)
            ctga = 1;
        else if (strcasecmp("-doOld", *argv) == 0)
            doold = 1;
        else if (strcasecmp("-countout", *argv) == 0)
            countout = 1;
        else if (strcasecmp("-nodes", *argv) == 0)
            infile_nodes = strdup(*(++argv));
        else
            infile = strdup(*argv);
    }

    fprintf(stderr,
	    "infile: %s inbam: %s names: %s search: %d ctga: %d countout: %d nodes: %s\n",
	    infile ? infile : "NULL",
	    inbam ? inbam : "NULL",
	    acc2tax ? acc2tax : "NULL",
	    search,
	    ctga,
	    countout,
	    infile_nodes ? infile_nodes : "NULL"
	    );
    if (!infile) {
      fprintf(stderr, "\t-> Error: infile is NULL or could not be opened, will exit\n");
      exit(1);
    }
    int2char name_map;
    if (acc2tax != NULL)
        name_map = parse_names(acc2tax);

    // map of taxid -> taxid
    int2int parent;
    // map of taxid -> rank
    int2char rank;
    // map of parent -> child taxids
    int2intvec child;

    if (infile_nodes != NULL)
        parse_nodes(infile_nodes, rank, parent, child, 1);
    if (search != -1 && doold == 0) {
        std::map<int, double *> retmap = load_bdamage3(infile, howmany);
        double *dbl = getval(retmap, child, search, howmany);
        double *dbldbl = new double [3 * howmany + 1];  // 3 because ct,ga,other
        dbldbl[0] = dbl[0];
        for (int i = 0; i < 3 * howmany; i++)
            dbldbl[i + 1] = dbl[1 + i] / dbl[0];

        fprintf(stdout, "%d\t%.0f", search, dbldbl[0]);
        for (int i = 0; i < 3 * howmany; i++)
            fprintf(stdout, "\t%f", dbldbl[1 + i]);
        fprintf(stdout, "\n");
	delete [] dbldbl;
        return 0;
    }

    BGZF *bgfp = NULL;
    samFile *samfp = NULL;
    bam_hdr_t *hdr = NULL;

    if (((bgfp = bgzf_open(infile, "r"))) == NULL) {
        fprintf(stderr, "Could not open input bdamage.gz file: %s\n", infile);
        return 1;
    }

    if (inbam != NULL) {
        if (((samfp = sam_open_format(inbam, "r", NULL))) == NULL) {
            fprintf(stderr, "Could not open input BAM file: %s\n", inbam);
            return 1;
        }
        if (((hdr = sam_hdr_read(samfp))) == NULL) {
            fprintf(stderr, "Could not read header for: %s\n", inbam);
            return 1;
        }
    }

    int printlength;
    if (bgzf_read(bgfp, &printlength, sizeof(int)) != sizeof(int)) {
      fprintf(stderr, "\t-> Error: failed to read expected number of bytes for printlength, will exit\n");
      exit(1);
    }
    fprintf(stderr, "\t-> printlength(howmany) from inputfile: %d\n", printlength);

    int ref_nreads[2];
    char *type_name = NULL;
    if (hdr != NULL)
        type_name = strdup("Reference");
    else if (acc2tax != NULL)
        type_name = strdup("FunkyName");
    else
        type_name = strdup("taxid");

    if (ctga == 0) {
        fprintf(stdout, "%s\tNalignments\tDirection\tPos\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n", type_name);
    } else {
        fprintf(stdout, "%s\tNalignment", type_name);
        for (int i = 0; i < howmany; i++)
            fprintf(stdout, "\tCT_%d", i);
        for (int i = 0; i < howmany; i++)
            fprintf(stdout, "\tGA_%d", i);
        fprintf(stdout, "\n");
    }

    float data[16];

    while (1) {
        int nread = bgzf_read(bgfp, ref_nreads, 2 * sizeof(int));
	double *ctgas = (double *)malloc(2 * printlength * sizeof(double));
	if (ctgas == NULL) {
	  fprintf(stderr, "\t-> Error: failed to allocate memory for ctgas, will exit\n");
	  exit(1);
	}
        
        if (nread == 0)
            break;
        fprintf(stderr, "ref: %d nreads: %d\n", ref_nreads[0], ref_nreads[1]);
	if (nread != 2 * sizeof(int)) {
	  fprintf(stderr, "\t-> Error: unexpected number of bytes read (nread != 2*sizeof(int)), will exit\n");
	  exit(1);
	}
	for (int at = 0; at < printlength; at++) {
	  if (bgzf_read(bgfp, data, sizeof(float) * 16) != 16 * sizeof(float)) {
	    fprintf(stderr, "\t-> Error: failed to read expected number of bytes (16 floats) from bgzf file, will exit\n");
	    exit(1);
	  }
            if ((at + 1) > howmany)
                continue;
            if (search == -1 || search == ref_nreads[0]) {
                if (ctga == 0) {
                    if (hdr != NULL)
                        fprintf(stdout, "%s\t%d\t5\'\t%d", hdr->target_name[ref_nreads[0]], ref_nreads[1], at);
                    else if (acc2tax != NULL) {
                        int2char::iterator itt = name_map.find(ref_nreads[0]);
                        if (itt == name_map.end()) {
                            fprintf(stderr, "\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n", ref_nreads[0], acc2tax);
                            exit(0);
                        }
                        fprintf(stdout, "\"%s\"\t%d\t5\'\t%d", itt->second, ref_nreads[1], at);
                    } else
                        fprintf(stdout, "%d\t%d\t5\'\t%d", ref_nreads[0], ref_nreads[1], at);
                } else {
                    if (at == 0)
                        fprintf(stdout, "%d\t%d", ref_nreads[0], ref_nreads[1]);
                }
                if (countout == 1) {
                    for (int i = 0; i < 16; i++)
                        fprintf(stdout, "\t%f", data[i]);
                    fprintf(stdout, "\n");
                } else {
                    float flt[16];

                    for (int i = 0; i < 4; i++) {
                        double tsum = 0;
                        for (int j = 0; j < 4; j++) {
                            tsum += data[i * 4 + j];
                            flt[i * 4 + j] = data[i * 4 + j];
                        }
                        if (tsum == 0) tsum = 1;
                        for (int j = 0; j < 4; j++)
                            flt[i * 4 + j] /= tsum;
                    }

                    if (ctga == 0) {
                        for (int j = 0; j < 16; j++)
                            fprintf(stdout, "\t%f", flt[j]);
                        fprintf(stdout, "\n");
                    } else
                        ctgas[at] = flt[7];
                }
            }
        }
        if (search == -1 || search == ref_nreads[0]) {
            if (ctga == 1) {
                for (int i = 0; i < howmany; i++)
                    fprintf(stdout, "\t%f", ctgas[i]);
            }
        }

        for (int at = 0; at < printlength; at++) {
	  if (bgzf_read(bgfp, data, sizeof(float) * 16) != 16 * sizeof(float)) {
	    fprintf(stderr, "\t-> Error: failed to read expected number of bytes (16 floats) from bgzf file, will exit\n");
	    exit(1);
	  }
            if (at + 1 > howmany)
                continue;
            if (search == -1 || search == ref_nreads[0]) {
                if (ctga == 0) {
                    if (hdr != NULL)
                        fprintf(stdout, "%s\t%d\t3\'\t%d", hdr->target_name[ref_nreads[0]], ref_nreads[1], at);
                    else if (acc2tax != NULL) {
                        int2char::iterator itt = name_map.find(ref_nreads[0]);
                        if (itt == name_map.end()) {
                            fprintf(stderr, "\t-> Problem finding taxid: \'%d' in namedatabase: \'%s\'\n", ref_nreads[0], acc2tax);
                            exit(0);
                        }
                        fprintf(stdout, "\"%s\"\t%d\t3\'\t%d", itt->second, ref_nreads[1], at);
                    } else
                        fprintf(stdout, "%d\t%d\t3\'\t%d", ref_nreads[0], ref_nreads[1], at);
                }
                if (countout == 1) {
                    for (int i = 0; i < 16; i++)
                        fprintf(stdout, "\t%f", data[i]);
                    fprintf(stdout, "\n");
                } else {
                    float flt[16];

                    for (int i = 0; i < 4; i++) {
                        double tsum = 0;
                        for (int j = 0; j < 4; j++) {
                            tsum += data[i * 4 + j];
                            flt[i * 4 + j] = data[i * 4 + j];
                        }
                        if (tsum == 0) tsum = 1;
                        for (int j = 0; j < 4; j++)
                            flt[i * 4 + j] /= tsum;
                    }
                    if (ctga == 0) {
                        for (int j = 0; j < 16; j++)
                            fprintf(stdout, "\t%f", flt[j]);
                        fprintf(stdout, "\n");
                    } else
                        ctgas[at + printlength] = flt[8];
                }
            }
        }

        if (search == -1 || search == ref_nreads[0]) {
            if (ctga == 1) {
                for (int i = 0; i < howmany; i++)
                    fprintf(stdout, "\t%f", ctgas[printlength + i]);
                fprintf(stdout, "\n");
            }
        }
	free(ctgas);
    }
    //clean up
    for(int2char::iterator it=name_map.begin();it!=name_map.end();it++)
      free(it->second);
    for(int2char::iterator it=rank.begin();it!=rank.end();it++)
      free(it->second);

    if (bgfp)
        bgzf_close(bgfp);
    if (hdr)
        bam_hdr_destroy(hdr);
    if (samfp)
        sam_close(samfp);
    if(type_name)
      free(type_name);
    if(infile_nodes)
      free(infile_nodes);
    if(infile)
      free(infile);
    if(inbam)
      free(inbam);
    if(acc2tax)
      free(acc2tax);
    return 0;
}

int2int getlcadist(char *fname) {
    //  fprintf(stderr,"fname: %s\n",fname);
    int2int lcadist;
    if (fname == NULL) {
      fprintf(stderr, "\t-> Error: fname is NULL, will exit\n");
      exit(1);
    }
    
    char tmp[1024];
    
    size_t len = strlen(fname);
    if (len + 5 >= sizeof(tmp)) {  // ".stat" + '\0' = 5
      fprintf(stderr,"\t-> Ridiculus long filename: \'%s\'\n",fname);//<- harry potter    
      exit(1);
    }
    snprintf(tmp, sizeof(tmp), "%s.stat", fname);

    FILE *fp = NULL;
    fp = fopen(tmp, "rb");
    if (fp == NULL) {
      fprintf(stderr, "\t-> Error: failed to open file %s, will exit\n", tmp);
      exit(1);
    }
    char buf[4096];
    while (fgets(buf, 4096, fp)) {
      char *tok1 = strtok(buf, "\t\n ");
      char *tok2 = strtok(NULL, "\t\n ");
      if(!tok1||!tok2)
	continue;
      int key = atoi(tok1);
      int val = atoi(tok2);
      lcadist[key] = val;
    }
    fprintf(stderr, "\t-> Done reading: %lu entries from file: \'%s\'\n", lcadist.size(), tmp);
    fclose(fp);
    return lcadist;
}

std::map<int, double *> getcsv(char *fname) {
    std::map<int, double *> ret;
    FILE *fp = NULL;
    fp = fopen(fname, "rb");
    if (fp == NULL) {
      fprintf(stderr, "\t-> Error: failed to open file %s, will exit\n", fname);
      exit(1);
    }
    char buf[4096];
    if (fgets(buf, 4096, fp) == NULL) {
      fprintf(stderr, "\t-> Error: failed to read header line from file, will exit\n");
      exit(1);
    }
    while (fgets(buf, 4096, fp)) {
        int key = atoi(strtok(buf, "\t\n, "));
        double *valval = new double[2];
        valval[0] = atof(strtok(NULL, "\t\n, "));
        valval[1] = atof(strtok(NULL, "\t\n, "));
        ret[key] = valval;
    }
    fprintf(stderr, "\t-> Done reading: %lu entries from file: \'%s\'\n", ret.size(), fname);
    return ret;
}

//this very verbose function is translated by a thinking machine
int main_print_all(int argc, char **argv) {
    if (argc <= 2 || (argc == 2 && (!strcasecmp(argv[1], "-h") || !strcasecmp(argv[1], "--help"))))
        return usage_print_all(stderr);

    char *infile_bdamage = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;

    while (*(++argv)) {
        if (!strcasecmp("-h", *argv) || !strcasecmp("--help", *argv)) {
            return usage_print_all(stderr);
        } else if (strcasecmp("-names", *argv) == 0) {
            if (*(argv + 1) == NULL) {
                fprintf(stderr, "\t-> Error: -names requires an argument, will exit\n");
                exit(1);
            }
            infile_names = strdup(*(++argv));
        } else if (strcasecmp("-nodes", *argv) == 0) {
            if (*(argv + 1) == NULL) {
                fprintf(stderr, "\t-> Error: -nodes requires an argument, will exit\n");
                exit(1);
            }
            infile_nodes = strdup(*(++argv));
        } else {
            infile_bdamage = strdup(*argv);
        }
    }

    if (infile_bdamage == NULL) {
        fprintf(stderr, "\t-> Error: no input bdamage file provided, will exit\n");
        exit(1);
    }
    if (infile_names == NULL) {
        fprintf(stderr, "\t-> Error: no names file provided, will exit\n");
        exit(1);
    }

    fprintf(stderr, "infile_names: %s infile_bdamage: %s nodes: %s ",
            infile_names ? infile_names : "NULL",
            infile_bdamage ? infile_bdamage : "NULL",
            infile_nodes ? infile_nodes : "NULL");
    fprintf(stderr, "#VERSION:%s\n", METADAMAGE_VERSION);

    int2int parent;
    int2char rank;
    int2intvec child;

    if (infile_nodes != NULL)
        parse_nodes(infile_nodes, rank, parent, child, 1);

    int howmany = 5;
    std::map<int, double *> retmap = load_bdamage3(infile_bdamage, howmany);
    fprintf(stderr, "\t-> number of entries in damage pattern file: %lu\n", retmap.size());

    int2char name = parse_names(infile_names);

    size_t presize = retmap.size();

    getval(retmap, child, 1, howmany);

    size_t postsize = retmap.size();
    fprintf(stderr, "\t-> pre: %lu post:%lu grownbyfactor: %f\n",
            (unsigned long)presize,
            (unsigned long)postsize,
            presize ? (double)postsize / presize : 0.0);

    for (std::map<int, double *>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        double *dbl = it->second;

        char *myrank = NULL;
        char *myname = NULL;

        int2char::iterator itc = rank.find(taxid);
        if (itc != rank.end())
            myrank = itc->second;

        itc = name.find(taxid);
        if (itc != name.end())
            myname = itc->second;

        if (dbl[0] > 0) {
            double dbldbl[3 * 5 + 1];
            dbldbl[0] = dbl[0];
            for (int i = 0; i < 3 * howmany; i++)
                dbldbl[i + 1] = dbl[1 + i] / dbl[0];

            fprintf(stdout, "%d:\"%s\":\"%s\":%.0f",
                    taxid,
                    myname ? myname : "NULL",
                    myrank ? myrank : "NULL",
                    dbldbl[0]);
            for (int i = 0; i < 3 * howmany; i++)
                fprintf(stdout, "\t%f", dbldbl[1 + i]);
            fprintf(stdout, "\n");
        }
    }

    free(infile_bdamage);
    free(infile_nodes);
    free(infile_names);

    return 0;
}
int main_print_ugly(int argc, char **argv) {
  htsFormat *dingding2 = (htsFormat *)calloc(1, sizeof(htsFormat));
  fprintf(stderr,"\t-> print_ugly functionality will be removed\n");
  
    if (argc <= 1 || (argc == 2 && (!strcasecmp(argv[1], "-h") || !strcasecmp(argv[1], "--help"))))
        return usage_print_ugly(stderr);
    char *infile_bdamage = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;
    char *infile_lcastat = NULL;
    char *infile_bam = NULL;
    char *out_prefix = NULL;
    int howmany;
    int abe = 4;//<- monkey

    while (*(++argv)) {
        if (!strcasecmp("-h", *argv) || !strcasecmp("--help", *argv))
            return usage_print_ugly(stderr);
        else if (strcasecmp("--names", *argv) == 0)
            infile_names = strdup(*(++argv));
        else if (strcasecmp("--nodes", *argv) == 0)
            infile_nodes = strdup(*(++argv));
        else if (strcasecmp("--lcastat", *argv) == 0)
            infile_lcastat = strdup(*(++argv));
        else if (strcasecmp("--bam", *argv) == 0)
            infile_bam = strdup(*(++argv));
        else if (strcasecmp("--out_prefix", *argv) == 0)
            out_prefix = strdup(*(++argv));
        else
            infile_bdamage = strdup(*argv);
    }
    // Use input file name as default prefix
    if (out_prefix == NULL) {
      out_prefix = strdup(infile_bdamage);
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
    snprintf(buf, 1024, "%s.uglyprint.mismatch.gz", out_prefix);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    gzFile fpfpfp = gzopen(buf, "wb");
    gzprintf(fpfpfp, "taxidStr\tdirection\tposition\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
    // map of taxid -> taxid
    int2int parent;
    // map of taxid -> rank
    int2char rank;
    // map of parent -> child taxids
    int2intvec child;
    kstring_t *kstr = new kstring_t;
    kstr->s = NULL; kstr->l =kstr->m = 0;

    if (infile_nodes != NULL)
        parse_nodes(infile_nodes, rank, parent, child, 1);

    std::map<int, mydataD> retmap = load_bdamage_full(infile_bdamage, howmany);
    fprintf(stderr, "\t-> Number of entries in damage pattern file: %lu printlength(howmany):%d\n", retmap.size(), howmany);

    int2char name_map;

    if (infile_names)
        name_map = parse_names(infile_names);

    float presize = retmap.size();
    //    getval_full(retmap, child,1, howmany);  // this will do everything
    std::map<int, mydataD>  results =   getval_full_norec(retmap, parent, howmany);  // this will do everything

    retmap = results;
    float postsize = retmap.size();
    fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);

    for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        //mydataD md = it->second;
        if (it->second.nal == 0)
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
    snprintf(buf, 1024, "%s.uglyprint.stat.gz", out_prefix);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    fpfpfp = gzopen(buf, "wb");
    gzprintf(fpfpfp, "taxid\tname\trank\tnalign\tnreads\tmean_rlen\tvar_rlen\tmean_gc\tvar_gc\tlca\ttaxa_path\n");
    std::map<int, mydata2> stats;
    if (infile_lcastat)
        stats = load_lcastat(infile_lcastat,1);
#if 1
    if(child.size()>0)
       getval_stats(stats, child, 1);  // this will do everything
#endif
    void aggr_stat3000(std::map<int, mydata2> &stats,int2int &parent);
    if(0&&parent.size()>0)
      aggr_stat3000(stats,parent);
    for (std::map<int, mydata2>::iterator it = stats.begin(); 1 && it != stats.end(); it++) {
        std::map<int, mydataD>::iterator itold = retmap.find(it->first);
        size_t nalign = 0;
        if (itold == retmap.end()) {
	  if(abe>0)
	    fprintf(stderr, "[%s]\t-> Problem finding taxid: %d will stop print this in: %d\n",__FUNCTION__, it->first,abe--);
	  
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
            gzprintf(fpfpfp, "%d\t\"%s\"\t\"%s\"\t%d\t%d\t%f\t%f\t%f\t%f", it->first, myname, myrank, nalign, it->second.nreads, it->second.data[0], it->second.data[1], it->second.data[2], it->second.data[3]);

	    print_chain(kstr, it->first, parent, rank, name_map,1);
	    gzwrite(fpfpfp,kstr->s,kstr->l);
	    kstr->l = 0;
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
    free(kstr->s);
    delete kstr;
    free(dingding2);
    return 0;
}
