int mod_in[] = {1649555, 1401172, 1582271, 374764, 242716, 1793725, 292451, 38298, 63403, 67357, 163247, 328623, 356150, 502130, 545877, 996989, 996990, 1086724, 1169024, 1576640, 1769757, 1802981, 1811974, 1955118};
int mod_out[] = {1333996, 1333996, 1582270, 1914213, 1917265, 1915309, 263865, 2801, 1916091, 285450, 1762941, 1916091, 157727, 1932322, 376133, 1762939, 1762946, 430531, 1169025, 1247960, 1769758, 1708715, 1708715, 1925741};

#include <htslib/hts.h>  // for htsFormat, seq_nt16_str
#include <htslib/sam.h>  // for bam1_t, bam1_core_t, sam_write1
#include <signal.h>      // for sigaction, sigemptyset
#include <stdint.h>      // for uint8_t
#include <strings.h>     // for strcasecmp
#include <sys/signal.h>  // for sigaction, SIGINT, SIGPIPE, sa_...
#include <time.h>        // for time, time_t
#include <zlib.h>        // for gzprintf, gzFile

#include <algorithm>
#include <cassert>  // for assert
#include <cmath>    // for pow
#include <cstdio>   // for fprintf, stderr, NULL, FILE
#include <cstdlib>  // for free, exit, calloc
#include <cstring>  // for strdup, strlen, strcmp
#include <map>      // for __map_iterator, operator!=, ope...
#include <vector>   // for vector, __vector_base<>::value_...

#include "ngsLCA_cli.h"  // for pars, get_pars, pars_free, prin...
#include "profile.h"     // for damage, destroy_damage, bam_is_...
#include "shared.h"      // for queue, bamRefId2tax, destroy_queue
#include "types.h"       // for int2int, int2char, char2int
#include "version.h"     // for METADAMAGE_VERSION

extern int SIG_COND;  // if we catch signal then quit program nicely
int VERBOSE = 1;
int really_kill = 3;
int2int errmap;
void handler(int s) {
    if (VERBOSE)
        fprintf(stderr, "\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");

    if (--really_kill != 3)
        fprintf(stderr, "\n\t-> If you really want ./ngsLCA to exit uncleanly ctrl+c: %d more times\n", really_kill + 1);
    fflush(stderr);
    if (!really_kill)
        exit(0);
    VERBOSE = 0;
    SIG_COND = 0;
}
// bit of a hang and counter intuitive
/*
  We associate a intrank to each char*rank, the value will be the intrank for the key which is the lcarank that we use of the classification.
*/
char2int setlevels(int norank2species, char *key, int &value) {
    const char *names[48] = {"superkingdom", "domain", "lineage", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "clade", "cohort", "subcohort", "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "infratribe", "genus", "subgenus", "section", "series", "subseries", "subsection", "species", "species group", "species subgroup", "subspecies", "varietas", "morph", "subvariety", "forma", "forma specialis", "biotype", "genotype", "isolate", "pathogroup", "serogroup", "serotype", "strain"};
    int values[48] = {36, 37, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    char2int c2i;
    value = -1;
    char *spec = strdup("species");
    for (int i = 0; i < 48; i++) {
        c2i[strdup(names[i])] = values[i];
        if (strcmp(key, names[i]) == 0) {
            value = values[i];
        }
    }
    if (value == -1) {
        fprintf(stderr, "\t-> Unknown level defined for -lca_rank: %s\n", key);
        exit(0);
    }
    if (norank2species) {
      char2int::iterator it = c2i.find(spec);
      assert(it != c2i.end());
      c2i[strdup("no rank")] = it->second;
    }
    fprintf(stderr, "\t-> Number of entries with level information: %lu \n", c2i.size());
    free(spec);
    return c2i;
}

int2int rank2level(int2char &i2c, int norank2species, char *key, int &value) {
    char2int c2i = setlevels(norank2species, key, value);
    int2int i2i;
    for (int2char::iterator it = i2c.begin(); it != i2c.end(); it++) {
        int key = -1;
        char2int::iterator it2 = c2i.find(it->second);
        if (it2 != c2i.end())
            key = it2->second;
        else {
            if (strcmp(it->second, "no rank") != 0)
                fprintf(stderr, "\t-> Problem finding level for rank: %s\n", it->second);
        }
        i2i[it->first] = key;
    }
    for (char2int::iterator it = c2i.begin(); it != c2i.end(); it++)
        free(it->first);
    return i2i;
}

// we are threading so we want make a nice signal handler for ctrl+c
void catchkill() {
    struct sigaction sa;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sa.sa_handler = handler;
    sigaction(SIGPIPE, &sa, 0);
    sigaction(SIGINT, &sa, 0);
}

int2int specWeight;    // Number of reads that map uniquely to a species.
int2int i2i_missing;   // contains counter of missing hits for each taxid that doesnt exists in acc2taxid
char2int c2i_missing;  // contains counter of missing hits for each taxid that doesnt exists in acc2taxid

void mod_db(int *in, int *out, int2int &parent, int2char &rank, int2char &name_map) {
    int2int::iterator iti;
    int2char::iterator itc;
    for (int i = 0; i < 24; i++) {
        if (parent.count(out[i]) != 1) {
            fprintf(stderr, "\t-> Problem \"fixing\" database entries with known issues, consider add -fix_ncbi 0 when running program\n");
            exit(0);
        }
        iti = parent.find(in[i]);
        if (iti != parent.end()) {
            int oldval_int = iti->second;
            parent.erase(iti);
            parent[out[i]] = oldval_int;
        }
        itc = rank.find(in[i]);
        if (itc != rank.end()) {
            char *oldval_char = itc->second;
            rank.erase(itc);
            rank[out[i]] = oldval_char;
        }
        name_map[in[i]] = strdup("satan");
    }
}

int nodes2root(int taxa, int2int &parent) {
    int dist = 0;
    while (taxa != 1) {
        int2int::iterator it = parent.find(taxa);
        taxa = it->second;
        dist++;
    }
    return dist;
}

int2int dist2root;

int satan(int taxid, int2int &parant) {
    int2int::iterator it = dist2root.find(taxid);

    return 0;
}

float mean(std::vector<float> &vec) {
    float tmp = vec[0];
    for (int i = 1; i < vec.size(); i++)
        tmp += vec[i];
    return tmp / vec.size();
}

int varinfo = 1;
float var(std::vector<float> &vec) {
    if (vec.size() <= 1) {
        if (varinfo == 1) {
            fprintf(stderr, "\t-> Calculation of variance is only defined for >1 datapoints (this message is only printed once)\n");
            varinfo = 0;
        }
        return 0.0;
    }
    float mea = mean(vec);
    float tmp = 0;
    for (int i = 1; i < vec.size(); i++)
        tmp += pow(vec[i] - mea, 2);

    return tmp / (vec.size() - 1);
}

// int2int lcadist;
typedef struct {
    int nalignments;
    std::vector<float> readlengths;
    std::vector<float> gccontents;
} lcatriplet;

std::map<int, lcatriplet> lcastat;

void adder(int taxid, int readlengths, float gccontent) {
    std::map<int, lcatriplet>::iterator it = lcastat.find(taxid);
    if (it == lcastat.end()) {
        lcatriplet tmp;
        tmp.nalignments = 1;
        tmp.readlengths.push_back(readlengths);
        tmp.gccontents.push_back(gccontent);
        lcastat[taxid] = tmp;
    } else {
        it->second.nalignments = it->second.nalignments + 1;
        it->second.readlengths.push_back(readlengths);
        it->second.gccontents.push_back(gccontent);
    }
}
int do_lca(std::vector<int> &taxids, int2int &parent) {
    //  fprintf(stderr,"\t-> [%s] with number of taxids: %lu\n",__func__,taxids.size());
    assert(taxids.size() > 0);
    if (taxids.size() == 1) {
        int taxa = taxids[0];
        if (parent.count(taxa) == 0) {
            fprintf(stderr, "\t-> Problem finding taxaid: %d will skip\n", taxa);
            taxids.clear();
            return -1;
        }

        taxids.clear();
        return taxa;
    }

    int2int counter;
    for (int i = 0; i < taxids.size(); i++) {
        int taxa = taxids[i];
        while (1) {
            //      fprintf(stderr,"taxa:%d\n",taxa);
            int2int::iterator it = counter.find(taxa);
            if (it == counter.end()) {
                //	fprintf(stderr,"taxa: %d is new will plugin\n",taxa);
                counter[taxa] = 1;
            } else
                it->second = it->second + 1;
            it = parent.find(taxa);

            if (it == parent.end()) {
                int2int::iterator it = errmap.find(taxa);
                if (it == errmap.end()) {
                    fprintf(stderr, "\t-> Problem finding parent of :%d\n", taxa);

                    errmap[taxa] = 1;
                } else
                    it->second = it->second + 1;
                taxids.clear();
                return -1;
            }

            if (taxa == it->second)  //<- catch root
                break;
            taxa = it->second;
        }
    }
    //  fprintf(stderr,"counter.size():%lu\n",counter.size());
    // now counts contain how many time a node is traversed to the root
    int2int dist2root;
    for (int2int::iterator it = counter.begin(); it != counter.end(); it++)
        if (it->second == taxids.size())
            dist2root[nodes2root(it->first, parent)] = it->first;
    for (int2int::iterator it = dist2root.begin(); 0 && it != dist2root.end(); it++)
        fprintf(stderr, "%d\t->%d\n", it->first, it->second);
    taxids.clear();
    if (!dist2root.empty())
        return (--dist2root.end())->second;
    else {
        fprintf(stderr, "\t-> Happens\n");
        return -1;
    }
}

void print_chain1(kstring_t *kstr, int taxa, int2char &rank, int2char &name_map,char sep) {
    int2char::iterator it1 = name_map.find(taxa);
    int2char::iterator it2 = rank.find(taxa);
    if (it1 == name_map.end()) {
        fprintf(stderr, "\t-> Problem finding taxaid:%d\n", taxa);
    }
    assert(it2 != rank.end());
    if (it1 == name_map.end() || it2 == rank.end()) {
        fprintf(stderr, "taxa: %d %s doesnt exists will exit\n", taxa, it1->second);
        exit(1);
    }
    ksprintf(kstr, "%c%d:\"%s\":\"%s\"",sep, taxa, it1->second, it2->second);
}

void print_rank(FILE *fp, int taxa, int2char &rank) {
    int2char::iterator it2 = rank.find(taxa);
    assert(it2 != rank.end());
    fprintf(stderr, "taxa: %d rank %s\n", taxa, it2->second);
}

void print_chain(kstring_t *kstr, int taxa, int2int &parent, int2char &rank, int2char &name_map,int donewline) {
  int first = 0;
  char sep = '\t';
  print_chain1(kstr, taxa, rank, name_map,sep);
  while (1) {
    print_chain1(kstr, taxa, rank, name_map,sep);
    int2int::iterator it = parent.find(taxa);
    assert(it != parent.end());
    if (taxa == it->second)  //<- catch root
      break;
    taxa = it->second;

    if(sep=='\t'){
      if(++first > 0)
	sep = ';';
    }
  }
  if(donewline)
    ksprintf(kstr,"\n");
}

int isuniq(std::vector<int> &vec) {
    int ret = 1;
    for (int i = 1; i < vec.size(); i++)
        if (vec[0] != vec[i])
            return 0;
    return ret;
}

int get_species1(int taxa, int2int &parent, int2char &rank) {
    //  fprintf(stderr,"\t-> taxa:%d\n",taxa);
    int2char::iterator it;
    while (1) {
        it = rank.find(taxa);
        if (it == rank.end())
            return -1;
        char *val = it->second;
        if (val == NULL)
            fprintf(stderr, "taxa:%d\n", taxa);
        assert(val);
        if (!strcmp(val, "species"))
            return taxa;
        int next = parent[taxa];
        if (next == taxa)
            break;
        taxa = next;
    }
    return taxa;
}

int2int get_species(int2int &i2i, int2int &parent, int2char &rank, int2char &names, FILE *fp) {
    int2int ret;
    for (int2int::iterator it = i2i.begin(); it != i2i.end(); it++) {
        // fprintf(stderr,"%d\t%d\n",it->first,it->second);
        int asdf = get_species1(it->second, parent, rank);
        if (asdf == -1) {
            fprintf(fp, "\t-> Removing pair(%d,%d) accnumber:%s since doesnt exists in node list\n", it->first, it->second, names[it->second]);
            i2i.erase(it);
        } else
            ret[it->second] = asdf;  //
    }
    return ret;
}
//
char *make_seq(bam1_t *aln) {
    int len = aln->core.l_qseq;
    char *qseq = new char[len + 1];
    uint8_t *q = bam_get_seq(aln);
    for (int i = 0; i < len; i++)
        qseq[i] = seq_nt16_str[bam_seqi(q, i)];
    qseq[len] = '\0';
    //  fprintf(stderr,"seq:%s\n",qseq);
    // exit(0);
    return qseq;
}

// a,c,g,t,n
// A,C,G,T,N
// 0,1,2,3,4
int refToInt[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 47
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 79
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 95
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 111
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 127
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 175
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 207
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4   // 255
};

float gccontent(char *seq) {
    int counts[5] = {0, 0, 0, 0, 0};
    for (int i = 0; i < strlen(seq); i++)
        counts[refToInt[seq[i]]]++;

    float gcs = counts[1] + counts[2];
    float tot = counts[0] + counts[1] + counts[2] + counts[3];
    return gcs / tot;
}

float gccontent(bam1_t *aln) {
    int counts[5] = {0, 0, 0, 0, 0};

    int len = aln->core.l_qseq;
    uint8_t *q = bam_get_seq(aln);
    for (int i = 0; i < len; i++)
        counts[refToInt[seq_nt16_str[bam_seqi(q, i)]]]++;

    float gcs = counts[1] + counts[2];
    float tot = counts[0] + counts[1] + counts[2] + counts[3];
    return gcs / tot;
}

int printonce = 1;

std::vector<int> purge(std::vector<int> &taxids, std::vector<int> &editdist) {
    if (printonce == 1)
        fprintf(stderr, "\t-> purging taxids oldsize:%lu\n", taxids.size());
    assert(taxids.size() == editdist.size());
    std::vector<int> tmpnewvec;
    int mylow = *std::min_element(editdist.begin(), editdist.end());
    for (int i = 0; i < taxids.size(); i++)
        if (editdist[i] <= mylow)
            tmpnewvec.push_back(taxids[i]);
    if (printonce-- == 1)
        fprintf(stderr, "\t-> purging taxids newsize:%lu this info is only printed once\n", tmpnewvec.size());
    return tmpnewvec;
}

void hts(gzFile fp, samFile *fp_in, int2int &i2i, int2int &parent, bam_hdr_t *hdr, int2char &rank, int2char &name_map, int minmapq, int discard, int editMin, int editMax, double scoreLow, double scoreHigh, int minlength, int lca_rank, char *prefix, int howmany, samFile *fp_usedreads, int skipnorank, int2int &rank2level, int nthreads, int weighttype,long maxreads,samFile *fp_famout) {
  fprintf(stderr, "[%s] \t-> editMin:%d editmMax:%d scoreLow:%f scoreHigh:%f minlength:%d discard: %d prefix: %s howmany: %d skipnorank: %d weighttype: %d maxreads: %ld\n", __FUNCTION__, editMin, editMax, scoreLow, scoreHigh, minlength, discard, prefix, howmany, skipnorank, weighttype,maxreads);
    assert(fp_in != NULL);
    damage *dmg = new damage(howmany, nthreads, 13);
    bam1_t *aln = bam_init1();  // initialize an alignment
    int comp;

    char *last = NULL;
    char *seq = NULL;
    std::vector<int> taxids;
    std::vector<int> specs;
    std::vector<int> editdist;
    queue *myq = init_queue(500);

    int lca;
    int2int closest_species;
    int skip = 0;
    int inc = 0;
    kstring_t *kstr = new kstring_t;
    kstr->s = NULL;kstr->l = kstr->m =0;
    long nreads = 0;
    while (sam_read1(fp_in, hdr, aln) >= 0) {
      if(maxreads!=-1&&nreads>=maxreads)
	break;
      if (bam_is_unmapped(aln)) {
	// fprintf(stderr,"skipping: %s unmapped \n",bam_get_qname(b));
	continue;
      }
        if (bam_is_failed(aln)) {
            // fprintf(stderr,"skipping: %s failed: flags=%d \n",bam_get_qname(b),b->core.flag);
            continue;
        }
        char *qname = bam_get_qname(aln);
        int chr = aln->core.tid;  // contig name (chromosome)
        //    fprintf(stderr,"%d %d\n",aln->core.qual,minmapq);
        static int ntimes = 3;
        if (aln->core.qual < minmapq && ntimes > 0) {
            ntimes--;
            fprintf(stderr, "Discarding due to low mapq, this message will only be printed three times\n");
            continue;
        }

        if (discard > 0 && (aln->core.flag & discard)) {
            fprintf(stderr, "Discard: %d coreflag:%d OR %d\n", discard, aln->core.flag, aln->core.flag & discard);
            fprintf(stderr, "Discarding due to core flag\n");
            continue;
        }
        if (last == NULL) {
            last = strdup(qname);
            seq = make_seq(aln);
        }
        if (minlength != -1 && (aln->core.l_qseq < minlength))
            continue;
        // change of ref
        if (strcmp(last, qname) != 0) {
            if (taxids.size() > 0 && skip == 0) {
                //	fprintf(stderr,"length of taxids:%lu and other:%lu minedit:%d\n",taxids.size(),editdist.size(),*std::min_element(editdist.begin(),editdist.end()));

                int size = taxids.size();
                assert(size == myq->l);
                if (editMin == -1 && editMax == -1)
                    taxids = purge(taxids, editdist);
		nreads++;
                lca = do_lca(taxids, parent);
                //	fprintf(stderr,"myq->l: %d\n",myq->l);
                if (lca != -1) {
                    gzprintf(fp, "%s\t%s\t%lu\t%d\t%f", last, seq, strlen(seq), size, gccontent(seq));
		   
		    print_chain(kstr, lca, parent, rank, name_map,1);
		    gzwrite(fp,kstr->s,kstr->l);
		    kstr->l =0;

                    int varisunique = isuniq(specs);
                    if (varisunique) {
                        int2int::iterator it = specWeight.find(specs[0]);
                        if (it == specWeight.end())
                            specWeight[specs[0]] = 1;
                        else
                            it->second = it->second + 1;
                    }

                    // lca_rank is integer and is different from minus one
                    int2int::iterator myit = rank2level.find(lca);
                    assert(myit != rank2level.end());
		    //write to family out
		    if (myit->second != -1 && (myit->second <= 16)) {
		      //family is 16, see rank2level.txt
		       for (int i = 0; i < myq->l; i++)
			 if(fp_famout)
			   assert(sam_write1(fp_famout, hdr, myq->ary[i]) >= 0);
		    }
		    //standard analyses
		    if (myit->second != -1 && (myit->second <= lca_rank)) {
                        adder(lca, strlen(seq), gccontent(seq));
			            //fprintf(stderr,"Looping through alignments we have :%d \n",myq->l);
                        for (int i = 0; i < myq->l; i++) {
			  
                            int2int::iterator it2k = i2i.find(myq->ary[i]->core.tid);
                            assert(it2k != i2i.end());
                            int2char::iterator ititit = rank.find(it2k->second);
                            if (ititit == rank.end()) {
                                fprintf(stderr, "\t-> Potential problem no rank for taxid: %d\n", it2k->second);
                                continue;
                            }
                            if (skipnorank == 1 && strcasecmp(ititit->second, "no rank") == 0)
                                continue;
                            dmg->damage_analysis(myq->ary[i], it2k->second, weighttype == 0 ? 1 : (1.0 / (float)myq->l));
                            if (fp_usedreads)
                                assert(sam_write1(fp_usedreads, hdr, myq->ary[i]) >= 0);
                        }
                    }
                }
            }
            skip = 0;
            specs.clear();
            editdist.clear();
            myq->l = 0;
            free(last);
            delete[] seq;
            last = strdup(qname);
            seq = make_seq(aln);
        }

        // filter by nm
        uint8_t *nm = bam_aux_get(aln, "NM");
        int thiseditdist;
        if (nm != NULL) {
            thiseditdist = (int)bam_aux2i(nm);
            //      fprintf(stderr,"[%d] nm:%d\t",inc++,val);
            if (editMin != -1 && thiseditdist < editMin) {
                skip = 1;
                // fprintf(stderr,"skipped1\n");
                continue;
            } else if (editMax != -1 && thiseditdist > editMax) {
                // fprintf(stderr,"continued1\n");
                continue;
            }
            double seqlen = aln->core.l_qseq;
            double myscore = 1.0 - (((double)thiseditdist) / seqlen);
            //      fprintf(stderr," score:%f\t",myscore);
            if (myscore > scoreHigh) {
                // fprintf(stderr,"skipped2\n");
                skip = 1;
                continue;
            } else if (myscore < scoreLow) {
                //	fprintf(stderr,"continued2\n");
                continue;
            }
        }
        int2int::iterator it = i2i.find(chr);
        // See if cloests speciest exists and plug into closests species
        int dingdong = -1;
        if (it != i2i.end()) {
            int2int::iterator it2 = closest_species.find(chr);
            if (it2 != closest_species.end())
                dingdong = it->second;
            else {
                dingdong = get_species1(it->second, parent, rank);
                if (dingdong != -1)
                    closest_species[it->second] = dingdong;
            }
            // fprintf(stderr,"\t-> closests size:%lu\n",closest_species.size());
        }

        if (it == i2i.end() || dingdong == -1) {
	  int2int::iterator it2missing = i2i_missing.find(chr);
            if (it2missing == i2i_missing.end()) {
                fprintf(stderr, "[log] \t-> problem finding chrid:%d chrname:%s (due to problem with missing entries in nodes.dmp or supportting files)\n", chr, hdr->target_name[chr]);
                i2i_missing[chr] = 1;
            } else
                it2missing->second = it2missing->second + 1;
        } else {
            assert(bam_copy1(myq->ary[myq->l], aln) != NULL);
            myq->l++;
            if (myq->l == myq->m)
                expand_queue(myq);
            taxids.push_back(it->second);
            specs.push_back(dingdong);
            editdist.push_back(thiseditdist);

            // fprintf(stderr,"it-.second:%d specs:%d thiseditdist:%d\n",it->second,dingdong,thiseditdist);
            //       fprintf(stderr,"EDIT\t%d\n",thiseditdist);
        }
    }
    if (taxids.size() > 0 && skip == 0) {
        int size = taxids.size();
        assert(size == myq->l);
        if (editMin == -1 && editMax == -1)
            taxids = purge(taxids, editdist);

        //  fprintf(stderr,"ntaxids: %lu\n",taxids.size());
	nreads++;
        lca = do_lca(taxids, parent);
        //    fprintf(stderr,"myq->l: %d lca: %d \n",myq->l,lca);
        if (lca != -1) {
            gzprintf(fp, "%s\t%s\t%lu\t%d\t%f", last, seq, strlen(seq), size, gccontent(seq));
            print_chain(kstr, lca, parent, rank, name_map,1);
	    gzwrite(fp,kstr->s,kstr->l);
	    kstr->l = 0;
            if (isuniq(specs)) {
                int2int::iterator it = specWeight.find(specs[0]);
                if (it == specWeight.end())
                    specWeight[specs[0]] = specs.size();
                else
                    it->second = it->second + specs.size();
            }
            // lca_rank is integer and is different from minus one
            int2int::iterator myit = rank2level.find(lca);
            assert(myit != rank2level.end());
	    //write to family out
	    if (myit->second != -1 && (myit->second <= 16)) {
	      //family is 16, see rank2level.txt
	      for (int i = 0; i < myq->l; i++)
		if(fp_famout)
		  assert(sam_write1(fp_famout, hdr, myq->ary[i]) >= 0);
	    }
            if (myit->second != -1 && (myit->second <= lca_rank)) {
                adder(lca, strlen(seq), gccontent(seq));
                //      if(correct_rank(lca_rank,lca,rank,norank2species)){
		        //fprintf(stderr,"Looping through alignments we have :%d \n",myq->l);
                for (int i = 0; i < myq->l; i++) {
                    // dmg->damage_analysis(myq->ary[i],myq->ary[i]->core.tid);
                    int2int::iterator ittt = i2i.find(myq->ary[i]->core.tid);
                    assert(ittt != i2i.end());
                    int2char::iterator ititit = rank.find(ittt->second);
                    if (ititit == rank.end()) {
                        fprintf(stderr, "\t-> Potential problem no rank for taxid: %d\n", ititit->first);
                        continue;
                    }
                    //	  fprintf(stderr,"uaua: %s taxid: %d\n",ititit->second,ititit->first);
                    if (skipnorank == 1 && strcasecmp(ititit->second, "no rank") == 0)
                        continue;
                    //	  fprintf(stderr,"uaua\n");
                    dmg->damage_analysis(myq->ary[i], ittt->second, weighttype == 0 ? 1 : (1.0 / (float)myq->l));
                    if (fp_usedreads)
                        assert(sam_write1(fp_usedreads, hdr, myq->ary[i]) >= 0);
                }
            }
        }
    }
    fprintf(stderr,"maxreads: %ld nreads: %ld\n",maxreads,nreads);
    dmg->bwrite(prefix);

    specs.clear();
    editdist.clear();
    bam_destroy1(aln);
    sam_close(fp_in);

    if (seq)
        delete[] seq;
    if (last)
        free(last);
    destroy_damage(dmg);
    destroy_queue(myq);
    free(kstr->s);delete kstr;
    return;  // 0;
}

void print_ref_rank_species(bam_hdr_t *h, int2int &i2i, int2char &names, int2char &rank) {
    for (int i = 0; i < h->n_targets; i++) {
        fprintf(stdout, "%d\t%s\t%s\n", i, names[i2i[i]], rank[i2i[i]]);
    }
}

int calc_valens(int2int &i2i, int2int &parent) {
    int2int ret;
    for (int2int::iterator it = i2i.begin(); it != i2i.end(); it++) {
        int2int::iterator pit = parent.find(it->second);
        if (pit == parent.end() || pit->second == 1)
            continue;
        //    fprintf(stdout,"%d %d\n",pit->first,pit->second);
        int2int::iterator rit = ret.find(pit->second);
        if (rit == ret.end())
            ret[pit->second] = 1;
        else
            rit->second = rit->second + 1;
    }

    for (int2int::iterator it = ret.begin(); it != ret.end(); it++)
        fprintf(stdout, "%d\t%d\n", it->first, it->second);
    return 0;
}

int calc_dist2root(int2int &i2i, int2int &parent) {
    int2int ret;
    for (int2int::iterator it = i2i.begin(); it != i2i.end(); it++) {
        int2int::iterator pit = parent.find(it->second);
        if (pit == parent.end() || pit->second == 1)
            continue;
        //    fprintf(stdout,"%d %d\n",pit->first,pit->second);
        int2int::iterator rit = ret.find(pit->second);
        if (rit == ret.end())
            ret[pit->second] = 1;
        else
            rit->second = rit->second + 1;
    }

    for (int2int::iterator it = ret.begin(); it != ret.end(); it++)
        fprintf(stdout, "%d\t%d\n", it->first, it->second);
    return 0;
}
#if 0
int2node makeNodes(int2int &parent){
  int2node ret;
  for(int2int::iterator it=parent.begin();it!=parent.end();it++){
    int2node::iterator it2=ret.find(it->second);
    if(it2==ret.end()){
      node nd;
      nd.up=it->second;
      ret[it->first] = nd;
    }
    
  }
  return ret;
}
#endif

int main_lca(int argc, char **argv) {
    htsFormat *dingding2 = (htsFormat *)calloc(1, sizeof(htsFormat));
    if (argc == 1) {
        fprintf(stderr, "\t-> ./metaDMG-cpp lca --names --nodes --acc2tax [-edit_dist_[min/max] --sim_score_[low/high] --min_mapq --bam --lca_rank  --filtered_acc2tax --used_reads [0,1] --weight_type [0,1] --famout [0,1] #%s \n", METADAMAGE_VERSION);
        return 0;
    }
    catchkill();
#if 0
  int2int ww = get_weight(argv[1]);
  return 0;
#endif
    time_t t2 = time(NULL);

    pars *p = get_pars(--argc, ++argv);
    print_pars(stderr, p);
    gzprintf(p->fp1,"#./metaDMG-cpp lca ");
    for(int i=0;i<argc;i++)
      gzprintf(p->fp1,"%s ",argv[i]);
    gzprintf(p->fp1,"\n");
    // map of bamref ->taxid
    int2int *i2i = NULL;
    //  fprintf(stderr,"p->header: %p\n",p->header);
    if (p->header)
      i2i = bamRefId2tax(p->header, p->acc2taxfile, p->htsfile, errmap, p->tempfolder, p->useDump, p->filteredAcc2taxfile,NULL);
    // map of taxid -> taxid
    int2int parent;
    // map of taxid -> rank
    int2char rank;
    // map of taxid -> name
    // map of parent -> child taxids
    int2intvec child;

    int2char name_map = parse_names(p->namesfile);

    parse_nodes(p->nodesfile, rank, parent, child, 0);
    int lca_rank;
    // converts each taxid to a intrepresentation of the rank, we call this level
    int2int tax2level = rank2level(rank, p->norank2species, p->lca_rank, lca_rank);
    //  calc_valens(i2i,parent);
    if (0) {
        print_ref_rank_species(p->header, *i2i, name_map, rank);
        return 0;
    }

    // closes species (direction of root) for a taxid
    //   int2int closest_species=get_species(i2i,parent,rank,name_map,p->fp3);
    //   fprintf(stderr,"\t-> Number of items in closest_species map:%lu\n",closest_species.size());

    if (p->fixdb) {
        fprintf(stderr, "\t-> Will add some fixes of the ncbi database due to merged names\n");
        mod_db(mod_in, mod_out, parent, rank, name_map);
    }
    samFile *usedreads_sam = NULL;
    if (p->usedreads_sam != NULL) {  // p->usedreads sam is const *, sorry this is confusing
        char out_mode[5] = "wb";
        if ((usedreads_sam = sam_open_format(p->usedreads_sam, out_mode, dingding2)) == 0) {
            fprintf(stderr, "Error opening file for writing\n");
            exit(0);
        }
        if (sam_hdr_write(usedreads_sam, p->header) < 0)
            fprintf(stderr, "writing headers to %s", p->usedreads_sam);
    }
    samFile *famout_sam = NULL;
    if (p->famout_sam != NULL) {  // p->usedreads sam is const *, sorry this is confusing
      char out_mode[5] = "wb";
      if ((famout_sam = sam_open_format(p->famout_sam, out_mode, dingding2)) == 0) {
	fprintf(stderr, "Error opening file for writing\n");
	exit(0);
      }
      if (sam_hdr_write(famout_sam, p->header) < 0)
	fprintf(stderr, "writing headers to %s", p->famout_sam);
    }
    gzprintf(p->fp1,"queryid\tseq\tlen\tnaln\tgc\tlca\ttaxa_path\n");
    hts(p->fp1, p->hts, *i2i, parent, p->header, rank, name_map, p->minmapq, p->discard, p->editdistMin, p->editdistMax, p->simscoreLow, p->simscoreHigh, p->minlength, lca_rank, p->outnames, p->howmany, usedreads_sam, p->skipnorank, tax2level, p->nthreads, p->weighttype,p->maxreads,famout_sam);

    fprintf(stderr, "\t-> Number of species with reads that map uniquely: %lu\n", specWeight.size());

    for (int2int::iterator it = errmap.begin(); it != errmap.end(); it++)
        fprintf(stderr, "[log] err\t%d\t%d\n", it->first, it->second);
#if 0
  for(int2int::iterator it=i2i_missing.begin();it!=i2i_missing.end();it++)
    fprintf(p->fp3,"missingtaxid \t%d\t%d\t%s\n",it->first,it->second,p->header[it->first]);
#endif
    // p->header points to bam_hdr_t what is expected here?
    for (int2int::iterator it = specWeight.begin(); 0 && it != specWeight.end(); it++)
      gzprintf(p->fp2, "%d\t%s\t%d\n", it->first, name_map[it->first], it->second);

    fprintf(stderr, "\t-> [ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));

    if (usedreads_sam != NULL)
        sam_close(usedreads_sam);
    if (famout_sam != NULL)
      sam_close(famout_sam);
    if (p->fp_lcadist) {
        gzprintf(p->fp_lcadist,"taxid\tnreads\tmea_len\tvar_len\tmean_gc\tvar_gc\tlca\trank\n");
        for (std::map<int, lcatriplet>::iterator it = lcastat.begin(); it != lcastat.end(); it++) {
            lcatriplet tmp = it->second;
            gzprintf(p->fp_lcadist, "%d\t%d\t%f\t%f\t%f\t%f", it->first, tmp.nalignments, mean(tmp.readlengths), var(tmp.readlengths), mean(tmp.gccontents), var(tmp.gccontents));
            int2char::iterator it1 = name_map.find(it->first);
            int2char::iterator it2 = rank.find(it->first);
            char *namnam, *rankrank;
            namnam = rankrank = NULL;
            if (it1 != name_map.end())
                namnam = it1->second;
            if (it2 != rank.end())
                rankrank = it2->second;
            gzprintf(p->fp_lcadist, "\t\"%s\"\t\"%s\"\n", namnam, rankrank);
        }
    }

    for (int2char::iterator it = name_map.begin(); it != name_map.end(); it++)
        free(it->second);
    for (int2char::iterator it = rank.begin(); it != rank.end(); it++)
        free(it->second);
    if (i2i)
        delete i2i;
    free(dingding2);
    pars_free(p);
    return 0;
}
