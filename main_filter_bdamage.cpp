#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <strings.h>
#include <zlib.h>

#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "main_filter_bdamage.h"
#include "shared.h"

struct filter_bdamage_args {
    char *infile;
    char *out_prefix;
    char *id_file;
    char *resolve_file;
    char *bam_file;
    char *acc2tax_file;
    char *nodes_file;
    char *stat_infile;
    char *rlens_infile;
    int min_reads;
    int max_reads;
    int nthreads;
    int exclude;
    int strict_resolve;
    std::set<int> ids;
    std::vector<std::string> selectors;
};

static int usage_filter_bdamage(FILE *fp) {
    fprintf(fp, "Usage: ./metaDMG-cpp filter_bdamage <in.bdamage.gz> [options]\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "  -o/--out_prefix STR     output prefix (default: filtered)\n");
    fprintf(fp, "  -n/--threads INT        BGZF read/write threads (default: 4)\n");
    fprintf(fp, "  --id INT                include this numeric id (repeatable)\n");
    fprintf(fp, "  --id_list STR           comma-separated numeric ids/selectors\n");
    fprintf(fp, "  --id_file FILE          newline separated ids/selectors (first column used)\n");
    fprintf(fp, "  --resolve STR           selector to resolve (repeatable)\n");
    fprintf(fp, "  --resolve_list STR      comma-separated selectors to resolve\n");
    fprintf(fp, "  --resolve_file FILE     newline separated selectors (first column used)\n");
    fprintf(fp, "  --bam FILE              BAM/CRAM/SAM for refname -> bam-offset (tid) resolve\n");
    fprintf(fp, "  --acc2tax FILE          accession2taxid map for accession -> taxid resolve\n");
    fprintf(fp, "  --nodes FILE            expand selected taxid to all descendants from nodes file\n");
    fprintf(fp, "  --stat FILE             companion stat.gz to filter (auto-inferred if omitted)\n");
    fprintf(fp, "  --rlens FILE            companion rlens.gz to filter (auto-inferred if omitted)\n");
    fprintf(fp, "  --strict_resolve INT    1: fail on unresolved selectors (default), 0: ignore\n");
    fprintf(fp, "  --exclude               exclude listed ids/selectors instead of including\n");
    fprintf(fp, "  --min_reads INT         keep entries with nreads >= INT\n");
    fprintf(fp, "  --max_reads INT         keep entries with nreads <= INT\n");
    fprintf(fp, "  -h/--help               show this help\n");
    fprintf(fp, "\nOutput:\n");
    fprintf(fp, "  <out_prefix>.bdamage.gz\n");
    fprintf(fp, "  <out_prefix>.stat.gz  (if companion found/specified)\n");
    fprintf(fp, "  <out_prefix>.rlens.gz (if companion found/specified)\n");
    fprintf(fp, "\nExamples:\n");
    fprintf(fp, "  ./metaDMG-cpp filter_bdamage in.bdamage.gz --id 11 --out_prefix subset\n");
    fprintf(fp, "  ./metaDMG-cpp filter_bdamage in.bdamage.gz --resolve ref2 --bam reads.bam --out_prefix subset\n");
    fprintf(fp, "  ./metaDMG-cpp filter_bdamage in.bdamage.gz --resolve ACC123 --acc2tax acc2taxid.map.gz --nodes nodes.dmp.gz --out_prefix subset\n");
    fprintf(fp, "  ./metaDMG-cpp filter_bdamage in.bdamage.gz --id_list 2,10,42 --exclude --min_reads 100 --out_prefix subset\n");
    return 0;
}

static filter_bdamage_args init_filter_bdamage_args() {
    filter_bdamage_args args;
    args.infile = NULL;
    args.out_prefix = strdup("filtered");
    args.id_file = NULL;
    args.resolve_file = NULL;
    args.bam_file = NULL;
    args.acc2tax_file = NULL;
    args.nodes_file = NULL;
    args.stat_infile = NULL;
    args.rlens_infile = NULL;
    args.min_reads = -1;
    args.max_reads = -1;
    args.nthreads = 4;
    args.exclude = 0;
    args.strict_resolve = 1;
    return args;
}

static int parse_int(const char *txt, int *out) {
    if (txt == NULL || out == NULL)
        return 1;

    errno = 0;
    char *end = NULL;
    long val = strtol(txt, &end, 10);
    if (end == txt || *end != '\0' || errno != 0 || val < INT_MIN || val > INT_MAX)
        return 1;
    *out = (int)val;
    return 0;
}

static int has_suffix(const std::string &txt, const std::string &suffix) {
    if (txt.size() < suffix.size())
        return 0;
    return txt.compare(txt.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static void add_selector_token(const char *tok, filter_bdamage_args &args) {
    if (tok == NULL || tok[0] == '\0')
        return;

    int id = 0;
    if (parse_int(tok, &id) == 0)
        args.ids.insert(id);
    else
        args.selectors.push_back(std::string(tok));
}

static int add_csv_tokens(const char *csv, filter_bdamage_args &args) {
    if (csv == NULL)
        return 1;

    char *tmp = strdup(csv);
    char *token = strtok(tmp, ",");
    while (token != NULL) {
        add_selector_token(token, args);
        token = strtok(NULL, ",");
    }
    free(tmp);
    return 0;
}

static int load_tokens_from_file(const char *path, filter_bdamage_args &args) {
    if (path == NULL)
        return 0;

    FILE *fp = fopen(path, "r");
    if (fp == NULL) {
        fprintf(stderr, "\t-> Error: could not open file: %s\n", path);
        return 1;
    }

    char line[4096];
    while (fgets(line, sizeof(line), fp) != NULL) {
        char *ptr = line;
        while (*ptr == ' ' || *ptr == '\t')
            ptr++;
        if (*ptr == '\0' || *ptr == '\n' || *ptr == '#')
            continue;

        char *tok = strtok(ptr, "\t\n ");
        if (tok == NULL)
            continue;

        add_selector_token(tok, args);
    }

    fclose(fp);
    return 0;
}

static int parse_filter_bdamage_args(int argc, char **argv, filter_bdamage_args &args) {
    for (int i = 1; i < argc; i++) {
        if (!strcasecmp(argv[i], "-h") || !strcasecmp(argv[i], "--help")) {
            usage_filter_bdamage(stdout);
            return 2;
        } else if (!strcasecmp(argv[i], "-o") || !strcasecmp(argv[i], "--out_prefix")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for %s\n", argv[i]);
                return 1;
            }
            free(args.out_prefix);
            args.out_prefix = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "-n") || !strcasecmp(argv[i], "--threads")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for %s\n", argv[i]);
                return 1;
            }
            args.nthreads = atoi(argv[++i]);
        } else if (!strcasecmp(argv[i], "--id")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --id\n");
                return 1;
            }
            int id = 0;
            if (parse_int(argv[++i], &id) != 0) {
                fprintf(stderr, "\t-> Error: invalid value for --id: %s\n", argv[i]);
                return 1;
            }
            args.ids.insert(id);
        } else if (!strcasecmp(argv[i], "--id_list")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --id_list\n");
                return 1;
            }
            if (add_csv_tokens(argv[++i], args) != 0)
                return 1;
        } else if (!strcasecmp(argv[i], "--id_file")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --id_file\n");
                return 1;
            }
            free(args.id_file);
            args.id_file = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "--resolve")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --resolve\n");
                return 1;
            }
            add_selector_token(argv[++i], args);
        } else if (!strcasecmp(argv[i], "--resolve_list")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --resolve_list\n");
                return 1;
            }
            if (add_csv_tokens(argv[++i], args) != 0)
                return 1;
        } else if (!strcasecmp(argv[i], "--resolve_file")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --resolve_file\n");
                return 1;
            }
            free(args.resolve_file);
            args.resolve_file = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "--bam")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --bam\n");
                return 1;
            }
            free(args.bam_file);
            args.bam_file = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "--acc2tax")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --acc2tax\n");
                return 1;
            }
            free(args.acc2tax_file);
            args.acc2tax_file = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "--nodes")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --nodes\n");
                return 1;
            }
            free(args.nodes_file);
            args.nodes_file = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "--stat")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --stat\n");
                return 1;
            }
            free(args.stat_infile);
            args.stat_infile = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "--rlens")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --rlens\n");
                return 1;
            }
            free(args.rlens_infile);
            args.rlens_infile = strdup(argv[++i]);
        } else if (!strcasecmp(argv[i], "--strict_resolve")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --strict_resolve\n");
                return 1;
            }
            args.strict_resolve = atoi(argv[++i]);
        } else if (!strcasecmp(argv[i], "--exclude")) {
            args.exclude = 1;
        } else if (!strcasecmp(argv[i], "--min_reads")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --min_reads\n");
                return 1;
            }
            args.min_reads = atoi(argv[++i]);
        } else if (!strcasecmp(argv[i], "--max_reads")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "\t-> Error: missing value for --max_reads\n");
                return 1;
            }
            args.max_reads = atoi(argv[++i]);
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "\t-> Error: unknown option: %s\n", argv[i]);
            return 1;
        } else {
            if (args.infile == NULL)
                args.infile = strdup(argv[i]);
            else {
                fprintf(stderr, "\t-> Error: multiple input files provided: %s and %s\n", args.infile, argv[i]);
                return 1;
            }
        }
    }

    if (args.infile == NULL) {
        fprintf(stderr, "\t-> Error: missing input .bdamage.gz file\n");
        return 1;
    }

    if (args.nthreads < 1) {
        fprintf(stderr, "\t-> Error: --threads must be >= 1\n");
        return 1;
    }
    if (args.min_reads < -1) {
        fprintf(stderr, "\t-> Error: --min_reads must be -1 or greater\n");
        return 1;
    }
    if (args.max_reads < -1) {
        fprintf(stderr, "\t-> Error: --max_reads must be -1 or greater\n");
        return 1;
    }
    if (args.min_reads != -1 && args.max_reads != -1 && args.min_reads > args.max_reads) {
        fprintf(stderr, "\t-> Error: --min_reads cannot be greater than --max_reads\n");
        return 1;
    }
    if (!(args.strict_resolve == 0 || args.strict_resolve == 1)) {
        fprintf(stderr, "\t-> Error: --strict_resolve must be 0 or 1\n");
        return 1;
    }

    if (load_tokens_from_file(args.id_file, args) != 0)
        return 1;
    if (load_tokens_from_file(args.resolve_file, args) != 0)
        return 1;

    return 0;
}

static int maybe_infer_companion_inputs(filter_bdamage_args &args) {
    if (args.infile == NULL)
        return 0;

    std::string in(args.infile);
    std::string base;
    if (has_suffix(in, ".bdamage.gz"))
        base = in.substr(0, in.size() - strlen(".bdamage.gz"));
    else if (has_suffix(in, ".bdamage"))
        base = in.substr(0, in.size() - strlen(".bdamage"));
    else
        return 0;

    if (args.stat_infile == NULL) {
        std::vector<std::string> candidates;
        candidates.push_back(base + ".stat.gz");
        candidates.push_back(base + ".stat");
        for (size_t i = 0; i < candidates.size(); i++) {
            if (fexists(candidates[i].c_str())) {
                args.stat_infile = strdup(candidates[i].c_str());
                break;
            }
        }
    }

    if (args.rlens_infile == NULL) {
        std::vector<std::string> candidates;
        candidates.push_back(base + ".rlens.gz");
        candidates.push_back(base + ".rlens");
        for (size_t i = 0; i < candidates.size(); i++) {
            if (fexists(candidates[i].c_str())) {
                args.rlens_infile = strdup(candidates[i].c_str());
                break;
            }
        }
    }

    return 0;
}

static int load_bam_lookup(const char *bam_file, std::map<std::string, int> &lookup) {
    if (bam_file == NULL)
        return 0;

    samFile *samfp = sam_open_format(bam_file, "r", NULL);
    if (samfp == NULL) {
        fprintf(stderr, "\t-> Error: could not open BAM/SAM/CRAM for resolve: %s\n", bam_file);
        return 1;
    }
    bam_hdr_t *hdr = sam_hdr_read(samfp);
    if (hdr == NULL) {
        fprintf(stderr, "\t-> Error: could not read header for resolve: %s\n", bam_file);
        sam_close(samfp);
        return 1;
    }

    for (int i = 0; i < hdr->n_targets; i++)
        lookup[std::string(hdr->target_name[i])] = i;

    bam_hdr_destroy(hdr);
    sam_close(samfp);
    return 0;
}

static int load_acc2tax_lookup(const char *acc2tax_file, std::map<std::string, int> &lookup) {
    if (acc2tax_file == NULL)
        return 0;

    gzFile fp = gzopen(acc2tax_file, "rb");
    if (fp == Z_NULL) {
        fprintf(stderr, "\t-> Error: could not open acc2tax for resolve: %s\n", acc2tax_file);
        return 1;
    }

    char line[8192];
    int at = 0;
    while (gzgets(fp, line, sizeof(line)) != NULL) {
        at++;
        if (at == 1)
            continue;

        char *tok = strtok(line, "\t\n ");
        if (tok == NULL)
            continue;
        char *key = strtok(NULL, "\t\n ");
        if (key == NULL)
            continue;
        char *val = strtok(NULL, "\t\n ");
        if (val == NULL)
            continue;

        lookup[std::string(key)] = atoi(val);
    }

    gzclose(fp);
    return 0;
}

static int resolve_selectors(filter_bdamage_args &args) {
    if (args.selectors.empty())
        return 0;

    std::map<std::string, int> lookup;
    if (load_bam_lookup(args.bam_file, lookup) != 0)
        return 1;
    if (load_acc2tax_lookup(args.acc2tax_file, lookup) != 0)
        return 1;

    int unresolved = 0;
    for (size_t i = 0; i < args.selectors.size(); i++) {
        const std::string &sel = args.selectors[i];
        std::map<std::string, int>::iterator it = lookup.find(sel);
        if (it == lookup.end()) {
            unresolved++;
            fprintf(stderr, "\t-> resolve: selector not found: %s\n", sel.c_str());
            continue;
        }
        args.ids.insert(it->second);
    }

    if (unresolved > 0 && args.strict_resolve == 1) {
        fprintf(stderr, "\t-> Error: %d selector(s) could not be resolved\n", unresolved);
        return 1;
    }

    if (unresolved > 0)
        fprintf(stderr, "\t-> resolve: ignored unresolved selectors because --strict_resolve 0\n");

    return 0;
}

static void expand_descendants_from_seed(int seed, int2intvec &child, std::set<int> &ids) {
    std::vector<int> stack;
    stack.push_back(seed);

    while (!stack.empty()) {
        int cur = stack.back();
        stack.pop_back();

        int2intvec::iterator it = child.find(cur);
        if (it == child.end())
            continue;

        for (size_t i = 0; i < it->second.size(); i++) {
            const int kid = it->second[i];
            if (ids.insert(kid).second)
                stack.push_back(kid);
        }
    }
}

static int maybe_expand_ids_from_nodes(filter_bdamage_args &args) {
    if (args.nodes_file == NULL || args.ids.empty())
        return 0;

    int2char rank;
    int2int parent;
    int2intvec child;
    parse_nodes(args.nodes_file, rank, parent, child, 1);

    const size_t before = args.ids.size();
    std::vector<int> seeds(args.ids.begin(), args.ids.end());
    for (size_t i = 0; i < seeds.size(); i++)
        expand_descendants_from_seed(seeds[i], child, args.ids);

    for (int2char::iterator it = rank.begin(); it != rank.end(); it++)
        free(it->second);

    fprintf(stderr,
            "\t-> nodes expansion with %s: %zu seed ids -> %zu ids after descendants\n",
            args.nodes_file,
            before,
            args.ids.size());
    return 0;
}

static int keep_entry(int id, int nreads, const filter_bdamage_args &args) {
    int keep = 1;

    if (!args.ids.empty()) {
        const int in_set = args.ids.find(id) != args.ids.end();
        keep = args.exclude ? !in_set : in_set;
    }

    if (keep && args.min_reads != -1 && nreads < args.min_reads)
        keep = 0;
    if (keep && args.max_reads != -1 && nreads > args.max_reads)
        keep = 0;

    return keep;
}

static int filter_rlens_companion(const char *infile, const char *outfile, const std::set<int> &kept_ids) {
    if (infile == NULL || outfile == NULL)
        return 0;

    gzFile in = gzopen(infile, "rb");
    if (in == Z_NULL) {
        fprintf(stderr, "\t-> Error: could not open rlens companion: %s\n", infile);
        return 1;
    }
    gzFile out = gzopen(outfile, "wb");
    if (out == Z_NULL) {
        fprintf(stderr, "\t-> Error: could not open rlens output: %s\n", outfile);
        gzclose(in);
        return 1;
    }

    char line[8192];
    size_t at = 0;
    size_t kept = 0;
    while (gzgets(in, line, sizeof(line)) != NULL) {
        at++;
        if (at == 1) {
            gzputs(out, line);
            continue;
        }

        char copy[8192];
        strncpy(copy, line, sizeof(copy));
        copy[sizeof(copy) - 1] = '\0';
        char *tok = strtok(copy, "\t\n ");
        if (tok == NULL)
            continue;

        int id = 0;
        if (parse_int(tok, &id) != 0)
            continue;

        if (kept_ids.find(id) != kept_ids.end()) {
            gzputs(out, line);
            kept++;
        }
    }

    gzclose(in);
    gzclose(out);

    fprintf(stderr, "\t-> Companion rlens filtered: %s -> %s (kept lines: %zu)\n", infile, outfile, kept);
    return 0;
}

static int resolve_stat_id(const char *tok,
                           const std::set<int> &kept_ids,
                           const std::map<std::string, int> &bam_lookup,
                           int has_bam_lookup,
                           int &is_keep,
                           int &is_resolved) {
    is_keep = 0;
    is_resolved = 0;

    if (tok == NULL)
        return 0;

    int id = 0;
    if (parse_int(tok, &id) == 0) {
        is_resolved = 1;
        is_keep = kept_ids.find(id) != kept_ids.end();
        return 0;
    }

    if (strcmp(tok, "global") == 0) {
        is_resolved = 1;
        is_keep = kept_ids.find(0) != kept_ids.end();
        return 0;
    }

    if (has_bam_lookup) {
        std::map<std::string, int>::const_iterator it = bam_lookup.find(tok);
        if (it != bam_lookup.end()) {
            is_resolved = 1;
            is_keep = kept_ids.find(it->second) != kept_ids.end();
            return 0;
        }
    }

    return 0;
}

static int filter_stat_companion(const char *infile,
                                 const char *outfile,
                                 const std::set<int> &kept_ids,
                                 const std::map<std::string, int> &bam_lookup,
                                 int has_bam_lookup) {
    if (infile == NULL || outfile == NULL)
        return 0;

    gzFile in = gzopen(infile, "rb");
    if (in == Z_NULL) {
        fprintf(stderr, "\t-> Error: could not open stat companion: %s\n", infile);
        return 1;
    }
    gzFile out = gzopen(outfile, "wb");
    if (out == Z_NULL) {
        fprintf(stderr, "\t-> Error: could not open stat output: %s\n", outfile);
        gzclose(in);
        return 1;
    }

    char line[8192];
    size_t at = 0;
    size_t kept = 0;
    size_t unresolved = 0;

    while (gzgets(in, line, sizeof(line)) != NULL) {
        at++;
        if (at == 1) {
            gzputs(out, line);
            continue;
        }

        char copy[8192];
        strncpy(copy, line, sizeof(copy));
        copy[sizeof(copy) - 1] = '\0';
        char *tok = strtok(copy, "\t\n ");
        if (tok == NULL)
            continue;

        int is_keep = 0;
        int is_resolved = 0;
        resolve_stat_id(tok, kept_ids, bam_lookup, has_bam_lookup, is_keep, is_resolved);
        if (!is_resolved) {
            unresolved++;
            continue;
        }
        if (is_keep) {
            gzputs(out, line);
            kept++;
        }
    }

    gzclose(in);
    gzclose(out);

    fprintf(stderr, "\t-> Companion stat filtered: %s -> %s (kept lines: %zu)\n", infile, outfile, kept);
    if (unresolved > 0 && !has_bam_lookup)
        fprintf(stderr, "\t-> Companion stat: %zu unresolved id labels (tip: provide --bam for refname->tid mapping)\n", unresolved);
    else if (unresolved > 0)
        fprintf(stderr, "\t-> Companion stat: %zu unresolved lines were skipped\n", unresolved);

    return 0;
}

int main_filter_bdamage(int argc, char **argv) {
    if (argc == 1 || (argc == 2 && (!strcasecmp(argv[1], "-h") || !strcasecmp(argv[1], "--help"))))
        return usage_filter_bdamage(stderr);

    filter_bdamage_args args = init_filter_bdamage_args();
    const int parse_rc = parse_filter_bdamage_args(argc, argv, args);
    if (parse_rc == 2)
        return 0;
    if (parse_rc != 0)
        return 1;

    if (maybe_infer_companion_inputs(args) != 0)
        return 1;
    if (resolve_selectors(args) != 0)
        return 1;
    if (maybe_expand_ids_from_nodes(args) != 0)
        return 1;

    char out_bdamage[4096];
    char out_stat[4096];
    char out_rlens[4096];
    snprintf(out_bdamage, sizeof(out_bdamage), "%s.bdamage.gz", args.out_prefix);
    snprintf(out_stat, sizeof(out_stat), "%s.stat.gz", args.out_prefix);
    snprintf(out_rlens, sizeof(out_rlens), "%s.rlens.gz", args.out_prefix);

    fprintf(stderr,
            "\t-> filter_bdamage infile: %s out_prefix: %s ids: %zu selectors: %zu exclude: %d min_reads: %d max_reads: %d nthreads: %d\n",
            args.infile,
            args.out_prefix,
            args.ids.size(),
            args.selectors.size(),
            args.exclude,
            args.min_reads,
            args.max_reads,
            args.nthreads);
    fprintf(stderr,
            "\t-> companions stat:%s rlens:%s\n",
            args.stat_infile ? args.stat_infile : "(none)",
            args.rlens_infile ? args.rlens_infile : "(none)");

    BGZF *in = bgzf_open(args.infile, "r");
    if (in == NULL) {
        fprintf(stderr, "\t-> Error: failed to open input file: %s\n", args.infile);
        return 1;
    }
    BGZF *out = bgzf_open(out_bdamage, "w");
    if (out == NULL) {
        fprintf(stderr, "\t-> Error: failed to open output file: %s\n", out_bdamage);
        bgzf_close(in);
        return 1;
    }

    if (args.nthreads > 1) {
        bgzf_mt(in, args.nthreads, 256);
        bgzf_mt(out, args.nthreads, 256);
    }

    int printlength = 0;
    if (bgzf_read(in, &printlength, sizeof(int)) != sizeof(int)) {
        fprintf(stderr, "\t-> Error: failed to read printlength from input bdamage file\n");
        bgzf_close(in);
        bgzf_close(out);
        return 1;
    }
    if (printlength <= 0) {
        fprintf(stderr, "\t-> Error: invalid printlength in input bdamage file: %d\n", printlength);
        bgzf_close(in);
        bgzf_close(out);
        return 1;
    }

    if (bgzf_write(out, &printlength, sizeof(int)) != sizeof(int)) {
        fprintf(stderr, "\t-> Error: failed to write printlength to output bdamage file\n");
        bgzf_close(in);
        bgzf_close(out);
        return 1;
    }

    const size_t nvals = (size_t)printlength * 16 * 2;
    const size_t datab = nvals * sizeof(float);
    std::vector<float> data(nvals);

    std::set<int> kept_ids;
    size_t total = 0;
    size_t kept = 0;
    while (1) {
        int ref_nreads[2];
        const int nread = bgzf_read(in, ref_nreads, 2 * sizeof(int));
        if (nread == 0)
            break;
        if (nread != 2 * (int)sizeof(int)) {
            fprintf(stderr, "\t-> Error: truncated/corrupt bdamage file (record header)\n");
            bgzf_close(in);
            bgzf_close(out);
            return 1;
        }

        if (bgzf_read(in, data.data(), datab) != (ssize_t)datab) {
            fprintf(stderr, "\t-> Error: truncated/corrupt bdamage file (record payload)\n");
            bgzf_close(in);
            bgzf_close(out);
            return 1;
        }

        total++;
        if (!keep_entry(ref_nreads[0], ref_nreads[1], args))
            continue;

        if (bgzf_write(out, ref_nreads, 2 * sizeof(int)) != 2 * sizeof(int)) {
            fprintf(stderr, "\t-> Error: failed writing record header to output bdamage file\n");
            bgzf_close(in);
            bgzf_close(out);
            return 1;
        }
        if (bgzf_write(out, data.data(), datab) != (ssize_t)datab) {
            fprintf(stderr, "\t-> Error: failed writing record payload to output bdamage file\n");
            bgzf_close(in);
            bgzf_close(out);
            return 1;
        }
        kept_ids.insert(ref_nreads[0]);
        kept++;
    }

    bgzf_close(in);
    bgzf_close(out);

    std::map<std::string, int> bam_lookup;
    int has_bam_lookup = 0;
    if (args.bam_file != NULL) {
        if (load_bam_lookup(args.bam_file, bam_lookup) != 0)
            return 1;
        has_bam_lookup = 1;
    }

    if (args.stat_infile != NULL) {
        if (filter_stat_companion(args.stat_infile, out_stat, kept_ids, bam_lookup, has_bam_lookup) != 0)
            return 1;
    }
    if (args.rlens_infile != NULL) {
        if (filter_rlens_companion(args.rlens_infile, out_rlens, kept_ids) != 0)
            return 1;
    }

    fprintf(stderr,
            "\t-> filter_bdamage done. Input entries: %zu Kept entries: %zu Removed entries: %zu Output: %s\n",
            total,
            kept,
            total - kept,
            out_bdamage);

    free(args.infile);
    free(args.out_prefix);
    if (args.id_file)
        free(args.id_file);
    if (args.resolve_file)
        free(args.resolve_file);
    if (args.bam_file)
        free(args.bam_file);
    if (args.acc2tax_file)
        free(args.acc2tax_file);
    if (args.nodes_file)
        free(args.nodes_file);
    if (args.stat_infile)
        free(args.stat_infile);
    if (args.rlens_infile)
        free(args.rlens_infile);

    return 0;
}
