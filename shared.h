#include <htslib/sam.h>  // for bam1_t, bam_hdr_t
#include <stddef.h>      // for size_t

#include "types.h"  // for int2int, int2char, int2intvec
typedef struct {
    bam1_t **ary;
    size_t l;
    size_t m;
} queue;

queue *init_queue(size_t maxsize);
void expand_queue(queue *ret);
void destroy_queue(queue *q);

int2char parse_names(const char *fname);
void parse_nodes(const char *fname, int2char &rank, int2int &parent);
void strip(char *line);
void parse_nodes(const char *fname, int2char &rank, int2int &parent, int2intvec &child, int dochild);
int2int *bamRefId2tax(bam_hdr_t *hdr, char *acc2taxfile, char *bamfile, char *tempfolder, int reallyDump, char *filteredAcc2taxfile,char2int *acc2taxmap);

int fexists(const char *str);
