#include <zlib.h>
#include <htslib/sam.h>
#include "types.h"


typedef struct{
  bam1_t **ary;
  size_t l;
  size_t m;
}queue;

queue *init_queue(size_t maxsize);
void expand_queue(queue *ret);

int2char parse_names(const char *fname);
void parse_nodes(const char *fname,int2char &rank,int2int &parent);
void strip(char *line);
void parse_nodes(const char *fname,int2char &rank,int2int &parent,int2intvec &child,int dochild);
int2int *bamRefId2tax(bam_hdr_t *hdr,char *acc2taxfile,char *bamfile,int2int &errmap);
