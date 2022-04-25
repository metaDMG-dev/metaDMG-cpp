#include "types.h"
void print_chain(gzFile fp, int taxa, int2int &parent, int2char &rank, int2char &name_map);
float gccontent(char *seq);
char *make_seq(bam1_t *aln);
float gccontent(bam1_t *aln);
float mean(std::vector<float> &vec);
float var(std::vector<float> &vec);
