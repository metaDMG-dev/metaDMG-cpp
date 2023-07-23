#include <htslib/kstring.h>
#include "types.h"
void print_chain(kstring_t *kstr, int taxa, int2int &parent, int2char &rank, int2char &name_map);
float gccontent(char *seq);
char *make_seq(bam1_t *aln);
float gccontent(bam1_t *aln);
float mean(std::vector<float> &vec);
float var(std::vector<float> &vec);
