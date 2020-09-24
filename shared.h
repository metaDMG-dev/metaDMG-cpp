#include <zlib.h>
#include "types.h"

int2char parse_names(const char *fname);
void parse_nodes(const char *fname,int2char &rank,int2int &parent);
void strip(char *line);
void parse_nodes2(int2int &parent,int2intvec &child);
