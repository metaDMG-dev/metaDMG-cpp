// gpl thorfinn@binf.ku.dk
#pragma once

#include <map>

#include "profile.h"
#include "types.h"

double *getval(std::map<int, double *> &retmap, int2intvec &child, int taxid, int howmany);
std::map<int, mydataD> getval_full_norec(std::map<int, mydataD> &retmap, int2int &parent, int howmany);
mydata2 getval_stats(std::map<int, mydata2> &retmap, int2intvec &child, int taxid);

int main_print(int argc, char **argv);
int main_print2(int argc, char **argv);
int main_print_all(int argc, char **argv);
int main_print_ugly(int argc, char **argv);
