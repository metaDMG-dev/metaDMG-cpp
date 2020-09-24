#pragma once
#include <cstring>
#include <vector>
#include <map>
struct cmp_str
{
   bool operator()(char const *a, char const *b) const
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef struct{
  int up;
  std::vector<int> down;
  int dist2root;
}node;

typedef std::map<int,char *> int2char;
typedef std::map<int,int> int2int;
typedef std::map<int,node> int2node;
typedef std::map<char *,int,cmp_str> char2int;

