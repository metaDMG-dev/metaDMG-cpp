#include <cstdio>
#include <zlib.h>
#include <map>
#include <cstring>
struct cmp_str
{
   bool operator()(char const *a, char const *b) const
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef std::map<char *, int, cmp_str> amap;

int main(int argc,char**argv){
  amap cmap;
  if(argc!=3){
    fprintf(stderr,"./extract_reads contigfile file.sam.gz\n");
    return 0;
  }
  char *key = argv[1];
  char *hts = argv[2];
  fprintf(stderr,"\t-> key: %s hts: %s\n",key,hts);
  FILE *fp = NULL;

  if(((fp=fopen(key,"rb")))==NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",key);
    return 0;
  }
  char buf[4096];
  while(fgets(buf,4096,fp)){
    char *tok = strtok(buf,"\n\t ");
    if(cmap.find(tok)!=cmap.end()){
      fprintf(stderr,"\t-> key: %s already exist will skip\n",buf);
    }
    cmap[strdup(tok)] = 1;
  }
  fprintf(stderr,"\t-> Done reading keys from: \'%s\' nitems: %lu\n",key,cmap.size());
  fclose(fp);

  gzFile gz = Z_NULL;
  if(((gz=gzopen(hts,"rb")))==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",hts);
    return 0;
  }
  char orig[4096];
  while(gzgets(gz,buf,4096)){
    if(buf[0]=='@'){//header
      fprintf(stdout,"%s",buf);
      continue;
    }
    //fprintf(stderr,"Done reading header\n");
    //noheader
    char *orig = strncpy(orig,buf,4096);
    char *tok = strtok(buf,"\t\n ");
    tok = strtok(NULL,"\t\n ");
    tok = strtok(NULL,"\t\n ");//this is now chr/contig in sam
    if(cmap.find(tok)!=cmap.end()){
      fprintf(stdout,"%s",orig);
      continue;
    }
  }
  gzclose(gz);
  return 0;
}
