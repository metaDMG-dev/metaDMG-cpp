#include "shared.h"


//usefull little function to split
char *strpop(char **str,char split){
  char *tok=*str;
  while(**str){
    if(**str!=split)
      (*str)++;
    else{
      **str='\0'; (*str)++;
      break;
    }
  }
  return tok;
}
//usefull little function to remove tab and newlines
void strip(char *line){
  int at=0;
  //  fprintf(stderr,"%s\n",line);
  for(int i=0;i<strlen(line);i++)
    if(line[i]=='\t'||line[i]=='\n')
      continue;
    else
      line[at++]=line[i];
  //  fprintf(stderr,"at:%d line:%p\n",at,line);
  line[at]='\0';
  //fprintf(stderr,"%s\n",line);
}


int2char parse_names(const char *fname){

  gzFile gz= Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
    exit(0);
  }
  int2char name_map;
  char buf[4096];
  int at=0;
  char **toks = new char*[5];
  
  while(gzgets(gz,buf,4096)){
    strip(buf);//fprintf(stderr,"buf:%s\n",buf);
    char *saveptr = buf;
    toks[0]=strpop(&saveptr,'|');
    toks[1]= strpop(&saveptr,'|');
    toks[2]= strpop(&saveptr,'|');
    toks[3]= strpop(&saveptr,'|');
    for(int i=0;0&&i<4;i++)
      fprintf(stderr,"%d):\'%s\'\n",i,toks[i]);

    int key=atoi(toks[0]);
    //    fprintf(stderr,"key:%d\n",key);
    if(toks[3]&&strcmp(toks[3],"scientific name")==0){
      int2char::iterator it=name_map.find(key);
      
      if(it!=name_map.end())
	fprintf(stderr,"\t->[%s] duplicate name(column1): %s\n",fname,toks[0]);
      else
	name_map[key]=strdup(toks[1]);

    }
    if(0&&at++>10)
      break;
  }
  //  int2char::iterator it = name_map.find(61564);  assert(it!=name_map.end());
  fprintf(stderr,"\t-> [%s] Number of unique names (column1): %lu with third column 'scientific name'\n",fname,name_map.size());
  return name_map;
}



void parse_nodes(const char *fname,int2char &rank,int2int &parent){

  gzFile gz= Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fname);
    exit(0);
  }
  char buf[4096];
  int at=0;
  char **toks = new char*[5];

  while(gzgets(gz,buf,4096)){
    strip(buf);//fprintf(stderr,"buf:%s\n",buf);
    char *saveptr = buf;
    toks[0]= strpop(&saveptr,'|');
    toks[1]= strpop(&saveptr,'|');
    toks[2]= strpop(&saveptr,'|');
    for(int i=0;0&&i<3;i++)
      fprintf(stderr,"%d):\'%s\'\n",i,toks[i]);

    int2int::iterator it = parent.find(atoi(toks[0]));
    if(it!=parent.end())
      fprintf(stderr,"\t->[%s] duplicate name(column0): %s\n",fname,toks[0]);
    else{
      int key=atoi(toks[0]);
      int val= atoi(toks[1]);
      parent[key]=val;
      rank[key]=strdup(toks[2]);
    }
  }
  fprintf(stderr,"\t-> Number of unique names (column1): %lu from file: %s\n",rank.size(),fname);

}

//this generates a downtree, parent to childs taxid->vector<taxids>
void parse_nodes2(int2int &parent,int2intvec &child){
  for(int2int::iterator it=parent.begin();it!=parent.end();it++){
    int down=it->first;
    int up = it->second;
    int2intvec::iterator it2=child.find(up);
    if(it2==child.end()){
      std::vector<int> tmp;
      tmp.push_back(down);
      child[up] = tmp;
    }else
      it2->second.push_back(down);

  }

}
