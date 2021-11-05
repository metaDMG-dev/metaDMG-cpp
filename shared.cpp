#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <libgen.h>
#include <htslib/bgzf.h>
#include <cassert>
#include <ctime>
#include "shared.h"


BGZF *getbgzf(const char*str1,const char *mode,int nthreads){
  BGZF *fp = NULL;
  fp = bgzf_open(str1,mode);
  fprintf(stderr,"\t-> opening file: \'%s\' mode: \'%s\'\n",str1,mode);
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem opening file: \"%s\"\n",str1);
    exit(0);
  }
  if(nthreads>1){
    fprintf(stderr,"\t-> Setting threads to: %d \n",nthreads);
    bgzf_mt(fp,nthreads,64);
  }
  return fp;
}

BGZF *getbgzf2(const char*str1,const char *str2,const char *mode,int nthreads){
  unsigned tmp_l = strlen(str1)+strlen(str2)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s",str1,str2);
  return  getbgzf(tmp,mode,nthreads);
}

BGZF *getbgzf3(const char*str1,const char *str2,const char *str3,const char *mode,int nthreads){
  unsigned tmp_l = strlen(str1)+strlen(str2)+strlen(str3)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s%s",str1,str2,str3);
  return  getbgzf(tmp,mode,nthreads);
}

BGZF *getbgzf4(const char*str1,const char *str2,const char *str3,const char *str4,const char *mode,int nthreads){
  unsigned tmp_l = strlen(str1)+strlen(str2)+strlen(str3)+strlen(str4)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s%s%s",str1,str2,str3,str4);
  return  getbgzf(tmp,mode,nthreads);
}



int fexists(const char* str){///@param str Filename given as a string.
  fprintf(stderr,"\t-> Checking if exits: \'%s\'\n",str);
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

int fexists2(const char*str1,const char* str2){
  unsigned tmp_l = strlen(str1)+strlen(str2)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s",str1,str2);
  return fexists(tmp);
}

size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


int fexists3(const char*str1,const char* str2,const char *str3){
  unsigned tmp_l = strlen(str1)+strlen(str2)+strlen(str3)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s%s",str1,str2,str3);
  size_t fs=fsize(tmp);
  return fexists(tmp)&&fs>0;
}




int fexists4(const char*str1,const char* str2,const char *str3,const char *str4){
  unsigned tmp_l = strlen(str1)+strlen(str2)+strlen(str3)+strlen(str4)+5;
  char tmp[tmp_l];
  snprintf(tmp,tmp_l,"%s%s%s%s",str1,str2,str3,str4);
  //  fprintf(stderr,"\t-> checking if : %s exists\n",tmp );
  size_t fs=fsize(tmp);
  return fexists(tmp)&&fs>0;
}



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
  gzclose(gz);
  delete [] toks;
  return name_map;
}



void parse_nodes(const char *fname,int2char &rank,int2int &parent,int2intvec &child,int dochild){
  //  fprintf(stderr,"Parsing: %s\n",fname);
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
      if(dochild){
	if(key==val)//catch 1 <-> 1
	  continue;
	int2intvec::iterator it2 = child.find(val);
	if(it2==child.end()){
	  std::vector<int> tmp;
	  tmp.push_back(key);
	  child[val] = tmp;
	}else
	  it2->second.push_back(key);
      }
    }
  }
  fprintf(stderr,"\t-> Number of unique names (column1): %lu from file: %s parent.size():%lu child.size():%lu\n",rank.size(),fname,parent.size(),child.size());
  //int2int::iterator it=parent.find(9532);
  //fprintf(stderr,"%d->%d\n",it->first,it->second);
  //int2intvec::iterator it=child.find(1);
  //fprintf(stderr,"%d->%lu\n",it->first,it->second.size());
  //exit(0);
  gzclose(gz);
  delete [] toks;
}

//this generates a downtree, parent to childs taxid->vector<taxids>
void parse_nodes2(int2int &parent,int2intvec &child){
  fprintf(stderr,"\t-> Generating reverse node table\n");
  for(int2int::iterator it=parent.begin();it!=parent.end();it++){
    int down=it->first;
    int up = it->second;
    //    fprintf(stdout,"%d\t%d\n",down,up);
    int2intvec::iterator it2=child.find(up);
    if(it2==child.end()){
      std::vector<int> tmp;
      tmp.push_back(down);
      child[up] = tmp;
    }else
      it2->second.push_back(down);

  }
  fprintf(stderr,"\t-> Done generating reverse node table: contains: %lu\n",child.size());
  //    exit(0);
}

//bamfile is only used for making smaller filedump (using the filename)
int SIG_COND=1;
int2int *bamRefId2tax(bam_hdr_t *hdr,char *acc2taxfile,char *bamfile,int2int &errmap,char *tempfolder) { 
  fprintf(stderr,"\t-> Starting to extract (acc->taxid) from binary file: \'%s\'\n",acc2taxfile);
  fflush(stderr);
  int dodump = !fexists4(tempfolder,basename(acc2taxfile),basename(bamfile),".bin");
  
  fprintf(stderr,"\t-> Checking if bimnary file exists. dodump=%d \n", dodump);
  
  time_t t=time(NULL);
  BGZF *fp= NULL;
  if(dodump)
    fp = getbgzf4(tempfolder,basename(acc2taxfile),basename(bamfile),".bin","wb",4);
  else
    fp =  getbgzf4(tempfolder,basename(acc2taxfile),basename(bamfile),".bin","rb",4);
  //this contains refname(as int) -> taxid
  int2int *am= new int2int;
  
  
  if(dodump){
    
    char buf[4096];
    int at=0;
    char buf2[4096];
    
    kstring_t *kstr =(kstring_t *)malloc(sizeof(kstring_t));
    kstr->l=kstr->m = 0;
    kstr->s = NULL;
    BGZF *fp2 = getbgzf(acc2taxfile,"rb",2);
    bgzf_getline(fp2,'\n',kstr);//skip header
    kstr->l =0;
    while(SIG_COND&&bgzf_getline(fp2,'\n',kstr)){
      if(kstr->l==0)
	break;
      //fprintf(stderr,"at: %d = '%s\'\n",at,kstr->s);
      if(!((at++ %100000 ) ))
	if(isatty(fileno(stderr)))
	  fprintf(stderr,"\r\t-> At linenr: %d in \'%s\'      ",at,acc2taxfile);
      char *tok = strtok(kstr->s,"\t\n ");
      char *key =strtok(NULL,"\t\n ");
      tok = strtok(NULL,"\t\n ");
      int val = atoi(tok);
      //fprintf(stderr,"key: %d val: %d\n",key,val);exit(0);
      int valinbam = bam_name2id(hdr,key);
      if(valinbam==-1)
	continue;
      assert(bgzf_write(fp,&valinbam,sizeof(int))==sizeof(int));
      assert(bgzf_write(fp,&val,sizeof(int))==sizeof(int));
      //fprintf(stderr,"key: %s val: %d valinbam:%d\n",key,val,valinbam);
      
      if(am->find(valinbam)!=am->end())
	fprintf(stderr,"\t-> Duplicate entries found \'%s\'\n",key);
      (*am)[valinbam] = val;
      kstr->l =0;
    }
    bgzf_close(fp2);
  }else{
    int valinbam,val;
    while(bgzf_read(fp,&valinbam,sizeof(int))){
      assert(bgzf_read(fp,&val,sizeof(int))==sizeof(int));
      (*am)[valinbam] = val;
    }
  }

  bgzf_close(fp);
  fprintf(stderr,"\t-> Number of entries to use from accesion to taxid: %lu, time taken: %.2f sec\n",am->size(),(float)(time(NULL) - t));
  return am;
}


queue *init_queue(size_t maxsize){
  queue *ret = new queue;
  ret->l =0;
  ret->m =maxsize;
  ret->ary = new bam1_t*[ret->m];
  for(int i=0;i<ret->m;i++)
    ret->ary[i] = bam_init1();
  return ret;
}

//expand queue with 500 elements
void expand_queue(queue *ret){
  bam1_t **newary = new bam1_t*[ret->m+500];
  for(int i=0;i<ret->l;i++)
    newary[i] = ret->ary[i];
  for(int i=ret->l;i<ret->m+500;i++)
    newary[i] = bam_init1();
  delete [] ret->ary;
  ret->ary = newary;
  ret->m += 500;
}

void destroy_queue(queue *q){
  for(int i=0;i<q->m;i++)
    bam_destroy1(q->ary[i]);
  delete [] q->ary;
  delete q;
  q = NULL;
}
