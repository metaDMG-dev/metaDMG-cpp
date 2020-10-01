#pragma once
#include <vector>
#include <map>


#define bam_is_sec(b)         (((b)->core.flag&BAM_FSECONDARY)        != 0)
#define bam_is_supp(b)        (((b)->core.flag&BAM_FSUPPLEMENTARY)    != 0)
#define bam_is_paired(b)      (((b)->core.flag&BAM_FPAIRED)     != 0)
#define bam_is_rmdup(b)       (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_qcfailed(b)    (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_unmapped(b)    (((b)->core.flag&BAM_FUNMAP)      != 0)
#define bam_is_read1(b)       (((b)->core.flag&BAM_FREAD1)      != 0)
#define bam_is_failed(b)      ( bam_is_qcfailed(b) || bam_is_rmdup(b) || bam_is_supp(b) )

typedef struct{
  size_t nreads;
  unsigned **mm5p;
  unsigned **mm3p;
}triple;


class damage{
  char *reconstructedTemp;
public:
  int nthreads;
  int MAXLENGTH;
  int minQualBase;//currently not set; should be set in init_damage
  int nclass;
  unsigned **mm5p;//will point to first. Just a hack
  unsigned **mm3p;//will point to first. Just a hack
  std::pair< kstring_t*, std::vector<int> >  reconstructedReference;
  std::map<int,triple > assoc;
  void write(char *prefix,bam_hdr_t *hdr);
  void bwrite(char *prefix,bam_hdr_t *hdr);
  int damage_analysis( bam1_t *b,int whichclass);
  void printit(FILE *fp,int l);
  damage(int maxlen,int nthd,int minqb){
    MAXLENGTH = maxlen;
    minQualBase = minqb;
    nthreads = nthd;
    reconstructedTemp=(char*)calloc(256,1);
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    reconstructedReference.first = kstr;
    mm5p=mm3p=NULL;
  }
  ~damage(){free(reconstructedTemp);}
};

void destroy_damage(damage *dmg);

std::map<int,double *> load_bdamage(const char *fname,int howmany);
std::map<int,double *> load_bdamage2(const char *fname,int howmany);
