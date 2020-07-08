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
  int MAXLENGTH;
  int minQualBase;//currently not set; should be set in init_damage
  int nclass;
  unsigned **mm5p;//will point to first. Just a hack
  unsigned **mm3p;//will point to first. Just a hack
  std::pair< kstring_t*, std::vector<int> >  reconstructedReference;
  std::map<int,triple > assoc;
  void write(char *prefix,bam_hdr_t *hdr);
  int damage_analysis( bam1_t *b,int whichclass);
  void printit(FILE *fp,int l);
  damage(){
    reconstructedTemp=(char*)calloc(256,1);
  }
};

damage *init_damage(int MAXLENGTH);
void destroy_damage(damage *dmg);
