#ifndef BZ_IUPAC__H
#define BZ_IUPAC__H

Recognizer new_iupacrecogniser(bz_flags_t *bz_flags, int *highest_pattern_id, char** id_list);
void getCons(char ***Cons, char ***Ids, int *number,bz_flags_t *bz_flags, int *highest_pattern_id, char** id_list);

typedef struct iupac_information{
  char *id;
  char dist1,dist2;
  char ori;
} iupac_information;

iupac_information *clone_iupac_extra_info(iupac_information *info);

#endif

