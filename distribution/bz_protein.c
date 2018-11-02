#include "bz_all.h"

/* Protein characters */
#define PROTEIN_CHARS "ACDEFGHIKLMNPQRSTVWXY"

bool is_Protein(uchar *s, int len){
  int i;
  for (i = 0; i < len; ++i){
    if (!strchr(PROTEIN_CHARS, toupper(s[i])))
      return 0;
  }
  return 1;
}
