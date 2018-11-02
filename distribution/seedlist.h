#ifndef SEEDLIST__H
#define SEEDLIST__H
#include "bz_table.h"

typedef struct {
  int len;
  int offset;
  int *used;
  int nused;
  seed_t **seed, **lastseed;
} seedList_t;

void free_seedList();
void init_seedList(int len1, int len2);
seed_t *get_seed_from_list(int pos);
void insert_seed_into_list(int pos, seed_t *seed);
void print_seedList();
#endif
