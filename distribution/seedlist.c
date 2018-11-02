#include "seedlist.h"
#include <assert.h>

static seedList_t seedList;
static bool noSeeds=1;

void free_seedList(){
  int i;
  assert(!noSeeds);
  for(i=0;i<seedList.nused;i++){
    free_seeds(seedList.seed[seedList.used[i]]);
  }
  ckfree(seedList.seed);
  ckfree(seedList.lastseed);
  ckfree(seedList.used);
  seedList.len=-1;
  noSeeds=1;
}

void init_seedList(int len1, int len2){
  const int lowest=1-len2;
  const int highest=len1-1;

  seedList.len=highest-lowest+1;
  seedList.offset=-lowest;
  seedList.used=ckcalloc(seedList.len,sizeof(int));
  seedList.nused=0;
  seedList.seed=ckcalloc(seedList.len,sizeof(seed_t*));
  seedList.lastseed=ckcalloc(seedList.len,sizeof(seed_t*));
  noSeeds=0;
}

seed_t *get_seed_from_list(int pos){
  if(noSeeds) return NULL;
  pos+=seedList.offset;
  assert(pos>=0&&pos<seedList.len);
  return seedList.seed[pos];
}

void insert_seed_into_list(int pos, seed_t *seed){
  assert(!noSeeds);
  pos+=seedList.offset;
  assert(pos>=0&&pos<seedList.len);
  if(! seedList.seed[pos]){
    seedList.used[seedList.nused++] = pos;
    seedList.seed[pos]=seedList.lastseed[pos]=seed;
  }
  else{
    if(! seedList.lastseed[pos]) 
      fatalf("Internal Error. Exiting.\t\t(%s:%d)\n",__FILE__,__LINE__);
    seedList.lastseed[pos]->next=seed;
  }
  for(;seedList.lastseed[pos]->next;seedList.lastseed[pos]=seedList.lastseed[pos]->next);
}

void print_seedList(){
  int i;
  for(i=0;i<seedList.nused;i++){
    printf("seedList[%d]:\n",seedList.used[i]-seedList.offset);
    print_seeds(seedList.seed[seedList.used[i]]);
  }
}

