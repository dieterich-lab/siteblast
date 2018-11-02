/* $Id: bz_table.h,v 4.0 2004/08/26 09:49:25 mmichael Exp mmichael $ */
#ifndef BZ_TABLE__H 
#define BZ_TABLE__H 
#include "int_types.h"
#include <assert.h>
#include "recognize.h"
#include "seq.h"
#include "dna.h"

#ifndef ZFREE
#define ZFREE(p) /*CONSTCOND*/do{ckfree(p);(p)=0;}while(0)
#endif

static const char *ori_names[] ={"sense","rev","compl.","rev. cmpl"};

enum Extratype { 
  UNDEF_EXTRA_INFO=0,
  IUPAC_EXTRA_INFO=1,
  WEIGHT_EXTRA_INFO=2
}; 

typedef struct seed {
  struct seed *next;
  /*  char *id, dist, ori;*/
  void* extra_information;
  enum Extratype extra_type;
  char len;
  int pos1, pos2;
} seed_t;

seed_t *cloneseed(seed_t *s, char recursiv);
seed_t *sort_seed(seed_t *s);
seed_t *rev_seed(seed_t *s);
void free_seeds(seed_t *seed);

typedef u32_t blast_ecode_t;

enum { BLAST_POS_TAB_SIZE = 65536UL };  /* 65536 = 4**8 */

typedef struct blast_table_w {
        int W;
        int num;
        const signed char *encoding;
        struct blast_epos_node **epos_tab;
} blast_table_w;

typedef struct blast_table_p { /*profile search*/
  const signed char *encoding;
  seed_t **occur1, **occur2;
  void *pssm;
  char orintations;
  int highest_pattern_id;
  char* id_list;
} blast_table_p;

typedef struct blast_table_i { /*Consensus seq search using Aho-Corasick*/
  const signed char *encoding;
  seed_t **occur1, **occur2;
  Recognizer r;
  int highest_pattern_id;
  char* id_list;
} blast_table_i;

typedef struct blast_table_I { /*Consensus seq search using primitiv search*/
const signed char *encoding;
  seed_t **occur1, **occur2;
  char **cons;
  char dist;
  char **ids;
  int type;
  int number;
  char orintations;
  int highest_pattern_id;
  char* id_list;
} blast_table_I;


typedef struct blast_table_n {
        struct hashpos **hashtab; /* array of pointers to positions in seq1, indexed by 12-mer */
        struct hashpos *pool;     /* pool of free nodes */
        int npool;      /* number of nodes available */
        int n_node;     /* index of next free node in pool */
} blast_table_n;

enum blast_table_type{ 
  blast_table_type_unknow=0,
  blast_table_type_w=1,
  blast_table_type_p=2,
  blast_table_type_I=3,
  blast_table_type_i=4,
  blast_table_type_n=5
}; 

typedef union {
  struct blast_table_w w;
  struct blast_table_p p;
  struct blast_table_I I;
  struct blast_table_i i;
  struct blast_table_n n;
} blast_table_t;

#ifndef SCORE_T
#define SCORE_T long
#endif
typedef SCORE_T score_t;

typedef struct {
	int len, pos1, pos2;
	score_t score, cum_score;
	int filter;
} msp_t;

typedef struct msp_table {
	int size;
	int num;
	msp_t *msp;
} msp_table_t;
#define MSP_TAB_FIRST(t) ((t)->msp)
#define MSP_TAB_NEXT(m) ((m)+1)
#define MSP_TAB_NUM(t) ((t)->num)
#define MSP_TAB_MORE(t,m) ((m-MSP_TAB_FIRST(t))<MSP_TAB_NUM(t))
#define MSP_TAB_NTH(t,n) ((t)->msp[n])

typedef int (*msp_cmp_t)(const void *, const void *);

blast_table_t *blast_table_new(SEQ *, int);
msp_table_t *blast_search(SEQ *seq1, SEQ *seq2, blast_table_t *bt, msp_table_t *mt, ss_t ss, int X, int K, int P, int T);
void blast_table_free(blast_table_t *bt);
void msp_free_table(msp_table_t *mt);
msp_table_t *msp_new_table(void);

msp_table_t *msp_compress(msp_table_t *mt);
msp_table_t *msp_sort_pos1(msp_table_t *mt);
msp_table_t *msp_sort(msp_table_t *mt);
int msp_add(msp_table_t *, int, int, int, score_t, int);

typedef int (*connect_t) (msp_t *, msp_t *, int);
int msp_make_chain (msp_table_t *, int, int, int, connect_t);

struct comb;
int msp_extend_hit(msp_table_t *mt,
                        SEQ *s1, SEQ *s2, ss_t ss, int X, int K, int W, int P,
                        int pos1, int pos2,
                        struct comb *diag_lev
                        );

typedef struct blast_epos_node {
	struct blast_epos_node *link;	/* next occurrence of the word */
	blast_ecode_t ecode;
	int pos;		/* position where word occurs in seq #1 */
} blast_epos_node_t;

#include "bz_all.h"
#endif




