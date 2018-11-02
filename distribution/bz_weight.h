#ifndef BZ_WEIGHT__H
#define BZ_WEIGHT__H
#ifndef BZ_WEIGHT__C
extern char alphabet[];
#endif

Recognizer new_weightrecogniser(bz_flags_t *bz_flags, int seq_length, int *highest_pattern_id, char** id_list);

typedef int row[4];
typedef double frow[4];

typedef struct PSSM{
  char *id;
  int length;
  int *max_remain_score1, *max_remain_score2;
  int *rev_max_remain_score1, *rev_max_remain_score2;
  double granularity;
  char gran_width;
  double t1, t2; /*threshold*/
  row *rows1, *rows2;
  char *order1, *order2;
  struct PSSM *next;
} pssm_t;

typedef struct SinglePSSM{
  int length;
  int *max_remain_score;
  int *rev_max_remain_score;
  double granularity;
  double t; /*threshold*/
  row *rows;
  char *order;
} singlepssm_t;

typedef struct Profile{
  char *id;
  int length;
  frow* rows;
  struct Profile *next;
} profile_t;


typedef struct {
  char ori;
  int score1,score2;
  char *id;
  double granularity;
  char gran_width;
} weight_information;

weight_information *clone_weight_extra_info(weight_information *info);
pssm_t *get_PSSM(bz_flags_t *bz_flags, char** id_list, int seq_length, int *highest_pattern_id);
#endif








