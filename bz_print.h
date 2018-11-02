/* $Id: bz_print.h,v 4.0 2004/08/26 09:49:25 mmichael Exp mmichael $ */

#include "bz_main.h"

void print_job_footer(void);
void print_job_header (ss_t ss, gap_scores_t *ds, int K, int L, int M, int protein);
void print_align_header(SEQ *seq1, SEQ *seq2, ss_t ss, gap_scores_t *ds,
  int K, int L);

void print_align_list(align_t *a);
void free_align_list(align_t *a);
void print_msp_list(SEQ *sf1, SEQ *sf2, msp_table_t *mt);
