#include "bz_all.h"

#ifndef __lint
static const char rcsid[] =
"$Id: bz_print.c,v 4.0 2004/08/26 09:49:25 mmichael Exp mmichael $";
#endif

static int isrev(SEQ *s)
{
	return (s->flags & SEQ_IS_REVCOMP);
}

static const char* revflag(SEQ *s)
{
	return (s->flags & SEQ_IS_REVCOMP) ? "-" : "";
}

static const char* revlabel(SEQ *s)
{
	return (s->flags & SEQ_IS_REVCOMP) ? " (reverse complement)" : "";
}

static int align_match_percent(int run, int match)
{
	return (run > 0) ? ((100*match + run/2)/run) : 0; /* round */
}

void print_job_footer()
{
	printf("#:eof\n");
}

void print_job_header (ss_t ss, gap_scores_t *ds, int K, int L, int M, int protein)
{
	int i, j;

	printf("#:lav\nd {\n  \"");
	ck_argc("print_align_header");
	fprintf_argv(stdout);
	
	if(protein){
	  uchar b[] = "ACDEFGHIKLMNPQRSTVWY";
	  printf("\n     A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y\n ");
	  for (i = 0; i < 20; ++i) {
	    for (j = 0; j < 20; ++j)
	      printf("%5d", ss[b[i]][b[j]]);
	    printf("\n ");
	  }
	}
	else{
	  uchar b[] = "ACGT";
	  printf("\n     A    C    G    T\n ");
	  for (i = 0; i < 4; ++i) {
	    for (j = 0; j < 4; ++j)
	      printf("%5d", ss[b[i]][b[j]]);
	    printf("\n ");
	  }
	}
	if (ds->O)
	    printf(" O = %d, E = %d, K = %d, L = %d, M = %d",
		    ds->O, ds->E, K, L, M);
	printf("\"\n}\n");
}

/* print_align_header  ------------  print the top part of an alignment file */
void print_align_header (SEQ * seq1, SEQ * seq2, ss_t ss, gap_scores_t * ds,
   int K, int L)
{
	printf("#:lav\n");
	printf("s {\n");
	printf("  \"%s%s\" %d %d %d %d\n",
		SEQ_NAME(seq1), revflag(seq1), SEQ_FROM(seq1), SEQ_TO(seq1),
		(seq1->flags & SEQ_IS_REVCOMP) > 0, seq1->count);
	printf("  \"%s%s\" %d %d %d %d\n",
		SEQ_NAME(seq2), revflag(seq2), SEQ_FROM(seq2), SEQ_TO(seq2),
		(seq2->flags & SEQ_IS_REVCOMP) > 0, seq2->count);
	printf("}\n");
	
	printf("h {\n   \"%s%s\"\n   \"%s%s\"\n}\n", 
		SEQ_HEAD(seq1), revlabel(seq1), SEQ_HEAD(seq2), revlabel(seq2));
	// pacify lint
	if (ss == 0 || ds == 0 || K == 0 || L == 0) { L = 0; }
}

static void print_align_lav(int score, const uchar *seq1, const uchar *seq2,
 int beg1, int end1, int beg2, int end2, edit_script_t *S, seed_t* seed)
{
        int M, N, i, j, o, nG = 0, lG = 0;
	int total_run = 0, total_match = 0;
 
        M = end1 - beg1 + 1;
        N = end2 - beg2 + 1;
	
	printf("a {\n  s %d\n  b %d %d\n  e %d %d\n",
		score, beg1, beg2, end1, end2);
	while(seed){
	  if(seed->extra_type==IUPAC_EXTRA_INFO){
	    iupac_information *info=(iupac_information*) seed->extra_information;
	    printf("  I %d %d %d %d %d|%d %s (%s)\n",
		   seed->pos1, seed->pos2, 
		   seed->pos1+seed->len-1, seed->pos2+seed->len-1,
		   info->dist1,info->dist2, 
		   info->id, 
		   ori_names[info->ori]);
	  }
	  else if(seed->extra_type==WEIGHT_EXTRA_INFO){
	    weight_information *info=(weight_information*) seed->extra_information;
	    printf("  P %d %d %d %d %.*f|%.*f %s (%s)\n",
		   seed->pos1, seed->pos2, 
		   seed->pos1+seed->len-1, seed->pos2+seed->len-1,
		   info->gran_width,info->granularity*info->score1,
		   info->gran_width,info->granularity*info->score2,
		   info->id, ori_names[info->ori]);
	  }
	  seed=seed->next;
	}
        for (o = i = j = 0; i < M || j < N; ) {
                int start_i = i;
                int start_j = j; 
                int match = 0;
		int run;

	        run = es_rep_len(S, &o, seq1+beg1+i-1, seq2+beg2+j-1, &match);
		i += run; j += run;
		total_match += match;
		total_run += run;

		printf("  l %d %d %d %d %d\n",
			beg1+start_i, beg2+start_j, beg1+i-1, beg2+j-1, 
		      align_match_percent(run, match));

                if (i < M || j < N) {
		  nG++;
		  lG += es_indel_len(S, &o, &i, &j);
		}
		 
        }
	printf("}\n");
}

void prettyprint_align(align_t *a) {
  printf("p {\n  s %d\n  b %d %d\n  e %d %d\n",
	 a->score, a->beg1, a->beg2, a->end1, a->end2);
  edit_script_prettyprint(a->seq1+a->beg1-1, a->seq2+a->beg2-1, 
			  a->beg1, a->beg2, a->script, a->seed, 60);
  printf("}\n\n");
}

void print_align_list(align_t *a) {
  int output;


  ck_argc("print_align_list");
  if(!get_argval('A',&output)) output=1;
  while (a) {
    select_seeds_for_alignment(a);
    if(output!=23){
    if(output!=2)
      print_align_lav(a->score, a->seq1, a->seq2, a->beg1, a->end1,
		      a->beg2, a->end2, a->script, a->seed);
    if(output!=1)
      prettyprint_align(a);
    }
    free_seeds(a->seed);
    a->seed=NULL;
    a = a->next_align;
  }
}

void free_align_list(align_t *a)
{
	align_t *b;

	while (a) {
		edit_script_free(a->script);
		b = a->next_align;
		ckfree(a);
		a = b;
	}
}

static int percent_match(SEQ *sf1, SEQ *sf2, int len, int pos1, int pos2)
{
	uchar *p = SEQ_CHARS(sf1) + pos1;
	uchar *q = SEQ_CHARS(sf2) + pos2;
	int m = 0, i = len;

	while (i-- > 0)
		m += (*p++ == *q++);
	return align_match_percent (len, m);
}

void print_msp_list(SEQ *sf1, SEQ *sf2, struct msp_table *mt)
{
	int len, beg1, beg2, end1, end2, score, pm;
	msp_t *mp;

	/* print each maximal segment pair */
	for (mp=MSP_TAB_FIRST(mt); MSP_TAB_MORE(mt,mp); mp=MSP_TAB_NEXT(mp)) {
		beg1 = mp->pos1 + 1;
		beg2 = mp->pos2 + 1;
		len = mp->len - 1;
		end1 = beg1 + len;
		end2 = beg2 + len;
		score = mp->score;
		pm = percent_match(sf1, sf2, len, beg1 - 1, beg2 - 1);
		printf("a {\n  s %d\n  b %d %d\n  e %d %d\n",
		  score, beg1, beg2, end1, end2);
		printf("  l %d %d %d %d %d\n}\n", beg1, beg2, end1, end2, pm);
	}
}
