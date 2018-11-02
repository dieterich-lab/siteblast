#ifndef __lint
static const char rcsid[] =
"$Id: bz_site.c,v 4.0 2004/08/26 09:49:25 mmichael Exp mmichael $";
#endif

#ifdef DEBUG_VALGRIND
#include "valgrind.h"
#endif

#include "bz_all.h"
#include "astack.c"
#include "bz_iupac.decode"

typedef struct blast_table_s bt_t;

typedef struct {
  SEQ *seq1, *seq2;
  blast_table_t *bt;
  bool fst;
  enum blast_table_type btt;
} callback_params;

static blast_table_t *cons_blast_table_p(
				       pssm_t *pssm,
				       char orintations,
				       int highest_pattern_id, 
				       const signed char *encoding){
  blast_table_t *bt = ckallocz(sizeof(*bt));

  bt->p.encoding = encoding;
  bt->p.occur1   = ckcalloc(highest_pattern_id+1, sizeof(seed_t*));
  bt->p.occur2   = ckcalloc(highest_pattern_id+1, sizeof(seed_t*));
  bt->p.pssm     = pssm;
  bt->p.orintations=orintations;
  bt->p.highest_pattern_id=highest_pattern_id;
  return bt;
}
static blast_table_t *cons_blast_table_i(
				       Recognizer recognizer,
				       int highest_pattern_id, 
				       const signed char *encoding){
  blast_table_t *bt = ckallocz(sizeof(*bt));

  bt->i.encoding = encoding;
  bt->i.occur1   = ckcalloc(highest_pattern_id+1, sizeof(seed_t*));
  bt->i.occur2   = ckcalloc(highest_pattern_id+1, sizeof(seed_t*));
  bt->i.r        = recognizer;
  bt->i.highest_pattern_id=highest_pattern_id;
  return bt;
}
static blast_table_t *cons_blast_table_I(char **cons, char **ids, int type, char dist,
				       int number,
				       char orintations,
				       int highest_pattern_id, 
				       const signed char *encoding){
  blast_table_t *bt = ckallocz(sizeof(*bt));

  bt->I.encoding = encoding;
  bt->I.occur1   = ckcalloc(highest_pattern_id+1, sizeof(seed_t*));
  bt->I.occur2   = ckcalloc(highest_pattern_id+1, sizeof(seed_t*));
  bt->I.cons     = cons;
  bt->I.ids      = ids;
  bt->I.type     = type;
  bt->I.dist     = dist;
  bt->I.number   = number;
  bt->I.orintations=orintations;
  bt->I.highest_pattern_id=highest_pattern_id;
  return bt;
}


seed_t* consseed(void* extra_information, char extra_type, 
		 char len, int pos1, int pos2){
  seed_t *seed=ckalloc(sizeof(seed_t));
  seed->next=NULL;
  seed->extra_information=extra_information;
  seed->extra_type=extra_type;
  seed->len=len;
  seed->pos1=pos1;
  seed->pos2=pos2;
  return seed;
}

static void callback_build(void *closure, int pattern_length, char *instance, 
			   void *extra_information, char extra_type, 
			   int pattern_id, void *params){
  
  seed_t *seed;
  blast_table_t *bt =((callback_params*) params)->bt;
  bool fst=((callback_params*) params)->fst;
  enum blast_table_type btt=((callback_params*) params)->btt;
  seed_t **occur;


  switch(btt){
     case blast_table_type_p:
	if(fst)
	  occur=bt->p.occur1;
	else
	  occur=bt->p.occur2;
	break;
     case blast_table_type_i:
	if(fst)
	  occur=bt->i.occur1;
	else
	  occur=bt->i.occur2;
	break;
     case blast_table_type_I:
	if(fst)
	  occur=bt->I.occur1;
	else
	  occur=bt->I.occur2;
	break;
     default:
	fatalf("Internal Error: Don't know what to do with a blast_table of type %d.\n"
	       "Exiting.\t\t\t(%s:%d)\n",btt,__FILE__,__LINE__);
  }
  /*dprint
  printf("Found %s(%s) at pos %d (%*.*s) in seq. 1.\t\t(%s:%d)\n",
	 ((weight_information*) extra_information)->id,
	 ori_names[((weight_information*) extra_information)->ori],
	 instance-(char*)closure+1,pattern_length,pattern_length,instance,
	 __FILE__,__LINE__);*/

  if(!occur[pattern_id]) 
    occur[pattern_id] = 
      consseed(extra_information, extra_type, 
	       pattern_length,instance-(char*)closure+1,-1);
  else{
    seed = occur[pattern_id];
    occur[pattern_id]=consseed(extra_information, extra_type, 
			       pattern_length,instance-(char*)closure+1,-1);
    occur[pattern_id]->next=seed;
    /* Later, we must reverse all occurences*/
  }
} 

bool matched_by_pssm(weight_information *extra, singlepssm_t *pssm, char *input, char rev, char cmpl){
  int pos=0;
  int curScore;
  char alphabetToPos[256];
  int i;
  
  for(i=0;i<256;i++) alphabetToPos[i]=-1;
  if(cmpl)
    for(i=0;i<4;i++) alphabetToPos[dna_cmpl(alphabet[i])]=i;
  else 
    for(i=0;i<4;i++) alphabetToPos[alphabet[i]]=i;

  if(alphabetToPos[*input]==-1)return 0;
  if(rev){
    for(curScore=pssm->rows[pssm->length-1][alphabetToPos[*input++]];
	curScore+pssm->rev_max_remain_score[pos]>=pssm->t/pssm->granularity;
	curScore+=pssm->rows[pssm->length-1-(++pos)][alphabetToPos[*input++]]){
      if(pos==pssm->length-1){
	if(curScore<pssm->t/pssm->granularity){
	  return 0;
	}
	else{
	  extra->score1=curScore;
	  return 1;
	}
      }
      if(alphabetToPos[*input]==-1)return 0;
    }
    return 0;
  }
  else{ /*not rev*/
    for(curScore=pssm->rows[0][alphabetToPos[*input++]];
	pos<pssm->length && 
	  curScore+pssm->max_remain_score[pos]>=pssm->t/pssm->granularity;
	curScore+=pssm->rows[++pos][alphabetToPos[*input++]]){
      if(pos==pssm->length-1){
	if(curScore<pssm->t/pssm->granularity){
	  return 0;
	}
	else{
	  extra->score1=curScore;
	  return 1;
	}
      }
      if(alphabetToPos[*input]==-1)return 0;
    }
    return 0;
  }
  fatalf("Internal error: Should not reach this point %s:%l\n",__FILE__,__LINE__);
}

search_by_pssm(pssm_t *pssm, bool fst, SEQ *seq, CallbackExt f, void* closure, callback_params *params){
  char *input = SEQ_CHARS(seq);
  char *current = input;
  int length=SEQ_LEN(seq);
  char rev, cmpl;
  weight_information *extra=ckalloc(sizeof(weight_information));
  int num;
  pssm_t *curPssm;
  singlepssm_t *sPssm=ckalloc(sizeof(singlepssm_t));
  char orintations= params->bt->p.orintations;

  if(params->fst!=fst)
    fatalf("Internal Error: contradiction if to search 1st or 2nd sequenze\nExiting.\t\t(%s:%d)\n",
	   __FILE__, __LINE__);
  /*printpssm(pssm); dprint*/
  /*dprint
  printf("alternative_search: SEQ_LEN(seq)=%d\n"
	 "SEQ_CHARS(seq)=%s\n"
	 "closure       =%s\n"
	 "current-closure=%d\n\n",
	 length,input,(char*) closure,current-(char*) closure);
  */

  for(;*current;current++){
    for(curPssm=pssm,num=0;curPssm;curPssm=curPssm->next,num+=4){
      if(input+length-current>=curPssm->length){
	sPssm->length=curPssm->length;
	sPssm->granularity=curPssm->granularity;
	if(fst){
	  sPssm->max_remain_score=curPssm->max_remain_score1;
	  sPssm->rev_max_remain_score=curPssm->rev_max_remain_score1;
	  sPssm->t=curPssm->t1;
	  sPssm->rows=curPssm->rows1;
	  sPssm->order=curPssm->order1;
	}
	else{
	  sPssm->max_remain_score=curPssm->max_remain_score2;
	  sPssm->rev_max_remain_score=curPssm->rev_max_remain_score2;
	  sPssm->t=curPssm->t2;
	  sPssm->rows=curPssm->rows2;
	  sPssm->order=curPssm->order2;
	}

	for(rev=0;rev<2;rev++){
	  for(cmpl=0;cmpl<2;cmpl++){
	    extra->ori=rev+(2*cmpl);
	    if(orintations&(1<<extra->ori)){
	      if(matched_by_pssm(extra, sPssm, current, rev, cmpl)){
		extra->id=curPssm->id;
		extra->score2=-1;
		extra->granularity=curPssm->granularity;
		extra->gran_width=curPssm->gran_width;
		f(closure,curPssm->length,current,extra,WEIGHT_EXTRA_INFO, num+extra->ori, params);
		extra=ckalloc(sizeof(weight_information));
	      }
	    }
	  }
	}
      }
    }
  }
  free(extra);
  free(sPssm);
}

static inline bool matched_by_consensus(iupac_information *extra, char *con, int len, bool protein, char dist, char *input, char rev, char cmpl){
  int i;
  int d=0;


#define MATCH(decode)                        \
  for(i=0;i<len;i++){                        \
    d+=decode[con[rev?len-i-1:i]][input[i]]; \
    if(d>dist) return 0;                     \
  }

  if(protein){
    if(cmpl){
      fatalf("Internal Error. Cannot build complement for protein seq.\n Exiting\t\t(%s:%d)\n",
	    __FILE__,__LINE__);
    }
    else{
      MATCH(IUPAC_pro_decode);
    }
  }
  else{
    if(cmpl){
      MATCH(IUPAC_dna_compl_decode);
    }
    else{
      MATCH(IUPAC_dna_decode);
    }
  }
#undef MATCH

  extra->dist1=d;
  return 1;
}

search_by_consensus(char **cons, char **ids, bool protein, char dist, int number, SEQ *seq, 
			 CallbackExt f, void* closure, callback_params *params){
  char *input = SEQ_CHARS(seq);
  char *current = input;
  int length=SEQ_LEN(seq);
  char rev, cmpl;
  iupac_information *extra=ckalloc(sizeof(iupac_information));
  int num;
  char *con, *id;
  int conLen;
  int *lengths=ckalloc(number*sizeof(int));
  char orintations= params->bt->I.orintations;

  for(num=0;num<number;num++)
    lengths[num]=strlen(cons[num]);

  for(;*current;current++){
    for(con=cons[0],id=ids[0],conLen=lengths[0],num=0;
	num<number;
	num++){
      con=cons[num],id=ids[num],conLen=lengths[num];
      if(input+length-current>=conLen){
	for(rev=0;rev<2;rev++){
	  for(cmpl=0;cmpl<2;cmpl++){
	    extra->ori=rev+(2*cmpl);
	    if(orintations&(1<<extra->ori)){
	      if(matched_by_consensus(extra, con, conLen, protein, dist, current, rev, cmpl)){
		extra->id=id;
		extra->dist2=-1;
		f(closure, conLen, current, extra, IUPAC_EXTRA_INFO, 4*num+extra->ori, params);
		extra=ckalloc(sizeof(iupac_information));
	      }
	    }
	  }
	}
      }
    }
  }
  free(extra);
  free(lengths);
}

static blast_table_t *siteblast_table_p_enc_new(SEQ *seq, 
					      pssm_t *pssm,
					      int highest_pattern_id,
					      bz_flags_t *bz_flags,
					      const signed char *enc){
  blast_table_t *const bt = 
    cons_blast_table_p(pssm, bz_flags->o, highest_pattern_id, enc);

  void* closure=SEQ_CHARS(seq);
  callback_params *params = ckalloc(sizeof(callback_params));
  int i;

  params->bt=bt;
  params->fst=1;
  params->btt=blast_table_type_p;

#ifdef TRACK_SEARCH
  /*dprint*/
    printf("Starting search in 1st seq.\n");
#endif

    /*dprint
    printf("%s:%d\n",__FILE__, __LINE__);*/

    bt->p.pssm=pssm;
    search_by_pssm(pssm, 1/*fst*/, seq, 
		   & callback_build, closure, params);
    ckfree(params);params=NULL;

  /* Reverse all occurences*/
  for(i=0;i<=highest_pattern_id;i++){
    if(bt->p.occur1[i]){
      bt->p.occur1[i]=rev_seed(bt->p.occur1[i]);
    }
  }
#ifdef TRACK_SEARCH
#define MAX_OUTPUT_FOR_PRINT_BLAST_TABLE 5
  /*dprint*/
  printf("Search in 1st seq done.\n");

  /*dprint
  printf("The blast table:\n");
  for(i=0;i<=highest_pattern_id;i++){
    int j;
    seed_t *s;
	
    char *id;
    char ori;
    s=bt->s.occur[i];
    if(!s) printf("%10d %5d\n",i,0);
    else{
      if(s->extra_type==IUPAC_EXTRA_INFO){
	iupac_information *info=(iupac_information*) s->extra_information;
	id=info->id;
	ori=info->ori;
      }
      else{
	weight_information *info=(weight_information*) s->extra_information;
	id=info->id;
	ori=info->ori;
      }
      for(j=0;s;j++,s=s->next);
      printf("%10d %5d %s (%s)\n",i,j,id,ori_names[ori]);
      for(j=0,s=bt->s.occur[i];s&&j<MAX_OUTPUT_FOR_PRINT_BLAST_TABLE;j++,s=s->next)
	printf("\t%3d %5d;%5d (%2d) \"%s\"\n",j,s->pos1,s->pos2,s->len,
	       s->extra_type==WEIGHT_EXTRA_INFO?((weight_information*) s->extra_information)->id:"");
    }
    }*/
#endif
  

  return bt;
}

static blast_table_t *siteblast_table_i_enc_new(SEQ *seq, 
					      Recognizer recognizer,
					      int highest_pattern_id,
					      bz_flags_t *bz_flags,
					      const signed char *enc){
  blast_table_t *const bt = 
    cons_blast_table_i(recognizer, highest_pattern_id, enc);

  void* closure=SEQ_CHARS(seq);
  callback_params *params = ckalloc(sizeof(callback_params));
  int i;

  params->bt=bt;
  params->fst=1;
  params->btt=blast_table_type_i;

#ifdef TRACK_SEARCH
  /*dprint*/
  printf("Starting search in 1st seq.\n");
#endif

  search_for_ident_ext(recognizer, (char*) (SEQ_CHARS(seq)), 
			 & callback_build, closure, params);
  ckfree(params);params=NULL;

  /* Reverse all occurences*/
  for(i=0;i<=highest_pattern_id;i++){
    if(bt->i.occur1[i]){
      bt->i.occur1[i]=rev_seed(bt->i.occur1[i]);
    }
  }
#ifdef TRACK_SEARCH
#define MAX_OUTPUT_FOR_PRINT_BLAST_TABLE 5
  /*dprint*/
    printf("Search in 1st seq done.\n");

  /*dprint
  printf("The blast table:\n");
  for(i=0;i<=highest_pattern_id;i++){
    int j;
    seed_t *s;
	
    char *id;
    char ori;
    s=bt->i.occur[i];
    if(!s) printf("%10d %5d\n",i,0);
    else{
      if(s->extra_type==IUPAC_EXTRA_INFO){
	iupac_information *info=(iupac_information*) s->extra_information;
	id=info->id;
	ori=info->ori;
      }
      else{
	weight_information *info=(weight_information*) s->extra_information;
	id=info->id;
	ori=info->ori;
      }
      for(j=0;s;j++,s=s->next);
      printf("%10d %5d %s (%s)\n",i,j,id,ori_names[ori]);
      for(j=0,s=bt->s.occur[i];s&&j<MAX_OUTPUT_FOR_PRINT_BLAST_TABLE;j++,s=s->next)
	printf("\t%3d %5d;%5d (%2d) \"%s\"\n",j,s->pos1,s->pos2,s->len,
	       s->extra_type==WEIGHT_EXTRA_INFO?((weight_information*) s->extra_information)->id:"");
    }
    }*/
#endif

  return bt;
}

static blast_table_t *siteblast_table_I_enc_new(SEQ *seq, 
						char **cons, char **ids,
						int number,
						int highest_pattern_id,
						bz_flags_t *bz_flags,
						const signed char *enc){
  blast_table_t *const bt = 
    cons_blast_table_I(cons, ids, bz_flags->t, bz_flags->D, number, bz_flags->o, highest_pattern_id, enc);

  void* closure=SEQ_CHARS(seq);
  callback_params *params = ckalloc(sizeof(callback_params));
  int i;

  params->bt=bt;
  params->fst=1;
  params->btt=blast_table_type_I;

#ifdef TRACK_SEARCH
  /*dprint*/
  printf("Starting search in 1st seq.\n");
#endif

  search_by_consensus(cons, ids, bz_flags->t, bz_flags->D, number, seq, 
			& callback_build, closure, params);
  ckfree(params);params=NULL;

  /* Reverse all occurences*/
  for(i=0;i<=highest_pattern_id;i++){
    if(bt->I.occur1[i]){
      bt->I.occur1[i]=rev_seed(bt->I.occur1[i]);
    }
  }
#ifdef TRACK_SEARCH
#define MAX_OUTPUT_FOR_PRINT_BLAST_TABLE 5
  /*dprint*/
    printf("Search in 1st seq done.\n");

  /*dprint
  printf("The blast table:\n");
  for(i=0;i<=highest_pattern_id;i++){
    int j;
    seed_t *s;
	
    char *id;
    char ori;
    s=bt->s.occur[i];
    if(!s) printf("%10d %5d\n",i,0);
    else{
      if(s->extra_type==IUPAC_EXTRA_INFO){
	iupac_information *info=(iupac_information*) s->extra_information;
	id=info->id;
	ori=info->ori;
      }
      else{
	weight_information *info=(weight_information*) s->extra_information;
	id=info->id;
	ori=info->ori;
      }
      for(j=0;s;j++,s=s->next);
      printf("%10d %5d %s (%s)\n",i,j,id,ori_names[ori]);
      for(j=0,s=bt->s.occur[i];s&&j<MAX_OUTPUT_FOR_PRINT_BLAST_TABLE;j++,s=s->next)
	printf("\t%3d %5d;%5d (%2d) \"%s\"\n",j,s->pos1,s->pos2,s->len,
	       s->extra_type==WEIGHT_EXTRA_INFO?((weight_information*) s->extra_information)->id:"");
    }
    }*/
#endif
  

  return bt;
}

blast_table_t *siteblast_table_new(SEQ *seq, int max_seq_length, bz_flags_t *bz_flags){
  static int highest_pattern_id;
  static char *id_list=NULL;
  blast_table_t *bt;
  static Recognizer recognizer=NULL; /*For param I*/
  static char **cons=NULL;           /*For param i*/
  static char **ids=NULL;            /*For param i*/
  static int number=0;               /*For param i*/

  /*dprint
  printf("%s:%d\n",__FILE__, __LINE__);*/

  if(bz_flags->I){
    if(bz_flags->i) fatalf("Use only I=<file> or i=<file>. Exiting.\t\t(%s:%d)\n",
			   __FILE__,__LINE__);
    if(bz_flags->P) fatalf("Use only I=<file> or P=<file>. Exiting.\t\t(%s:%d)\n",
			   __FILE__,__LINE__);
    if(!(recognizer&&id_list)){
      if(recognizer||id_list)
	fatalf("Internal Error. Exiting.\t\t(%s:%d)\n",__FILE__,__LINE__);
      else{
	recognizer = new_iupacrecogniser(bz_flags, &highest_pattern_id, &id_list);
      }
    }
    bt= siteblast_table_i_enc_new(seq, recognizer, highest_pattern_id, bz_flags, fasta_encoding);
    bt->i.id_list=id_list;
    return bt;
  }
  if(bz_flags->i){
    if(bz_flags->P) fatalf("Use only i=<file> or P=<file>. Exiting.\t\t(%s:%d)\n",
			   __FILE__,__LINE__);
    if(!(cons&&ids&&id_list&&number)){
      if(cons||ids||id_list||number)
	fatalf("Internal Error. Exiting.\t\t(%s:%d)\n",__FILE__,__LINE__);
      else{
	getCons(&cons, &ids, &number,bz_flags, &highest_pattern_id, &id_list);
      }
    }
    bt= siteblast_table_I_enc_new(seq, cons, ids, number, highest_pattern_id, bz_flags, fasta_encoding);
    bt->I.id_list=id_list;
    return bt;
  }
  if(bz_flags->P){
    /*The PSSMs are depending on GC content must be recomputed for each seq. pair*/
    pssm_t *pssm=get_PSSM(bz_flags, &id_list, max_seq_length, &highest_pattern_id);
    bt= siteblast_table_p_enc_new(seq, pssm, highest_pattern_id, bz_flags, fasta_encoding);
    bt->p.id_list=id_list;
    return bt;
  }
  else
    fatalf("siteblast_table_new called with I==i==P==NULL. Exiting.\t\t(%s:%d)\n",
	   __FILE__,__LINE__);
}

/* -----------------------   search the second sequence   --------------------*/

static double entropy(const uchar *s, int n)
{
	double pA, pC, pG, pT, qA, qC, qG, qT, e;
	int count[128];
	int i, cA, cC, cG, cT;

	for (i = 0; i < 128; ++i)
		count[i] = 0;
	for (i = 0; i < n; ++i)
		++count[s[i]];
	cA = count['A'];
	cC = count['C'];
	cG = count['G'];
	cT = count['T'];
	pA = ((double)count['A']) / ((double)n);
	pC = ((double)count['C']) / ((double)n);
	pG = ((double)count['G']) / ((double)n);
	pT = ((double)count['T']) / ((double)n);
	qA = (cA ? log(pA) : 0.0);
	qC = (cC ? log(pC) : 0.0);
	qG = (cG ? log(pG) : 0.0);
	qT = (cT ? log(pT) : 0.0);
	e = -(pA*qA + pC*qC + pG*qG + pT*qT)/log(4.0);

	return e;
}

msp_table_t *process_siteblast_hits(SEQ *seq1, SEQ *seq2, blast_table_t *bt, 
				    enum blast_table_type btt,
				    msp_table_t *msp_tab,
				    ss_t ss, int X, int K, int P, int T){
  const int len1 = SEQ_LEN(seq1);
  const int len2 = SEQ_LEN(seq2);
  static comb_t *diag_lev = NULL;
  seed_t **occur1, **occur2;
  int highest_pattern_id, i;

  if (seq1 == 0) return 0;
  if (seq2 == 0) return 0;
  if (bt == 0) return 0;
  if (T != 0) fprintf(stderr,"Using T!=0 and P, I or i is not supported. Ignoring T!\n");

#ifdef TRACK_SEARCH
  printf("Computing MSPs started\n");
#endif

  { 
    unsigned int n = (len1+len2+1)*sizeof(int);
    n = roundup(n, 64*1024);
    diag_lev = comb_resize(diag_lev,n,len1);
  }

  switch(btt){
     case blast_table_type_p:
	occur1=bt->p.occur1;
	occur2=bt->p.occur2;
	highest_pattern_id=bt->p.highest_pattern_id;
	break;
     case blast_table_type_i:
	occur1=bt->i.occur1;
	occur2=bt->i.occur2;
	highest_pattern_id=bt->i.highest_pattern_id;
	break;
     case blast_table_type_I:
	occur1=bt->I.occur1;
	occur2=bt->I.occur2;
	highest_pattern_id=bt->I.highest_pattern_id;
	break;
     default:
	fatalf("Internal Error: Don't know what to do with a blast_table of type %d.\n"
	       "Exiting.\t\t\t(%s:%d)\n",btt,__FILE__,__LINE__);
  }


  /*dprint
  for(i=0;i<=highest_pattern_id;i++){
    if(occur1[i]&&occur2[i]){
      printf("occur1[%d]:\t\t(%s:%d)\n",i,__FILE__,__LINE__);
      print_seeds(occur1[i]);
    }
  }
  */

  for(i=0;i<=highest_pattern_id;i++){
    seed_t *occured1, *occured2, *seed;
    for(occured2=occur2[i];occured2;occured2=occured2->next)
      for(occured1=occur1[i];occured1;occured1=occured1->next){
	seed = cloneseed(occured1,0);
	seed->next=NULL;
	/*dprint
	printf("%s:%d\n",__FILE__,__LINE__);
	print_seeds(seed);
	*/
	seed->pos2=occured2->pos1;
	/*dprint
	printf("%s:%d\n",__FILE__,__LINE__);
	print_seeds(seed);
	*/
	if(seed->extra_type==WEIGHT_EXTRA_INFO){
	  weight_information *infoFrom=(weight_information*) occured2->extra_information;
	  weight_information *infoTo  =(weight_information*) seed->extra_information;
	  if(infoTo->score2!=-1){
	    fatalf("Information for 2nd Seq. must not be set. (%s:%d)\ni=%d\n",
		   __FILE__,__LINE__,i);
	  }
	  infoTo->score2=infoFrom->score1;
	}
	else if(seed->extra_type==IUPAC_EXTRA_INFO){
	  iupac_information *infoFrom=(iupac_information*) occured2->extra_information;
	  iupac_information *infoTo  =(iupac_information*) seed->extra_information;
	  if(infoTo->dist2!=-1){
	    fatalf("Information for 2nd Seq. must not be set. (%s:%d)\ni=%d\n",
		   __FILE__,__LINE__,i);
	  }
	  infoTo->dist2=infoFrom->dist1;
	}
	else
	  fatalf("Internal Error: %s:%d\n",__FILE__,__LINE__);

	insert_seed_into_list(seed->pos1-seed->pos2,seed);

	msp_extend_hit(msp_tab, seq1, seq2, 
			    ss, X, K, seed->len, P, 
			    seed->pos1-1+seed->len, seed->pos2-1+seed->len,
			    diag_lev);
      }
  }

  comb_clear(diag_lev);

#ifdef TRACK_SEARCH
  printf("Computing MSPs done.\n");
#endif

  return msp_tab;
}
void correct_bt_for_2nd_search(int length , blast_table_t *bt, enum blast_table_type btt){
  seed_t **occur2, *tmp;
  int highest_pattern_id, i;
  char *ori;

  /*dprint
  printf("correct_bt_for_2nd_search(%d, %p) called.\n", length, bt);*/
  
  switch(btt){
     case blast_table_type_p:
	occur2=bt->p.occur2;
	highest_pattern_id=bt->p.highest_pattern_id;
	break;
     case blast_table_type_i:
	occur2=bt->i.occur2;
	highest_pattern_id=bt->i.highest_pattern_id;
	break;
     case blast_table_type_I:
	occur2=bt->I.occur2;
	highest_pattern_id=bt->I.highest_pattern_id;
	break;
     default:
	fatalf("Internal Error: Don't know what to do with a blast_table of type %d.\n"
	       "Exiting.\t\t\t(%s:%d)\n",btt,__FILE__,__LINE__);
  }
  if(length){
    /*dprint
    printf("Have to changed stored seeds:\n");
    printf("1st: index in occur2\n");*/
    for(i=0;i<=highest_pattern_id;i+=4){
      tmp=occur2[i]; 
      occur2[i]=occur2[i+3]; /*sense <-> rev. cmpl.*/
      occur2[i+3]=tmp;
      tmp=occur2[i+1]; 
      occur2[i+1]=occur2[i+2]; /*rev. <-> cmpl.*/
      occur2[i+2]=tmp;
    }
    /*dprint
      printf("2nd: stored pos2\n");*/
    for(i=0;i<=highest_pattern_id;i++)
      for(tmp=occur2[i];tmp;tmp=tmp->next){
	if(tmp->extra_type==WEIGHT_EXTRA_INFO){
	  weight_information *info =(weight_information*) tmp->extra_information;
	  ori=&(info->ori);
	}
	else if(tmp->extra_type==IUPAC_EXTRA_INFO){
	  iupac_information *info =(iupac_information*) tmp->extra_information;
	  ori=&(info->ori);
	}
	else
	  fatalf("Internal Error: %s:%d\n",__FILE__,__LINE__);
	
	*ori=3-*ori;
	tmp->pos1=length-tmp->pos1-tmp->len+2;
      }
  }
  else
  {
    /*dprint
    printf("Have to delete stored seeds:\n");*/
    for(i=0;i<=highest_pattern_id;i++){
      free_seeds(occur2[i]);
      occur2[i]=NULL;
    }
  }
  /*dprint
  printf("%s:%d\n",__FILE__, __LINE__);*/
}

msp_table_t * siteblast_search(SEQ *seq1, SEQ *seq2, blast_table_t *bt, enum blast_table_type btt,
			      msp_table_t *msp_tab,
			      ss_t ss, int X, int K, int P, int T){
  const int len1 = SEQ_LEN(seq1);
  const int len2 = SEQ_LEN(seq2);
  callback_params *params;
  int i;

  if (seq1 == 0) return 0;
  if (seq2 == 0) return 0;
  if (bt == 0) return 0;
  if (T != 0) fprintf(stderr,"Using T!=0 and P, I or i is not supported. Ignoring T!\n");

  params = ckalloc(sizeof(callback_params));

  params->seq1=seq1;
  params->seq2=seq2;
  params->bt=bt;
  params->fst=0;
  params->btt=btt;

#ifdef TRACK_SEARCH
  /*dprint*/
  printf("Starting search\n");
#endif

  if(btt==blast_table_type_i){

    /*dprint
    printf("%s:%d\n",__FILE__, __LINE__);
    for(i=0;i<=bt->i.highest_pattern_id;i++){
      printf("bt->i.occur2[%d]:\n",i);
      print_seeds(bt->i.occur2[i]);
    }
    printf("%s:%d\n",__FILE__, __LINE__);
*/
    search_for_ident_ext(bt->i.r,              (char*) (SEQ_CHARS(seq2)), 
			 & callback_build, SEQ_CHARS(seq2),params);
    /* Reverse all occurences*/
    for(i=0;i<=bt->i.highest_pattern_id;i++){
      if(bt->i.occur2[i]){
	bt->i.occur2[i]=rev_seed(bt->i.occur2[i]);
      }
    }
  }
  else if(btt==blast_table_type_p){

    /*dprint
    printf("%s:%d\n",__FILE__, __LINE__);
    for(i=0;i<=bt->p.highest_pattern_id;i++){
      printf("bt->p.occur2[%d]:\n",i);
      print_seeds(bt->p.occur2[i]);
    }
    printf("%s:%d\n",__FILE__, __LINE__);
*/
    search_by_pssm ((pssm_t *)(bt->p.pssm), 0 /*!fst*/, seq2, 
		    & callback_build, SEQ_CHARS(seq2), params);
    /* Reverse all occurences*/
    for(i=0;i<=bt->p.highest_pattern_id;i++){
      if(bt->p.occur2[i]){
	bt->p.occur2[i]=rev_seed(bt->p.occur2[i]);
      }
    }
  }
  else if(btt==blast_table_type_I){

    /*dprint
    printf("%s:%d\n",__FILE__, __LINE__);
    for(i=0;i<=bt->I.highest_pattern_id;i++){
      printf("bt->I.occur2[%d]:\n",i);
      print_seeds(bt->I.occur2[i]);
    }
    printf("%s:%d\n",__FILE__, __LINE__);
*/
    search_by_consensus(bt->I.cons, bt->I.ids, bt->I.type, bt->I.dist, bt->I.number, seq2, 
			& callback_build, SEQ_CHARS(seq2),params);
    /* Reverse all occurences*/
    for(i=0;i<=bt->I.highest_pattern_id;i++){
      if(bt->I.occur2[i]){
	bt->I.occur2[i]=rev_seed(bt->I.occur2[i]);
      }
    }
  }
  else {
    fatalf("Internal Error: Don't know what to do with a blast_table of type %d.\n"
	   "Exiting.\t\t\t(%s:%d)\n",btt,__FILE__,__LINE__);
  }
  ckfree(params);params=NULL;

  

#ifdef TRACK_SEARCH
  /*dprint*/
  printf("Search done\n");
#endif

  return process_siteblast_hits(seq1, seq2, bt, btt, msp_tab, ss, X, K, P, T);
}



void siteblast_table_p_free(blast_table_t *bt){
  register int i;
  if (bt == 0) return;

  for(i=0;i<=bt->p.highest_pattern_id;i++){
    free_seeds(bt->p.occur1[i]);
    free_seeds(bt->p.occur2[i]);
  }
  ZFREE(bt->p.occur1);
  ZFREE(bt->p.occur2);
  free_pssms(bt->p.pssm);
  ZFREE(bt);
}
void siteblast_table_i_free(blast_table_t *bt){
  register int i;
  if (bt == 0) return;

  for(i=0;i<=bt->i.highest_pattern_id;i++){
    free_seeds(bt->i.occur1[i]);
    free_seeds(bt->i.occur2[i]);
  }
  ZFREE(bt->i.occur1);
  ZFREE(bt->i.occur2);
  ZFREE(bt);
}
void siteblast_table_I_free(blast_table_t *bt){
  register int i;
  if (bt == 0) return;

  for(i=0;i<=bt->I.highest_pattern_id;i++){
    free_seeds(bt->I.occur1[i]);
    free_seeds(bt->I.occur2[i]);
  }
  ZFREE(bt->I.occur1);
  ZFREE(bt->I.occur2);
  ZFREE(bt);
}
