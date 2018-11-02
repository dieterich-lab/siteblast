#include "util.h"
#include "seq.h"
#include "bz_all.h"
#include <assert.h>

static const char rcsid[]=
"$Id: edit.c,v 4.0 2004/08/26 09:49:25 mmichael Exp mmichael $";

static edit_op_t *ed_ops_realloc(edit_op_t *op, unsigned int n);
static edit_op_t edit_op_cons(unsigned int op, unsigned int val);
static edit_op_t edit_op_inc(edit_op_t op, unsigned int n);
static edit_op_t edit_op_inc_last(edit_script_t *es, unsigned int n);
static edit_script_t *edit_script_concat(edit_script_t *es, edit_script_t *et);
static edit_script_t *edit_script_copy(edit_script_t *);
static edit_script_t *edit_script_fin(edit_script_t *data);
static edit_script_t *edit_script_init(edit_script_t *es);
static int edit_script_more(edit_script_t *data, unsigned int op, unsigned int k);
static int edit_script_put(edit_script_t *es, unsigned int op, unsigned int n);
static int edit_script_ready(edit_script_t *es, unsigned int n);
static int edit_script_readyplus(edit_script_t *es, unsigned int n);

unsigned int edit_opc_get(edit_op_t op)
{
	return op & EDIT_OP_MASK;
}

unsigned int edit_val_get(edit_op_t op)
{
	return op >> 2;
}

static edit_op_t edit_op_cons(unsigned int op, unsigned int val)
{
	return (val << 2) | (op & EDIT_OP_MASK);
}

static edit_op_t edit_op_inc(edit_op_t op, unsigned int n)
{
	return edit_op_cons(edit_opc_get(op), edit_val_get(op) + n);
}

static edit_op_t edit_op_inc_last(edit_script_t *es, unsigned int n)
{
	edit_op_t *last;
	assert (es->num > 0);
	last = &(es->op[es->num-1]);
	*last = edit_op_inc(*last, n);
	return *last;
}

void edit_script_prnt(edit_script_t *es)
{
	unsigned int i;

	printf("num=%d size=%d\n", es->num, es->size);
	for (i = 0; i < es->num; ++i) {
		edit_op_t this = es->op[i];
		assert(this != 0);
		printf("%d:%c:%d ", i, "?IDR"[edit_opc_get(this)],
		  edit_val_get(this));
		if (i%8 == 7) putchar('\n');
	}
	putchar('\n');
}

static edit_script_t *edit_script_init(edit_script_t *es)
{
	es->op = 0;
	es->size = es->num = 0;
	es->last = 0;
	edit_script_ready(es, 8);
	return es;
}

static edit_op_t *ed_ops_realloc(edit_op_t *op, unsigned int n)
{
	return ckrealloc(op, n*sizeof(edit_op_t)); 
	/* XXX - assumes posix realloc */
}

static int edit_script_ready(edit_script_t *es, unsigned int n)
{
	edit_op_t *p;
	unsigned int m = n + n/2;

	if (es->size <= n) {
		p = ed_ops_realloc(es->op, m);
		if (p == 0) {
			return 0;
		} else {
			es->op = p;
			es->size = m;
		}
	}
	return 1;
}

static int edit_script_readyplus(edit_script_t *es, unsigned int n)
{
	if (es->size - es->num <= n)
		return edit_script_ready(es, n + es->num);
	return 1;
}

static int edit_script_put(edit_script_t *es, unsigned int op, unsigned int n)
{
	if (!edit_script_readyplus(es, 2))
		return 0;
	es->last = op;
	assert(op != 0);
	es->op[es->num] = edit_op_cons(op, n);
	es->num += 1;
	es->op[es->num] = 0; /* sentinal */
	return 1;
}

static edit_script_t *edit_script_fin(edit_script_t *es)
{
	edit_op_t *p = ed_ops_realloc(es->op, es->num);
	if (!p)
		return 0;
	es->op = p;
	es->size = es->num;
	return es;
}

static int edit_script_more(edit_script_t *data, unsigned int op, unsigned int k)
{
	if (op == EDIT_OP_ERR)
		fatalf("edit_script_more: bad opcode %d:%d", op, k);
	if (edit_opc_get(data->last) == op)
		edit_op_inc_last(data, k);
	else
		edit_script_put(data, op, k);
	return 0;
}

int edit_script_del(edit_script_t *data, unsigned int k)
{
	return edit_script_more(data, EDIT_OP_DEL, k);
}

int edit_script_ins(edit_script_t *data, unsigned int k)
{
	return edit_script_more(data, EDIT_OP_INS, k);
}

int edit_script_rep(edit_script_t *data, unsigned int k)
{
	return edit_script_more(data, EDIT_OP_REP, k);
}

edit_script_t *edit_script_reverse_inplace(edit_script_t *es)
{
	unsigned int i;
	const unsigned int num = es->num;
	const unsigned int mid = num/2;
	const unsigned int end = num-1;

	for (i = 0; i < mid; ++i) {
		const edit_op_t t = es->op[i];
		es->op[i] = es->op[end-i];
		es->op[end-i] = t;
	}
	return es;
}

edit_script_t *edit_script_new(void)
{
	edit_script_t *es = ckallocz(sizeof(*es));
	if (!es)
		return 0;
	return edit_script_init(es);
}

edit_script_t *edit_script_free(edit_script_t *es)
{
	if (es) {
		if (es->op)
			ckfree(es->op);
		memset(es, 0, sizeof(*es));
		ckfree(es);
	}
	return 0;
}

/* deep copy of es */
static edit_script_t *edit_script_copy(edit_script_t *es)
{
	edit_script_t *nes = ckallocz(sizeof(*nes));
	unsigned int size = sizeof(*nes->op) * es->num;
	nes->op = ckallocz(size);
	if (!nes->op)
		return 0; /* XXX - leak */
	memcpy(nes->op, es->op, size);
	nes->size = nes->num = es->num;
	nes->last = 0;
	return nes;
}

edit_op_t *edit_script_first(edit_script_t *es)
{
	return es->num > 0 ? &es->op[0] : 0;
}

edit_op_t *edit_script_next(edit_script_t *es, edit_op_t *op)
{
	/* XXX - assumes flat address space */
	if (&es->op[0] <= op && op < &es->op[es->num-1])
		return op+1;
	else
		return 0;
}

/* build a new script from es and et */
static edit_script_t *edit_script_concat(edit_script_t *es, edit_script_t *et)
{
	edit_op_t *op;
	edit_script_t *eu;

	if (!(eu = edit_script_new()))
		return 0;
	for (op = edit_script_first(es); op; op = edit_script_next(es, op))
		edit_script_more(eu, edit_opc_get(*op), edit_val_get(*op));
	for (op = edit_script_first(et); op; op = edit_script_next(et, op))
		edit_script_more(eu, edit_opc_get(*op), edit_val_get(*op));
	return edit_script_fin(eu);
}

/* add et to es */
edit_script_t *edit_script_append(edit_script_t *es, edit_script_t *et)
{
	edit_op_t *op;

	for (op = edit_script_first(et); op; op = edit_script_next(et, op))
		edit_script_more(es, edit_opc_get(*op), edit_val_get(*op));
	return es;
}

int es_rep_len(edit_script_t *S, int *n, const uchar *p, const uchar *q, int *match)
{
	int len;

	len = *match = 0;
	while (((unsigned int)*n < S->num) && (edit_opc_get(S->op[*n]) == EDIT_OP_REP)) {
		int num = edit_val_get(S->op[*n]);
		len += num;
		while (num-- > 0)
		    /* masked regions are lower case */
                    if (toupper(*p++) == toupper(*q++)) ++(*match);
		*n += 1;
	}
	return len;
}

int es_indel_len(edit_script_t *S, int *n, int *i, int *j)
{
	if (S->num <= (unsigned int)*n)
		return 0;
	else {
		edit_op_t op = S->op[*n];
		int len = edit_val_get(op);

		switch (edit_opc_get(op)) {
		case EDIT_OP_INS:
			*j += len;
			break;
		case EDIT_OP_DEL:
			*i += len;
			break;
		default:
			fatalf("es_indel_len: cannot happen!");
		}
		*n += 1;
		return len;
	}
}

int es_len(edit_script_t *S){
  int i,len;
  for(i=len=0;i<S->num;i++) {
    len+=edit_val_get(S->op[i]);
  }
  return len;
}

#ifdef TESTING
int main()
{
	int i;
	edit_script_t *es, *et, *eu;

	es = edit_script_new();
	for (i = 0; i < 13; ++i) {
		unsigned int r = random() % 4;
		unsigned int s = random() % 100;
		if (1 <= r && r <= 3) edit_script_more(es, r, s);
	}
	edit_script_prnt(es);
	et = edit_script_reverse_inplace(edit_script_copy(es));
	edit_script_prnt(et);
	eu = edit_script_concat(es, et);
	edit_script_prnt(eu);
	eu = edit_script_append(eu, es);
	edit_script_prnt(eu);

	es = edit_script_free(es);
	et = edit_script_free(et);
	eu = edit_script_free(eu);
	exit(0);
}
#endif


bool seed_in_alignment(seed_t *seed, align_t *a){
  edit_script_t *es;
  int current_script,len,op,i;
  int pos1,pos2;
  int offset=seed->pos1-seed->pos2;

  if(a->beg1>seed->pos1||a->end1<seed->pos1+seed->len-1||
     a->beg2>seed->pos2||a->end2<seed->pos2+seed->len-1)
    return 0;
  es=a->script;
  pos1=a->beg1-1;
  pos2=a->beg2-1;

  /*dprint
#define DPRINTF\
  printf("seed->pos1=%4d,    offset=%d\n      pos1=%4d, pos1-pos2=%d\t\t(%s:%d)\n",\
	 seed->pos1, offset, pos1, pos1-pos2, __FILE__, __LINE__);
  DPRINTF;
  */
  for(current_script=0;current_script<es->num;current_script++){
    len=edit_val_get(es->op[current_script]);
    op =edit_opc_get(es->op[current_script]);
    if(op==EDIT_OP_REP){
      if(pos1-pos2==offset){
	for(i=0;i<len;i++){
	  pos1++;
	  pos2++;
	  if((pos1>=seed->pos1)&&(pos1<seed->pos1+seed->len)){ 
	    return 1;
	  }
	}
      }
      else{
	pos1+=len;
	pos2+=len;
      }
    }
    else if(op==EDIT_OP_DEL){
      pos1+=len;
    }
    else if(op==EDIT_OP_INS){
      pos2+=len;
    }
  }
  return 0;
}

void select_seeds_for_alignment(align_t *a){
  edit_script_t *es=a->script;
  int *offsets;
  int noffsets=0;
  int pos1=a->beg1, pos2=a->beg2;
  int current_script,len,op,i;
  bool flag;


  seed_t *seed=a->seed, *cor=NULL, *corEnd=NULL;
  char fst_found=0;

  offsets=ckcalloc(es->num,sizeof(int));
  offsets[noffsets++]=pos1-pos2;
  for(current_script=0;current_script<es->num;current_script++){
    len=edit_val_get(es->op[current_script]);
    op =edit_opc_get(es->op[current_script]);
    if(op==EDIT_OP_DEL){
      pos1+=len;
    }
    else if(op==EDIT_OP_INS){
      pos2+=len;
    }
    if(op!=EDIT_OP_REP){
      for(i=0,flag=1;flag&&i<noffsets;i++)
	if(offsets[i]==pos1-pos2) flag=0;
      if(flag) 
	offsets[noffsets++]=pos1-pos2;
    }
  }

  /*dprint
  edit_script_prnt(es);
  printf("offsets:\n");
  for(i=0;i<es->num;i++)
    printf("%6d ",i);
  printf("\n");
  for(i=0;i<es->num;i++)
    printf("%6d ",offsets[i]);
  printf("\n");
  */
  
  for(i=0;i<noffsets;i++){
    seed=get_seed_from_list(offsets[i]);
    while(seed){
      if(seed_in_alignment(seed,a)){
	if(fst_found){
	  corEnd->next=cloneseed(seed,0);
	  corEnd=corEnd->next;
	}
	else{
	  corEnd=cor=cloneseed(seed,0);
	  fst_found=1;
	}
	corEnd->next=NULL;
      }
      seed=seed->next;
    }
  }
  ckfree(offsets);offsets=NULL;
  /*dprint
    print_seeds(cor);*/
  a->seed=sort_seed(cor);
}


void edit_script_prettyprint(uchar* seq1, uchar* seq2, int start1, int start2,
			     edit_script_t* es,
			     seed_t *seed, int linelength){

  const size_t minlength=2*linelength;
  int currentlength=minlength, oldneededlength, neededlength=0;
  int current_script,j,k;
  uchar *pseq1, *pseq2, *newseq1, *newseq2;
  int  *pos1,  *pos2,  *newpos1, *newpos2;

  pseq1=ckcalloc(minlength+1,sizeof(char));
  pseq2=ckcalloc(minlength+1,sizeof(char));
  pos1=ckcalloc(minlength,sizeof(int));
  pos2=ckcalloc(minlength,sizeof(int));
  for(current_script=0;current_script<es->num;current_script++){
    int len=edit_val_get(es->op[current_script]);
    int op=edit_opc_get(es->op[current_script]);
    newseq1=ckcalloc(len+1,sizeof(char));
    newseq2=ckcalloc(len+1,sizeof(char));
    newpos1=ckcalloc(len,sizeof(int));
    newpos2=ckcalloc(len,sizeof(int));

    if(op==EDIT_OP_REP){
      strncpy(newseq1,seq1,len);seq1+=len;
      strncpy(newseq2,seq2,len);seq2+=len;
      for(j=0;j<len;j++,start1++,start2++)
	newpos1[j]=start1,
	  newpos2[j]=start2;
    }
    else if(op==EDIT_OP_DEL){
      uchar *tmp=newseq2;
      strncpy(newseq1,seq1,len);seq1+=len;
      for(j=0;j<len;j++,start1++)
	*tmp++='-',
	  newpos1[j]=start1,
	  newpos2[j]=-1;
    }
    else if(op==EDIT_OP_INS){
      uchar *tmp=newseq1;
      strncpy(newseq2,seq2,len);seq2+=len;
      for(j=0;j<len;j++,start2++)
	*tmp++='-',
	  newpos1[j]=-1,
	  newpos2[j]=start2;
    }
    oldneededlength=neededlength;
    neededlength+=len;
    if(neededlength>currentlength){
      pos2 =ckrealloc(pos2,  neededlength   *sizeof(int));
      pos1 =ckrealloc(pos1,  neededlength   *sizeof(int));
      pseq2=ckrealloc(pseq2,(neededlength+1)*sizeof(char));
      pseq1=ckrealloc(pseq1,(neededlength+1)*sizeof(char));
    }
    memcpy(pos1 +oldneededlength,newpos1,len*sizeof(int)); free(newpos1);
    memcpy(pos2 +oldneededlength,newpos2,len*sizeof(int)); free(newpos2);
    memcpy(pseq1+oldneededlength,newseq1,len+1);           free(newseq1);
    memcpy(pseq2+oldneededlength,newseq2,len+1);           free(newseq2);
    /* Output? */
    while(neededlength>0&&(neededlength>=linelength||
			   current_script==(es->num-1))){
      int start1,end1,start2,end2,i;
      uchar *line =ckcalloc(linelength+20,sizeof(char));
      if(current_script==es->num-1&&linelength>neededlength) 
	linelength=neededlength;

      for(end1=pos1[j=linelength-1];end1==-1&&j>0;               end1=pos1[--j]);
      for(start1=pos1[j=0];         start1==-1&&j<linelength-1;start1=pos1[++j]);
      for(end2=pos2[j=linelength-1];end2==-1&&j>0;               end2=pos2[--j]);
      for(start2=pos2[j=0];         start2==-1&&j<linelength-1;start2=pos2[++j]);

      /*mark seed in 1st seq., but only if there is a gap...*/
      {
	seed_t *tmpseed=seed;
	uchar *b, *a, blen;
	uchar *ext_id;
	int i;
	bool mis=0;
	for(;tmpseed;tmpseed=tmpseed->next,mis=0){
	  if(((tmpseed->pos1+tmpseed->len>start1)&&(tmpseed->pos1<=end1))){
	    int markstart=-1;
	    for(i=0;i<linelength&&i<neededlength;i++){
	      if(markstart!=-1&&pos1[i]==-1){
		strcat(line," ");
	      }
	      if((tmpseed->pos1+tmpseed->len>pos1[i])&&(tmpseed->pos1<=pos1[i])){
		if(markstart==-1) markstart=i;
		if((tmpseed->pos2+tmpseed->len>pos2[i])&&(tmpseed->pos2<=pos2[i])){
		  strcat(line,"^");
		}
		else{
		  strcat(line,"~");
		  mis=1;
		}
	      }
	    }
	    if(mis){
	      if(tmpseed->extra_type==IUPAC_EXTRA_INFO){
		iupac_information *info=(iupac_information*) tmpseed->extra_information;
		ext_id=ckalloc(sizeof(char)*(1+strlen(info->id)+
					     strlen(ori_names[info->ori])+7+
					     ((info->dist1||info->dist2)?20:0)));
		strcpy(ext_id,info->id);
		sprintf(ext_id,"%s (%s)",ext_id,ori_names[info->ori]);
		if(info->dist1||info->dist2) 
		  sprintf(ext_id,"%s dist:%d|%d",ext_id,info->dist1,info->dist2);
	      }
	      else if(tmpseed->extra_type==WEIGHT_EXTRA_INFO){
		weight_information *info=(weight_information*) tmpseed->extra_information;
		ext_id=ckalloc(sizeof(char)*(1+strlen(info->id)+
					     strlen(ori_names[info->ori])+7+
					     10+2*info->gran_width));
		strcpy(ext_id,info->id);
		sprintf(ext_id,"%s (%s)",ext_id,ori_names[info->ori]);
		sprintf(ext_id,"%s score:%.*f|%.*f",ext_id,
			info->gran_width,info->granularity*info->score1,
			info->gran_width,info->granularity*info->score2);
	      }
	      else ext_id="Unknown";

	      if((blen=strlen(ext_id))<markstart) 
		b=ext_id,a="";
	      else 
		a=ext_id,b="",blen=0;
	      for(k=blen;k<markstart;k++) printf(" ");
	      printf("%s",b);
	      printf("%s",line);
	      printf("%s\n",a);
	      free(ext_id);
	    }
	    free(line);
	    line =ckcalloc(linelength+20,sizeof(char));
	  }
	}
      }

      /*1st pos line*/
      for(end1=pos1[j=linelength-1];end1==-1&&j>0;               end1=pos1[--j]);
      for(start1=pos1[j=0];         start1==-1&&j<linelength-1;start1=pos1[++j]);
      for(i=0;i<j;i++) strcat(line," ");
      if(start1!=-1) sprintf(line,"%s%d ",line,start1);
      for(j=strlen(line);j<linelength;j++)
	if(pos1[j]==-1) strcat(line," ");
	else if(!((pos1[j]-start1)%10)) strcat(line,":");
	else if(!((pos1[j]-start1)%5)) strcat(line,".");
	else strcat(line," ");
      printf("%s\n",line);
      /* seq lines*/
      *line=0;
      strncat(line,pseq1,linelength);
      printf("%s\n",line);
      *line=0;
      strncat(line,pseq2,linelength);
      printf("%s\n",line);
      /*2nd pos line*/
      *line=0;
      for(end2=pos2[j=linelength-1];end2==-1&&j>0;               end2=pos2[--j]);
      for(start2=pos2[j=0];         start2==-1&&j<linelength-1;start2=pos2[++j]);
      for(i=0;i<j;i++) strcat(line," ");
      if(start2!=-1) sprintf(line,"%s%d ",line,start2);
      for(j=strlen(line);j<linelength;j++)
	if(pos2[j]==-1) strcat(line," ");
	else if(!((pos2[j]-start2)%10)) strcat(line,":");
	else if(!((pos2[j]-start2)%5)) strcat(line,".");
	else strcat(line," ");
      printf("%s\n",line);

      /* mark seeds */
      {
	seed_t *tmpseed=seed;
	uchar *b, *a, blen;
	uchar *ext_id;
	*line=0;
	for(;tmpseed;tmpseed=tmpseed->next){
	  if(((tmpseed->pos2+tmpseed->len>start2)&&(tmpseed->pos2<=end2))){
	    int markstart=-1;
	    for(i=0;i<linelength&&i<neededlength;i++){
	      if(markstart!=-1&&pos2[i]==-1){
		strcat(line," ");
	      }
	      if((tmpseed->pos2+tmpseed->len>pos2[i])&&(tmpseed->pos2<=pos2[i])){
		if(markstart==-1) markstart=i;
		if((tmpseed->pos1+tmpseed->len>pos1[i])&&(tmpseed->pos1<=pos1[i])){
		  strcat(line,"^");
		}
		else{
		  strcat(line,"~");
		}
	      }
	    }
	    if(tmpseed->extra_type==IUPAC_EXTRA_INFO){
	      iupac_information *info=(iupac_information*) tmpseed->extra_information;
	      ext_id=ckalloc(sizeof(char)*(1+strlen(info->id)+
					   strlen(ori_names[info->ori])+7+
					   ((info->dist1||info->dist2)?20:0)));
	      strcpy(ext_id,info->id);
	      sprintf(ext_id,"%s (%s)",ext_id,ori_names[info->ori]);
	      if(info->dist1||info->dist2) 
		sprintf(ext_id,"%s dist:%d|%d",ext_id,info->dist1,info->dist2);
	    }
	    else if(tmpseed->extra_type==WEIGHT_EXTRA_INFO){
	      weight_information *info=(weight_information*) tmpseed->extra_information;
	      ext_id=ckalloc(sizeof(char)*(1+strlen(info->id)+
					   strlen(ori_names[info->ori])+7+
					   10+2*info->gran_width));
	      strcpy(ext_id,info->id);
	      sprintf(ext_id,"%s (%s)",ext_id,ori_names[info->ori]);
	      sprintf(ext_id,"%s score:%.*f|%.*f",ext_id,
		      info->gran_width,info->granularity*info->score1,
		      info->gran_width,info->granularity*info->score2);
	    }
	    else ext_id="Unknown";

	    if((blen=strlen(ext_id))<markstart) 
	      b=ext_id,a="";
	    else 
	      a=ext_id,b="",blen=0;
	    for(k=blen;k<markstart;k++) printf(" ");
	    printf("%s",b);
	    printf("%s",line);
	    printf("%s\n",a);
	    free(ext_id);
	  }
	  free(line);
	  line =ckcalloc(linelength+20,sizeof(char));
	}
      }
      printf("\n");
      /* cleaning up*/
      free(line);
      neededlength-=linelength;
      if(neededlength>0){
	memmove(pos1, pos1 +linelength,  neededlength*sizeof(int));
	memmove(pos2, pos2 +linelength,  neededlength*sizeof(int));
	memmove(pseq1,pseq1+linelength,1+neededlength);
	memmove(pseq2,pseq2+linelength,1+neededlength);
      }
    }
  }
  free(pos1);
  free(pos2);
  free(pseq1);
  free(pseq2);
}

#ifdef TESTPRETTYPRINT

int main(int argc, char **argv){
  int i, len1, len2, seedsnumber, start1, start2;
  uchar *seq1, *seq2, *ali1, *ali2, *poss1, *poss2, **seeds;
  edit_script_t *es;
  seed_t *seed, *seedtmp, *seedtemp;
  uchar** pseudoid;

  pseudoid=ckalloc(5*sizeof(char*));
  
  es = edit_script_new();
  for(i=len1=len2=0;i<13;i++){
    int op = random()%3 +1;
    int len = random()%10 +11;
    edit_script_more(es,op,len);
    if(op%2) len2+=len; /*INS or REP*/
    if(op>>1) len1+=len; /*DEL or REP*/
  }
  seq1=ckalloc(len1);
  seq2=ckalloc(len2);
  for(i=0;i<len1;i++) seq1[i]=random()%26+'A';
  for(i=0;i<len2;i++) seq2[i]=random()%26+'A';
  seq1[len1]=0;
  seq2[len2]=0;

  seed=ckalloc(sizeof(seed_t));
  seedtmp=seed;
  for(i=0;i<5;i++){
    seedtmp->pos1=seedtmp->pos2=random()%((len1<len2?len1:len2)-10);
    seedtmp->dist=seedtmp->ori=0;
    seedtmp->len=random()%8+1;
    pseudoid[i]=ckalloc(30);
    sprintf(pseudoid[i],"\"pseudo%d{len=%2d,pos=%3d}\"",i,seedtmp->len,seedtmp->pos1);
    seedtmp->id=pseudoid[i];
    seedtemp=seedtmp;
    seedtmp->next=ckalloc(sizeof(seed_t));
    seedtmp=seedtmp->next;
  }
  free(seedtmp);
  seedtmp=NULL;
  seedtemp->next=NULL;
  seedtemp=NULL;


  printf("\nStupid alignment of Random sequences:\n\nseq1=\"%s\"\nseq2=\"%s\"\n\nAlignment:\n",
	 seq1,seq2);
  edit_script_prnt(es);

  if(argc<=1)start1=start2=1;
  else if(argc==2)start1=start2=atoi(argv[1]);
  else start1=atoi(argv[1]),start2=atoi(argv[2]);
  edit_script_prettyprint(seq1, seq2, start1, start2, es, seed, 60);
  exit(0);
}



#endif



