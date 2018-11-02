#define BZ_WEIGHT__C
#include "bz_all.h"
#define scale 1.
#define gran .01
#define isGC(i) (((i)&1) ^ (((i)&2)>1)) /*True if i==1 || i==2*/
#define GC(i,GCcontent)  ((isGC(i)?GCcontent:1-GCcontent)/2.)
typedef void (*Do_add) (Recognizer r, char *pattern, singlepssm_t* pssm, int score, 
			bz_flags_t *bz_flags);
static int pattern_id=0;
static int counter=0;
const char alphabet[]={'A','C','G','T'};

#if 0
static inline void do_add(Recognizer r, char *pattern, pssm_t* pssm, int score, 
		   bz_flags_t *bz_flags){
  char* mod_pattern;
  int i;
  int o=bz_flags->o;
  const int len=pssm->length;
  weight_information *extra;

#ifdef TEST_WEIGHT
  printf("Should add pattern=\"%s\" with score %f\n",
	 pattern,score*pssm->granularity);
  return;
#endif

#define build_extra(ORI) \
 extra=ckalloc(sizeof(weight_information));\
 extra->ori=ORI;\
 extra->score1=score;extra->score2=-1;\
 extra->id=pssm->id;\
 extra->granularity=pssm->granularity;\
 extra->gran_width=pssm->gran_width;\
 counter++;//count(ORI);

#define count(WITCH)\
  if(!(++counter%1000)) printf("%dth Pattern: %s\n",counter,WITCH?mod_pattern:pattern);

  if(o&1){
    build_extra(0);
    add_ident(r,pattern,pattern_id,extra,WEIGHT_EXTRA_INFO);
  }

  if(o&(2|4|8)){
    mod_pattern=ckalloc((len+1)*sizeof(char));
    mod_pattern[len]=pattern[len]; 
  
    if(o&2){ /*reverse*/
      for(i=0;i<len;i++) 
	mod_pattern[i]=pattern[len-i-1];
      build_extra(1);
      add_ident(r,mod_pattern,pattern_id+1,extra,WEIGHT_EXTRA_INFO);
    }

    if(o&4){ /*complement*/
      for(i=0;i<len;i++) 
	mod_pattern[i]=dna_cmpl(pattern[i]);
      build_extra(2);
      add_ident(r,mod_pattern,pattern_id+2,extra,WEIGHT_EXTRA_INFO);
    }

    if(o&8){ /*reverse complement*/
      for(i=0;i<len;i++) 
	mod_pattern[i]=dna_cmpl(pattern[len-i-1]);
      build_extra(3);
      add_ident(r,mod_pattern,pattern_id+3,extra,WEIGHT_EXTRA_INFO);
    }

    free(mod_pattern);
  }
#undef build_extra
}
#endif

void printpssm(pssm_t *pssm){/*dprint*/
  for(;pssm;pssm=pssm->next){
    int i,j;
    printf("pssm(%p): next=(%p) threshold=%lg|%lg id=\"%s\"\n1st:\n",
	   pssm,pssm->next,pssm->t1,pssm->t2,pssm->id);
    for(i=0;i<pssm->length;i++){
      for(j=0;j<4;j++) printf("  %5d",pssm->rows1[i][j]);
      printf("      %7d     ",pssm->max_remain_score1[i]);
      for(j=0;j<4;j++) printf(" %1d",3&(pssm->order1[i]>>(2*j)));
      printf("\n");
    }
    printf("\n2nd:\n");
    for(i=0;i<pssm->length;i++){
      for(j=0;j<4;j++) printf("  %5d",pssm->rows2[i][j]);
      printf("      %7d     ",pssm->max_remain_score2[i]);
      for(j=0;j<4;j++) printf(" %1d",3&(pssm->order2[i]>>(2*j)));
      printf("\n");
    }
    printf("\n");
  }
}

void printsinglepssm(pssm_t *pssm){/*dprint*/
  pssm_t *tmp;

  if(pssm){
    tmp=pssm->next;
    pssm->next=NULL;
    printpssm(pssm);
    pssm->next=tmp;
  }
}

void printprofile(profile_t *p){/*dprint*/
  for(;p;p=p->next){
    int i,j;
    printf("profile(%p): next=(%p) id=\"%s\"\n",p,p->next,p->id);
    for(i=0;i<p->length;i++){
      for(j=0;j<4;j++) printf("  %12f",p->rows[i][j]);
      printf("\n");
    }
    printf("\n");
  }
}

#if 0
void genarate_from_profile(char* pattern, int pos, int current_score,
			   singlepssm_t* pssm, Recognizer r, Do_add f,
			   bz_flags_t *bz_flags){
  int i,j;

  if(!pssm) 
    fatalf("Called genarate_from_profile with pssm eq. NULL\t(%s:%d)\n",
	   __FILE__,__LINE__);

  if(pos==pssm->length){
    f(r, pattern, pssm, current_score, bz_flags);
    return;
  }

  /*dprint printf("GCcontent: %f\n",bz_flags->GCcontent);*/
  if(!(pssm->t==pssm->t)){
    fatalf("There is no correct threshold. (internal error)\n Exiting.\t\t(%s:%l)\n",
	   __FILE__,__LINE__);
  }
  /*dprint
  printsinglepssm(pssm);
  printf("%s:%d\n",__FILE__,__LINE__);*/
  for(j=(pssm->order[pos])&3,i=0;
      i<4&&current_score+pssm->rows[pos][j]+pssm->max_remain_score[pos]
	>= 
	pssm->t/pssm->granularity;
      j=(pssm->order[pos]>>(2*++i))&3){
#ifdef TEST_WEIGHT   
    printf("pos=%d, j=%d, i=%d, pssm->order[%d]=#%X\n",pos,j,i,pos,pssm->order[pos]);
#endif
    pattern[pos]=alphabet[j];
    genarate_from_profile(pattern,pos+1,current_score+pssm->rows[pos][j],
			  pssm,	r, f, bz_flags);
  }
}
 
static void add_to_recognizer(Recognizer r, pssm_t* pssm, bz_flags_t *bz_flags){
  int counterold;
  for(;pssm;pssm=pssm->next){
    char *pattern = ckcalloc(1+pssm->length,sizeof(char));
#ifdef TRACK_SEARCH
    /*dprint*/
    printf("Genarating for PSSM with id:%s \t(%s:%d)\n",
	   pssm->id,__FILE__,__LINE__);
#endif
    counterold=counter;
    genarate_from_profile(pattern,0,0,pssm,r,&do_add,bz_flags);
#ifdef TRACK_SEARCH
    /*dprint*/
    printf("Genarated %d patterns (in total)\n",counter);
    printf("++ %d patterns for this (length=%d).\n",counter-counterold,pssm->length);
#endif
    free(pattern);
    pattern_id+=4;
  }
}
#endif

static void order (pssm_t *pssm){
  for(;pssm;pssm=pssm->next){
    int i,j,k,o1[4],o2[4];
    pssm->order1=ckalloc(pssm->length*sizeof(char));
    pssm->order2=ckalloc(pssm->length*sizeof(char));
    for(i=0;i<pssm->length;i++){
      for(j=0;j<4;j++) o1[j]=o2[j]=j;
      for(j=0;j<4;j++)
	for(k=j+1;k<4;k++){
	  if(pssm->rows1[i][o1[k]]>pssm->rows1[i][o1[j]]){
	    int h;
	    h=o1[j];
	    o1[j]=o1[k];
	    o1[k]=h;
	  }
	  if(pssm->rows2[i][o2[k]]>pssm->rows2[i][o2[j]]){
	    int h;
	    h=o2[j];
	    o2[j]=o2[k];
	    o2[k]=h;
	  }
	}
      pssm->order1[i]=o1[0]|(o1[1]<<2)|(o1[2]<<4)|(o1[3]<<6);
      pssm->order2[i]=o2[0]|(o2[1]<<2)|(o2[2]<<4)|(o2[3]<<6);
      /*dprint
      printf("Scores: %d %d %d %d. Ordering: %d %d %d %d\t(%s:%d)\n",
	     pssm->rows[i][0],pssm->rows[i][1],pssm->rows[i][2],pssm->rows[i][3],
	     o[0],o[1],o[2],o[3],__FILE__,__LINE__);
      */
    }
 }
}

static void calc_max_remain(pssm_t *pssm){
  for(;pssm;pssm=pssm->next){
    int i;
    pssm->max_remain_score1=ckalloc(pssm->length*sizeof(int));
    pssm->max_remain_score2=ckalloc(pssm->length*sizeof(int));
    pssm->max_remain_score1[pssm->length-1]= pssm->max_remain_score2[pssm->length-1]=0;
    for(i=pssm->length-2;i>=0;i--){
      pssm->max_remain_score1[i]=
	pssm->max_remain_score1[i+1]+
	pssm->rows1[i+1][3&pssm->order1[i+1]];
      pssm->max_remain_score2[i]=
	pssm->max_remain_score2[i+1]+
	pssm->rows2[i+1][3&pssm->order2[i+1]];
    }
    pssm->rev_max_remain_score1=ckalloc(pssm->length*sizeof(int));
    pssm->rev_max_remain_score2=ckalloc(pssm->length*sizeof(int));
    pssm->rev_max_remain_score1[pssm->length-1]=pssm->rev_max_remain_score2[pssm->length-1]=0;
    for(i=pssm->length-2;i>=0;i--){
      pssm->rev_max_remain_score1[i]=
	pssm->rev_max_remain_score1[i+1]+
	pssm->rows1[pssm->length-2-i][3&pssm->order1[pssm->length-2-i]];
      pssm->rev_max_remain_score2[i]=
	pssm->rev_max_remain_score2[i+1]+
	pssm->rows2[pssm->length-2-i][3&pssm->order2[pssm->length-2-i]];
    }
  }
}


#if 0 /*For threshold calculation, we need the Profile*/
pssm_t *read_pssm_file(FILE* infile,char* id_list){
  pssm_t *pssm_base, *pssm, *lastpssm;
  int length;
  char* line = ckalloc(201*sizeof(char));
  char* lineResult;
  char* pos;
  double *scorefield;
  int id_list_length=0;
  char* id_list_pointer;
  int id_length;
  register int j;

  pssm_base=pssm=ckalloc(sizeof(pssm_t));

  lineResult=fgets(line,200,infile);
  while(lineResult&&line[0]!='>')
    lineResult=fgets(line,200,infile);
  while(lineResult){
    for(pos=line;*pos;pos++);
    for(pos--;isspace(*pos);pos--);
    for(;isdigit(*pos);pos--);
    for(;isspace(*pos);pos--);
    pos++;
    *pos='\0';
    id_length=pos-line;
    id_list_length+=id_length;
    id_list_pointer=id_list+id_list_length-id_length;
    strncpy(id_list_pointer,line+1,id_length);
    pos++;
    sscanf(pos,"%d",&length);

    pssm->id=id_list_pointer;
    pssm->length=length;
    pssm->t=0./0.; /*NaN*/
    pssm->rows=ckalloc(length*sizeof(row));
    scorefield = ckalloc(length*4*sizeof(double));

    /* parsing rows */
    for(j=0,lineResult=fgets(line,200,infile);
	lineResult&&line[0]!='<';
	j++,lineResult=fgets(line,200,infile)){
      if (j<length) sscanf(line,"%lf %lf %lf %lf",
			   &scorefield[4*j],
			   &scorefield[4*j+1],
			   &scorefield[4*j+2],
			   &scorefield[4*j+3]);
    }
    if(j!=length)
      fatalf("Last Line:\"%s\"\nIn pssm \"%s\" are %d rows instead of %d. Exiting.\t(%s:%d)\n",
	     line,id_list_pointer, j, length,__FILE__,__LINE__);

    /*parsing granularity*/
    if(lineResult&&line[0]=='<'){
      double tmp;
      sscanf(line+2,"%lf", &(pssm->granularity));
      tmp=-log(gran)/log(10)+1.-1e-8;
      pssm->gran_width=tmp;
    }
    else
      fatalf("Error while Parsing pssm file. Exiting.\t(%s:%d)\n",__FILE__,__LINE__);
	     
    for(j=0;j<length;j++){
      int i;
      pssm->rows[j][0]=(int)(scorefield[4*j]/pssm->granularity);
      pssm->rows[j][1]=(int)(scorefield[4*j+1]/pssm->granularity);
      pssm->rows[j][2]=(int)(scorefield[4*j+2]/pssm->granularity);
      pssm->rows[j][3]=(int)(scorefield[4*j+3]/pssm->granularity);
    }
    free(scorefield);scorefield=NULL;
    pssm->next=ckalloc(sizeof(pssm_t));
    lastpssm=pssm;
    pssm=pssm->next;

    lineResult=fgets(line,200,infile);
    while(lineResult&&line[0]!='>'){
      lineResult=fgets(line,200,infile);
    }
  }
  rewind(infile);
  free(line);
  free(lastpssm->next);
  lastpssm->next=NULL;
  pssm=NULL;
  order(pssm_base);
  calc_max_remain(pssm_base);
  return pssm_base;
}
#endif

profile_t *read_profile_file(FILE* infile,char* id_list){
  profile_t *profile_base, *profile, *lastprofile;
  int length;
  char* line = ckalloc(501*sizeof(char));
  char* lineResult;
  char* pos;
  int id_list_length=0;
  char* id_list_pointer;
  int id_length;
  register int j;

  profile_base=profile=ckalloc(sizeof(profile_t));

  lineResult=fgets(line,500,infile);
  while(lineResult&&line[0]!='>')
    lineResult=fgets(line,500,infile);
  while(lineResult){
    for(pos=line;*pos;pos++);
    for(pos--;isspace(*pos);pos--);
    for(;isdigit(*pos);pos--);
    for(;isspace(*pos);pos--);
    pos++;
    *pos='\0';
    id_length=pos-line;
    id_list_length+=id_length;
    id_list_pointer=id_list+id_list_length-id_length;
    strncpy(id_list_pointer,line+1,id_length);
    pos++;
    sscanf(pos,"%d",&length);

    profile->id=id_list_pointer;
    profile->length=length;
    profile->rows=ckalloc(length*sizeof(frow));

    /* parsing rows */
    for(j=0,lineResult=fgets(line,500,infile);
	lineResult&&line[0]!='<';
	j++,lineResult=fgets(line,500,infile)){
      if (j<length) sscanf(line,"%lf %lf %lf %lf",
			   &(profile->rows[j][0]),
			   &(profile->rows[j][1]),
			   &(profile->rows[j][2]),
			   &(profile->rows[j][3]));
    }
    if(j!=length)
      fatalf("Last Line:\"%s\"\nIn pssm \"%s\" are %d rows instead of %d. Exiting.\t(%s:%d)\n",
	     line,id_list_pointer, j, length,__FILE__,__LINE__);

    profile->next=ckalloc(sizeof(profile_t));
    lastprofile=profile;
    profile=profile->next;

    lineResult=fgets(line,500,infile);
    while(lineResult&&line[0]!='>'){
      lineResult=fgets(line,500,infile);
    }
  }
  rewind(infile);
  free(line);
  free(lastprofile->next);
  lastprofile->next=NULL;
  profile=NULL;
  return profile_base;
}


pssm_t *profileToPSSM(profile_t *profile_base, double GCcontent1, double GCcontent2){
  profile_t *profile;
  pssm_t *pssm_base, *pssm;
  register int i,j;
  double tmp=-log(gran)/log(10)+1.-1e-8;

  pssm_base=pssm=ckalloc(sizeof(pssm_t));
  for(profile=profile_base;profile;profile=profile->next){
    pssm->id=profile->id;
    pssm->length=profile->length;
    pssm->granularity=gran;
    pssm->gran_width=tmp;
    pssm->t1=pssm->t2=0./0.; /*NaN*/
    pssm->rows1=ckalloc(pssm->length*sizeof(row));
    pssm->rows2=ckalloc(pssm->length*sizeof(row));
    for(j=0;j<profile->length;j++)
      for(i=0;i<4;i++){
	pssm->rows1[j][i]=(int) ((scale*log(profile->rows[j][i]/GC(i,GCcontent1)))/gran);
	pssm->rows2[j][i]=(int) ((scale*log(profile->rows[j][i]/GC(i,GCcontent2)))/gran);
      }
    if(profile->next){
      pssm->next=ckalloc(sizeof(pssm_t));
      pssm=pssm->next;
    }
    else pssm->next=NULL;
  }
  order(pssm_base);
  calc_max_remain(pssm_base);
  return pssm_base;
}
#undef isGC 
#undef GC 

void free_profiles(profile_t *p){
  if(!p) return;
  free_profiles(p->next);
  ZFREE(p->rows);
  ZFREE(p);
}
void free_pssms(pssm_t *p){
  if(!p) return;
  free_pssms(p->next);
  ZFREE(p->max_remain_score1);
  ZFREE(p->max_remain_score2);
  ZFREE(p->rev_max_remain_score1);
  ZFREE(p->rev_max_remain_score2);
  ZFREE(p->rows1);
  ZFREE(p->rows2);
  ZFREE(p->order1);
  ZFREE(p->order2);
  ZFREE(p);
}

pssm_t *get_PSSM(bz_flags_t *bz_flags, char** id_list, int seq_length, int *highest_pattern_id){
  char* line = ckalloc(201*sizeof(char));
  char* lineResult;
  int id_list_length=0;
  char* pos;
  FILE *infile;
  static profile_t *profile=NULL;
  pssm_t *pssm;

  profile_t *pr;
  pssm_t *ps, *pre_ps;


  /*dprint
  printf("%s:%d\n",__FILE__, __LINE__);*/

  if(! (profile&& *id_list)){
    if(profile|| *id_list)
      fatalf("Internal Error. Exiting.\t\t(%s:%d)\n",__FILE__,__LINE__);
    else{

  /*dprint
  printf("%s:%d\n",__FILE__, __LINE__);*/

      *highest_pattern_id=0;

      /*determine length of ids*/
      if(bz_flags->P) infile=bz_flags->P;
      else fatalf("Internal Error. Exiting.\t\t(%s:%d)\n",__FILE__,__LINE__);

      lineResult=fgets(line,200,infile);
      while(lineResult&&line[0]!='>') lineResult=fgets(line,200,infile);
      while(lineResult){
	(*highest_pattern_id)++;
	for(pos=line;*pos;pos++);
	for(pos--;isspace(*pos);pos--);
	for(;isdigit(*pos);pos--);
	for(;isspace(*pos);pos--);
	pos++;
	id_list_length+=pos-line;
	lineResult=fgets(line,200,infile);
	while(lineResult&&line[0]!='>') lineResult=fgets(line,200,infile);
      }
      *highest_pattern_id*=4;
      (*highest_pattern_id)--;
      rewind(infile);

      *id_list=ckalloc(id_list_length);


      /* obtain PSSMs */
      profile=read_profile_file(bz_flags->P,*id_list);
      /* printprofile(profile_base); dprint*/
    }
  }
  ckfree(line);

  /*dprint
  printf("%s:%d\n",__FILE__, __LINE__);*/

  pssm=profileToPSSM(profile,bz_flags->GCcontent1,bz_flags->GCcontent2);

  /* dprint
  printf("GCcontent=%f\n",bz_flags->GCcontent);
  printf("%s:%d\n",__FILE__,__LINE__);
  printpssm(pssm);*/

  ps=pssm;
  pre_ps=NULL;
  pr=profile;
  while(ps&&pr){
    if(!calc_threshold(ps, pr, bz_flags)){
      if(pre_ps){
	/*skip ps*/
	pre_ps->next=ps->next;
	/*free ps*/
	ps->next=NULL;
	free_pssms(ps);
	/*set ps for next step*/
	ps=pre_ps->next;
      }
      else{
	ps=pssm->next;
	pssm->next=NULL;
	free_pssms(pssm);
	pssm=ps;
      }
    }
    else{
      pre_ps=ps;
      ps=ps->next;
    }
    pr=pr->next;
  }

  /*
  for(ps=pssm, pre_ps=NULL, pr=profile; ps&&pr; ps=(pre_ps=ps)->next, pr=pr->next)
    calculate_threshold(ps, pr, seq_length, bz_flags);
  */

  /*
    free_profiles(profile);
  */
  return pssm;
}

weight_information *clone_weight_extra_info(weight_information *info){
  weight_information *ret=ckalloc(sizeof(weight_information));
  ret->ori=info->ori;
  ret->score1=info->score1;
  ret->score2=info->score2;
  ret->id=info->id;
  ret->granularity=info->granularity;
  ret->gran_width=info->gran_width;
  return ret;
}

