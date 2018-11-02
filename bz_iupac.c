#include "bz_all.h"

typedef void (*Do_add) (Recognizer r, char *pattern, char *id, int dist, 
			bz_flags_t *bz_flags);
static int pattern_id=0;

void print_iupac_information(iupac_information* info){	
  printf("mem:%p\n",info);
  printf("id:%p \"%s\"\n",info->id,info->id);
  printf("ori:%p\n",info->ori);
  printf("dist:%d|%d\n",info->dist1,info->dist2);
}

static inline void do_add(Recognizer r, char *pattern, char* id, int dist, 
		   bz_flags_t *bz_flags){
  char* mod_pattern;
  iupac_information *extra;
  int i;
  int o=bz_flags->o;
  const int len=strlen(pattern);

  /*dprint
    printf("Adding \"%s\" for \"%s\" with dist %d.\t(%s:%d)\n",pattern, id, dist,__FILE__,__LINE__);*/

#define build_extra(ORI) \
 extra=ckalloc(sizeof(iupac_information));\
 extra->id=id;\
 extra->dist1=dist;extra->dist2=-1;\
 extra->ori=ORI

  if(o&1){
    build_extra(0);
    add_ident(r,pattern,pattern_id,extra,IUPAC_EXTRA_INFO);
    /*dprint
    printf("Added pattern with extra:\t\t(%s:%d)\n",__FILE__,__LINE__);
    print_iupac_information(extra);*/
  }

  if(o&(2|4|8)){
    mod_pattern=ckalloc((len+1)*sizeof(char));
    mod_pattern[len]=pattern[len]; 
  
    if(o&2){ /*reverse*/
      for(i=0;i<len;i++) 
	mod_pattern[i]=pattern[len-i-1];
      build_extra(1);
      add_ident(r,mod_pattern,pattern_id+1,extra,IUPAC_EXTRA_INFO);
     /*dprint
      printf("Added pattern with extra:\t\t(%s:%d)\n",__FILE__,__LINE__);
      print_iupac_information(extra);*/
    }

    if(o&4){ /*complement*/
      for(i=0;i<len;i++) 
	mod_pattern[i]=dna_cmpl(pattern[i]);
      build_extra(2);
      add_ident(r,mod_pattern,pattern_id+2,extra,IUPAC_EXTRA_INFO);
     /*dprint
      printf("Added pattern with extra:\t\t(%s:%d)\n",__FILE__,__LINE__);
      print_iupac_information(extra);*/
    }

    if(o&8){ /*reverse complement*/
      for(i=0;i<len;i++) 
	mod_pattern[i]=dna_cmpl(pattern[len-i-1]);
      build_extra(3);
      add_ident(r,mod_pattern,pattern_id+3,extra,IUPAC_EXTRA_INFO);
     /*dprint
      printf("Added pattern with extra:\t\t(%s:%d)\n",__FILE__,__LINE__);
      print_iupac_information(extra);*/
    }

    free(mod_pattern);
  }
#undef build_extra
}

void mutate(char* pattern, int maxdist, int useddist, int fix, Recognizer r, 
	    Do_add f, char* id, bz_flags_t *bz_flags){
  char *tpattern;
  const int length=strlen(pattern);
  char c, cost;
  
  if(maxdist<useddist) return;

  if(fix==length){
    f(r,pattern,id,useddist,bz_flags);
    return;
  }
  tpattern=ckalloc((length+1)*sizeof(char));
  strcpy(tpattern,pattern);

  if(bz_flags->t){
    /*********************************************
     * http://virology.wisc.edu/acp/CommonRes/SingleLetterCode.html (31st Aug 2004)
     *
     *	Single, 3-letter and ambiguity codes for Amino Acids
     *	
     *	The 20 naturally occuring amino acids: 
     *	
     *	Alanine         Ala     A
     *	Cysteine        Cys     C
     *	Aspartic AciD   Asp     D     (or B)
     *	Glutamic Acid   Glu     E     (or Z)
     *	Phenylalanine   Phe     F
     *	Glycine         Gly     G
     *	Histidine       His     H
     *	Isoleucine      Ile     I
     *	Lysine          Lys     K
     *	Leucine         Leu     L
     *	Methionine      Met     M
     *	AsparagiNe      Asn     N     (or B)
     *	Proline         Pro     P       
     *	Glutamine       Gln     Q     (or Z)
     *	ARginine        Arg     R
     *	Serine          Ser     S
     *	Threonine       Thr     T
     *	Valine          Val     V
     *	Tryptophan      Trp     W
     *	TYrosine        Tyr     Y
     *	
     *	Amino acids ambiguity codes: 
     *	
     *	Asparagine/Aspartic Acid    Asx   B 
     *	Glutamine/Glutamic Acid     Glx   Z 
     *	
     *	Not assigned to amino acids: 
     *	
     *	J O U
     *	
     *********************************************/

    for(tpattern[fix]='A';tpattern[fix]<'Z';tpattern[fix]++){
      switch(tpattern[fix]){
	 case 'A':
	 case 'C':
	 case 'F':
	 case 'G':
	 case 'H':
	 case 'I':
	 case 'K':
	 case 'L':
	 case 'M':
	 case 'P':
	 case 'R':
	 case 'S':
	 case 'T':
	 case 'V':
	 case 'W':
	 case 'Y':
	    cost=((pattern[fix]==tpattern[fix])?0:1);
	    break;

	 case 'D':
	 case 'N':
	    cost=(((pattern[fix]==tpattern[fix])||pattern[fix]=='B')?0:1);
	    break;
	    
	 case 'E':
	 case 'Q':
	    cost=(((pattern[fix]==tpattern[fix])||pattern[fix]=='Z')?0:1);
	    break;

	 default:
	    cost=maxdist+1;
      }
      if(maxdist>=useddist+cost)
	mutate(tpattern,maxdist,useddist+cost,fix+1,r,f,id,bz_flags);
    }
  }
  else{

    /*********************************************
     *	IUPAC Code (http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html (10th May 2004))
     *
     *      A            A           Adenine
     *      C            C           Cytosine
     *      G            G           Guanine
     *      T            T           Thymine
     *      U            U           Uracil
     *      M          A or C
     *      R          A or G
     *      W          A or T
     *      S          C or G
     *      Y          C or T
     *      K          G or T
     *      V        A or C or G
     *      H        A or C or T
     *      D        A or G or T
     *      B        C or G or T
     *      X      G or A or T or C
     *      N      G or A or T or C
     *********************************************/

  tpattern[fix]='A';
  cost=(strchr("AMRWVHDXN",pattern[fix])?0:1);
  if(maxdist>=useddist+cost)
    mutate(tpattern,maxdist,useddist+cost,fix+1,r,f,id,bz_flags);

  tpattern[fix]='C';
  cost=(strchr("CMSYVHBXN",pattern[fix])?0:1);
  if(maxdist>=useddist+cost)
    mutate(tpattern,maxdist,useddist+cost,fix+1,r,f,id,bz_flags);
  
  tpattern[fix]='G';
  cost=(strchr("GRSKVDBXN",pattern[fix])?0:1);
  if(maxdist>=useddist+cost)
    mutate(tpattern,maxdist,useddist+cost,fix+1,r,f,id,bz_flags);

  tpattern[fix]='T';
  cost=(strchr("TWYKHDBXN",pattern[fix])?0:1);
  if(maxdist>=useddist+cost)
    mutate(tpattern,maxdist,useddist+cost,fix+1,r,f,id,bz_flags);
  }
  free(tpattern);
}

static void add_to_recognizer(Recognizer r, char* pattern, char* id, 
		       bz_flags_t *bz_flags){
#ifdef TRACK_SEARCH
  /*dprint*/
  printf("adding \"%s\" for %s.\n",pattern,id);
#endif
  mutate(pattern,bz_flags->D,0,0,r,&do_add,id,bz_flags);
  pattern_id+=4;
}

void getCons(char ***Cons, char ***Ids, int *number,bz_flags_t *bz_flags, int *highest_pattern_id, char** id_list){
  char* con = ckalloc(100*sizeof(char));
  char* con_id = ckalloc(100*sizeof(char));
  char* line = ckalloc(201*sizeof(char));
  char* lineResult;
  char **cons, **ids;

  register int i,j, num;

  pattern_id=0;

  /* Determing number of cons and length of ids */
  lineResult=fgets(line,200,bz_flags->i);
  for(j=0,num=0;
      lineResult;
      lineResult=fgets(line,200,bz_flags->i)){
    sscanf(line,"  %[A-Z]  %[^\n]",con,con_id);
    j+=strlen(con_id)+1;
    num++;
  }
  rewind(bz_flags->i);
  
  free(con_id);

  *id_list = ckcalloc(j,sizeof(char));
  cons= ckalloc(num*sizeof(char*));
  ids = ckalloc(num*sizeof(char*));
  con_id=*id_list;
  i=0;

  lineResult=fgets(line,200,bz_flags->i);
  while(lineResult){
    sscanf(line,"  %[A-Z]  %[^\n]",con,con_id);
    cons[i]=ckalloc(sizeof(char)*(strlen(con)+1));
    ids[i]=con_id;
    strcpy(cons[i],con);
    i++;
    con_id+=strlen(con_id)+1;
    lineResult=fgets(line,200,bz_flags->i);
  }
  rewind(bz_flags->i);

  for(con_id=*id_list;con_id<*id_list+j;con_id++) 
    if(*con_id=='\t') *con_id=' ';
  free(con);
  free(line);
  *Ids=ids;
  *Cons=cons;
  *number=num;

  *highest_pattern_id=num*4-1;

  /*dprint
  {
    int i;
    printf("%s:%d\n",__FILE__,__LINE__);
    printf("num=%d, *number=%d\ncons=%p, *Cons=%p\n",num, *number, cons, *Cons);
    printf("%s:%d\n",__FILE__,__LINE__);
    for(i=0;i<num;i++){
      printf("cons[%d]=%s\n",i,cons[i]);
    }
    printf("cons=*Cons...\n");
    printf("%s:%d\n",__FILE__,__LINE__);
    for(i=0;i<num;i++){
      printf("...cons[%d]=%p\n",i,cons[i]);
    }
  }*/
  /*dprintprintf("%s:%d\n",__FILE__,__LINE__);*/
}

Recognizer new_iupacrecogniser(bz_flags_t *bz_flags, int *highest_pattern_id, char** id_list){
  char* seed = ckalloc(100*sizeof(char));
  char* seed_id = ckalloc(100*sizeof(char));
  char* line = ckalloc(201*sizeof(char));
  char* lineResult;

  Recognizer r = new_recognizer("","");
  register int j;

  pattern_id=0;

  /*dprint*/
  printf("new_iupacrecogniser started\t\t(%s:%d)\n",__FILE__,__LINE__);

  /* Determing length of seed ids */
  lineResult=fgets(line,200,bz_flags->I);
  for(j=0;
      lineResult;
      lineResult=fgets(line,200,bz_flags->I)){
    sscanf(line,"  %[A-Z]  %[^\n]",seed,seed_id);
    j+=strlen(seed_id)+1;
  }
  rewind(bz_flags->I);
  
  free(seed_id);
  *id_list = ckcalloc(j,sizeof(char));
  seed_id=*id_list;

  lineResult=fgets(line,200,bz_flags->I);
  while(lineResult){
    sscanf(line,"  %[A-Z]  %[^\n]",seed,seed_id);
    add_to_recognizer(r,seed,seed_id,bz_flags);
    seed_id+=strlen(seed_id)+1;
    lineResult=fgets(line,200,bz_flags->I);
  }
  rewind(bz_flags->I);
  stop_adding(r);
#ifdef TRACK_SEARCH
  /*dprint*/
  printf("\nAdding done!\n\n");
#endif
  for(seed_id=*id_list;seed_id<*id_list+j;seed_id++) 
    if(*seed_id=='\t') *seed_id=' ';
  free(seed);
  free(line);
  *highest_pattern_id=pattern_id-1;

  /*dprint*/
  printf("new_iupacrecogniser finished\t\t(%s:%d)\n",__FILE__,__LINE__);
  return(r);
}

iupac_information *clone_iupac_extra_info(iupac_information *info){
  iupac_information *ret=ckalloc(sizeof(iupac_information));

  ret->ori=info->ori;
  ret->dist1=info->dist1;
  ret->dist2=info->dist2;
  ret->id=info->id;
  return ret;
}

