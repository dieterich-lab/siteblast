static const char rcsid[]=
"$Id: bz_main.c,v 4.0 2004/08/26 09:49:25 mmichael Exp mmichael $";

#include "bz_all.h"
#include <signal.h>
#include <math.h>
#ifdef TIME
#include <time.h>
#endif
enum
{ DEFAULT_E =  30
, DEFAULT_O = 400
, DEFAULT_W =   8
, DEFAULT_R =   0
, DEFAULT_G =   0
};
#define DEFAULT_pValue 1e-3
#ifdef TIME
time_t *starttime;
#endif
static const uchar nchars[] = "ACGT";
#define CLEN(s) (sizeof((s))-1)

static const char Usage[] =
"%s seq1 [seq2]\n\n\
        t(0) type of seq:  0: DNA 1: protein\n\
	m(80M) bytes of space for trace-back information\n\
	v(0) 0: quiet; 1: verbose progress reports to stderr\n\
        A(1) Output 1:compressed; 2:prettyprinted; 3:both\n\
	B(2) 0: single strand; >0: both strands\n\
        C(0) 0: no chaining; 1: just output chain; 2: chain and extend;\n\
		3: just output HSPs\n\
	E(30) gap-extension penalty.\n\
        G(0) diagonal chaining penalty.\n\
        H(0) interpolate between alignments at threshold K = argument.\n\
	K(3000) threshold for MSPs\n\
	L(K) threshold for gapped alignments\n\
	M(50) mask threshold for seq1, if a bp is hit this many times\n\
	N(1) 0: entropy not used; 1: entropy used; >1 entropy with feedback.\n\
	O(400) gap-open penalty.\n\
	Q load the scoring matrix from a file.\n\
        R(0) antidiagonal chaining penalty.\n\
        I <file> file with patterns for seeds (IUPAC Code). Search by key word tree.\n\
        i <file> file with patterns for seeds (IUPAC Code). Search trivial (one by one).\n\
        D(0) Distance for IUPAC patterns.\n\
        P <file> file with sequence profiles for seeds.\n\
        pValue(1e-3) for threshold determination (for Profile)\n\
          powerLimit(undef) if power<powerLimit ignore Profile\n\
        power(undef) for threshold determination (for Profile)\n\
          pValueLimit(undef) if pValue>pValueLimit ignore Profile\n\
        o(15) orientation for seeds: o&1 look for sense patterns\n\
                                     o&2 lock for reverse patterns\n\
                                     o&4 lock for complement patterns\n\
                                     o&8 lock for rev. compl. patterns\n\
        T(1) 0: W-bp words;  1: 12of19;  2: 12of19 without transitions.\n\
	W(8) word size.\n\
	Y(O+300E) X-drop parameter for gapped extension.\n\
        X(..) .....\n";

static bz_flags_t bz_flags;
static int Connect(msp_t *q, msp_t *p, int scale);
align_t *inner(align_t *a, SEQ *sf1, SEQ *sf2, gap_scores_t gs, ss_t ss,
  ss_t sss, int Y, int innerK, TBack tback, bz_flags_t bz_flags,
  connect_t Connect, int MSP_Scale);
static int verbose;

// static ss_t ss;
// static ss_t sss;
static int ss[NACHARS][NACHARS];
static int sss[NACHARS][NACHARS];

enum { MSP_Scale = 100 };

static void mkmask_ss(ss_t ss, ss_t sss)
{
	int i, j, bad = ss['A']['X'];

	for (i = 0; i < NACHARS; ++i)
		for (j = 0; j < NACHARS; ++j)
		   sss[i][j] = ((isupper(i) && isupper(j)) ? ss[i][j] : bad);
}
/*get next score line*/
nextscores(FILE *fp, ss_t ss, int *gapopen, int *gapext, int *K, int *L)
{
	char buf[1024];
  int i,j;
  long t;
  float ms;	
  int bad, a, b, A, B;
  int sumofmatch;
  sumofmatch=0;
	if (fgets(buf, sizeof buf, fp) == 0)
		fatal("cannot read score file");

	for (i = 0; i < NACHARS; ++i)
		for (j = 0; j < NACHARS; ++j)
			ss[i][j] = -100;

	for (i = 0; i < CLEN(nchars); ++i) {
		A = nchars[i];
		a = tolower(A);
		for (j = 0; j < CLEN(nchars); ++j) {
		  if(fscanf(fp, "%f", &ms)==0){fatal("CANNOT read score");}
			B = nchars[j];
			b = tolower(B);
		    ss[A][B] = ss[a][B] = ss[A][b] = ss[a][b] = t = 100 * ms;
		    if(i==j){sumofmatch+=(t/4);}

		    /*printf("Current match %d\n",t);*/
		  }
		}

	for (i = 0; i < NACHARS; ++i)
		ss['N'][i] = ss[i]['N'] = ss['n'][i] = ss[i]['n'] = -1000;
	/*printf("Sum of match %d\t\n",sumofmatch*11);	*/
	*gapopen=sumofmatch*11;
	*K=*L=sumofmatch*11;
	*gapext=sumofmatch*0.5;
	mkmask_ss(ss, sss);
}

/* Compute score matrix according to REV model*/

generate_matrix(ss_t ss, double (*freq1)[4],double (*freq2)[4], double pam_dist, int *gapopen, int *gapext, int *K, int *L)
{
        double sc[4][4];
	double freq[4];
	int i,j,k;
	int c1,c2,c3,c4;
	double rate_a_g, rate_g_a, rate_tv;
	double a,c,d,f,b,e;
	double factor, sum;
	double sumofmatch;

	sumofmatch = 0;

	/*printf("PAM: %f\n",pam_dist);*/

	for (i = 0; i < NACHARS; ++i)
		for (j = 0; j < NACHARS; ++j)
			ss[i][j] = -100;

	for (k=0;k<4;k++){freq[k]=(((*freq1)[k] + (*freq2)[k])/2); /*printf("%f\n",freq[k]);*/}
	
	rate_a_g=3;  // transition A:T -> G:C
	rate_g_a=5;	// transition G:C -> A:T
	rate_tv=1;	// transversion

	pam_dist/=100;		// rescale pam

	// set up probability matrix
	a=((rate_tv/freq[0])+(rate_tv/freq[1]))/2;
	c=((rate_tv/freq[0])+(rate_tv/freq[3]))/2;
	d=((rate_tv/freq[1])+(rate_tv/freq[2]))/2;
        f=((rate_tv/freq[2])+(rate_tv/freq[3]))/2;

	b=((rate_g_a/freq[0])+(rate_a_g/freq[2]))/2;
	e=((rate_a_g/freq[1])+(rate_g_a/freq[3]))/2;
	
	sc[1][0]=a*freq[0];
	sc[2][0]=b*freq[0];
	sc[3][0]=c*freq[0];
	sc[2][1]=d*freq[1];
	sc[3][1]=e*freq[1];
	sc[3][2]=f*freq[2];

	sc[0][1]=a*freq[1];
	sc[0][2]=b*freq[2];
	sc[1][2]=d*freq[2];
	sc[0][3]=c*freq[3];
	sc[1][3]=e*freq[3];
	sc[2][3]=f*freq[3];

	sc[0][0]=0;
	sc[1][1]=0;
	sc[2][2]=0;
	sc[3][3]=0;
	
	// scale  matrix to the right pam distance
	sum=0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			sum+=freq[i]*sc[i][j];

	factor=pam_dist/sum;

	
	for (i=0;i<4;i++){
		sc[i][i]=1;
		for (j=0;j<4;j++){
			if (i!=j){
				sc[i][j]*=factor;
				sc[i][i]-=sc[i][j];  // get diagonal right
			}
		}
	}	
	
	// evaluate log odds scores
	for (i=0;i<4;i++){
	  c1=nchars[i];
		for (j=0;j<4;j++){
		  c2=nchars[j];
		  
			ss[tolower(c1)][tolower(c2)] = ss[c1][c2] = (int) 100*log(sc[i][j]/(freq[j]+1e-9));
			//printf("SCORE: %c to %c is %d\n",c1,c2,ss[c1][c2]);
			if(i==j){sumofmatch+= (double)(ss[nchars[i]][nchars[j]] / 4);}		
		}
	}	

	*gapopen= (int)(sumofmatch*11);
	*K=*L= (int)(sumofmatch*11);
	*gapext= (int)(sumofmatch*0.5);
	mkmask_ss(ss, sss);
	//printf("Sum of match %f\t\n",sumofmatch);	

}

/* connect - compute penalty for connecting two fragements */
static int Connect(msp_t *q, msp_t *p, int scale)
{
	int q_xend, q_yend, diag_diff, substitutions;

	if (p->pos1 <= q->pos1 || p->pos2 <= q->pos2)
		fatal("HSPs improperly ordered for chaining.");
	q_xend = q->pos1 + q->len - 1;
	q_yend = q->pos2 + q->len - 1;
	diag_diff = (p->pos1 - p->pos2) - (q->pos1 - q->pos2);
	if (diag_diff >= 0)
		substitutions = p->pos2 - q_yend - 1;
	else {
		substitutions = p->pos1 - q_xend - 1;
		diag_diff = -diag_diff;
	}
	return diag_diff*bz_flags.G + (substitutions >= 0 ?
		substitutions*bz_flags.R :
		-substitutions*scale*ss[(uchar)'A'][(uchar)'A']);
}

/* strand1 -- process one sequence from the second file in one orientation */
static void strand1(SEQ *sf1, uchar *rf1, SEQ *sf2, blast_table_t *bt,
		    gap_scores_t gs, ss_t ss, ss_t sss, bz_flags_t *bf,
		    TBack tback, census_t census[], bool simple_2nd_search){
  static msp_table_t *mt = 0;
  int C = bf->C;
  int T = bf->T;
  int X = bf->X;
  int Y = bf->Y;
  int K = bf->K;
  int L = bf->L;
  int N = bf->N;
  int innerK = bf->H;

  if (mt == 0) mt = msp_new_table(); // XXX - not reentrant


  mt->num = 0;
  if(bf->I||bf->i||bf->P){
    init_seedList(SEQ_LEN(sf1),SEQ_LEN(sf2));
    if(simple_2nd_search){
      mt = process_siteblast_hits(sf1, sf2, bt,
				  bf->I?blast_table_type_i:(bf->i?blast_table_type_I:
							    blast_table_type_p),
				  mt, sss, X, K, N, T==1);
    }
    else{
      mt = siteblast_search(sf1, sf2, bt,
			    bf->I?blast_table_type_i:(bf->i?blast_table_type_I:
						      blast_table_type_p),
			    mt, sss, X, K, N, T==1);
    }
  }
  else {
    mt = (T?blast_1219_search:blast_search)(sf1, sf2, bt, mt, sss, X, K, N, T==1);
  }

#ifdef TIME
  printf("Search done after %5.1f seconds.\n",difftime(time(NULL),*starttime));
#endif

  if (C == 1 || C == 2) {
    msp_t *p;
    int f = 0;
	    
    (void)msp_make_chain(mt, bz_flags.G, bz_flags.G, MSP_Scale, Connect);
    for (p = MSP_TAB_FIRST(mt); MSP_TAB_MORE(mt,p); p = MSP_TAB_NEXT(p)) {
      p->filter = 1 - p->cum_score;
      f |= p->filter;
    }

    msp_compress(mt);
    msp_sort_pos1(mt);
  }
  if (C == 1 || C > 2) {
    print_align_header(sf1, sf2, ss, &gs, K, L);
    print_msp_list(sf1, sf2, mt);
  } else if (mt != 0 && MSP_TAB_NUM(mt) != 0) {
    align_t *a;
    a = bz_extend_msps(sf1, rf1, sf2, mt, &gs, ss, Y, L, tback);
    /* next two lines added Aug. 16, 2002 */
    if (a && innerK)
      a = inner(a, sf1, sf2, gs, ss, sss, Y, innerK, tback,
		bz_flags, Connect, MSP_Scale);
    if (a) {
      int n;
      print_align_header(sf1, sf2, ss, &gs, K, L);
      print_align_list(a);
      n = census_mask_align(a, SEQ_LEN(sf1), SEQ_CHARS(sf1), rf1, census, bf->M);
      printf("x {\n  n %d\n}\n", n);
    }
    free_align_list(a);
  }
  if(bf->I||bf->i||bf->P)
    free_seedList();
}

static void bye(int sig)
{
	exit(sig);
}

void reverse_inplace(uchar *s, int len)
{
        uchar *p = s + len - 1;

        while (s <= p) {
                register uchar t = *p;
                *p-- = *s;
                *s++ = t;
	}
}



#ifndef HIDE_bz_main_main
int main(int argc, char *argv[])
{
	SEQ *sf1, *sf2;
	uchar *rf1;
	blast_table_t *bt;
	gap_scores_t gs;
	int flag_strand, flag_size, flag_census, flag_reverse;
	char *scorefile;
	TBack tback;
	census_t *census;
	FILE *fp;
	double pi[2][4];
	// flush stdio buffers when killed
	signal(SIGHUP, bye);
	signal(SIGINT, bye);
	signal(SIGTERM, bye);


#ifdef TIME
	starttime=ckalloc(sizeof(time_t));
	time(starttime);
#endif
	if (argc < 3)
		fatalf(Usage, argv[0]);
	ckargs("ABCDEGHIiKLMNOPQRTWYXbcmortv", argc, argv, 2);

	get_argval_nonneg('t', &(bz_flags.t), 0);

	if (get_cargval('Q', &scorefile))
	  {
	  if(scorefile[0]!='@' && scorefile[0]!='~')
	    {
	  scores_from_file(scorefile, ss);
	      }
	  }
	else{
	  if(bz_flags.t){
	    fatal("For protein seq. you must specify a score file Q. Exiting\n");
	  }
	  else
	  {
	    DNA_scores(ss);
	  }
	}
	mkmask_ss(ss, sss);

	get_argval_pos('m', &flag_size, 80*1024*1024);
	get_argval_pos('b', &flag_reverse, 0);
	get_argval_pos('c', &flag_census, 0);
	get_argval_pos('v', &verbose, 0);
	get_argval_nonneg('B', &flag_strand, bz_flags.t?0:2);
	get_argval_nonneg('C', &(bz_flags.C), 0);
	get_argval_nonneg('D', &(bz_flags.D), 0);
	get_argval_nonneg('E', &(gs.E), DEFAULT_E);
	get_argval_nonneg('G', &(bz_flags.G), DEFAULT_G);
	get_argval_nonneg('H', &(bz_flags.H), 0);
	get_argval_pos('K', &(bz_flags.K), 3000);
	get_argval_pos('L', &(bz_flags.L), bz_flags.K);
	get_argval_nonneg('M', &(bz_flags.M), 50);
	get_argval_nonneg('O', &(gs.O), DEFAULT_O);
	get_argval_nonneg('N', &(bz_flags.N), 1);
	get_argval_nonneg('R', &(bz_flags.R), DEFAULT_R);
	get_argfile('I', &(bz_flags.I), NULL);
	get_argfile('i', &(bz_flags.i), NULL);
	get_argfile('P', &(bz_flags.P), NULL);
	get_argval_nonneg('T', &(bz_flags.T), (bz_flags.I||bz_flags.i||bz_flags.P)?0:1);
	get_argval_pos('W', &(bz_flags.W), (bz_flags.I||bz_flags.i||bz_flags.P)?0:DEFAULT_W);
	get_argval_pos('X', &(bz_flags.X), 10*ss[(uchar)'A'][(uchar)'A']);
	get_argval_pos('Y', &(bz_flags.Y), (int)(gs.O+300*gs.E));
	get_argval_pos('o', &(bz_flags.o), (bz_flags.I||bz_flags.i||bz_flags.P)?(bz_flags.t?3:15):0);
	get_long_fargval_limit("pValue", 0., .1, &(bz_flags.pValue), DEFAULT_pUNDEF);
	get_long_fargval_limit("power", 0., 1., &(bz_flags.power), DEFAULT_pUNDEF);
	get_long_fargval_limit("pValueLimit", 0., .1, &(bz_flags.pValueLimit), DEFAULT_pUNDEF);
	get_long_fargval_limit("powerLimit", 0., 1., &(bz_flags.powerLimit), DEFAULT_pUNDEF);
	if(bz_flags.power==DEFAULT_pUNDEF&&bz_flags.pValue==DEFAULT_pUNDEF&&bz_flags.P)
	  bz_flags.pValue=DEFAULT_pValue;

	if(bz_flags.t&&!(bz_flags.I||bz_flags.i))
	  fatal("Please process protein sequences (t=1) only with patterns (i or I). Exiting\n");

	if(bz_flags.t&&(bz_flags.o&12))
	  fatal("For protein seq. use only o=1, 2 or 3! Cannot build complement. Exiting\n");

	if(bz_flags.t&&flag_strand)
	  fatal("For protein seq. use only B=0! Cannot build complement. Exiting\n");

	if((bz_flags.power!=DEFAULT_pUNDEF||bz_flags.pValue!=DEFAULT_pUNDEF)&&
	   !bz_flags.P)
	  printf("Parameter power and pValue are only used together with P and are otherwise ignored.\n");
	if(bz_flags.power!=DEFAULT_pUNDEF && bz_flags.pValue!=DEFAULT_pUNDEF)
	    fatal("Please specify either pValue or power!\n"
		  "You can power together with pValuelimit or pValue with powerlimit. Exiting\n");
	if(bz_flags.power==DEFAULT_pUNDEF&&bz_flags.pValueLimit!=DEFAULT_pUNDEF)
	  printf("pValueLimit can only be used with power. Ignoring pValueLimit!\n");
	if(bz_flags.pValue==DEFAULT_pUNDEF&&bz_flags.powerLimit!=DEFAULT_pUNDEF)
	  printf("powerLimit can only be used with pValue. Ignoring powerLimit!\n");
	if(bz_flags.T&&(bz_flags.I||bz_flags.i||bz_flags.P))
	  printf("Using T and I, i or P is not supported. Ignoring T!\n"),
	    bz_flags.T=0;
	if(bz_flags.W&&(bz_flags.I||bz_flags.i||bz_flags.P))
	  printf("W is ignored. Instead we use the length of each seed in \"%s\".\n",
		 get_argfilename(bz_flags.I?'I':(bz_flags.i?'i':'P')));
	if(bz_flags.D&&!(bz_flags.I||bz_flags.i))
	  printf("Parameter D is only used together with I or i and is otherwise ignored.\n");

	if(bz_flags.o&&!(bz_flags.I||bz_flags.i||bz_flags.P))
	  printf("Parameter o is only used together with I, i or P and is otherwise ignored.\n");

	if (bz_flags.T) bz_flags.W = 12;

        sf1 = seq_open(argv[1]);
        if (!seq_read(sf1))
		fatalf("Cannot read sequence from %s.", argv[1]);
	if(bz_flags.t){
	  if (!is_DNA(SEQ_CHARS(sf1), SEQ_LEN(sf1)))
	    fatal("The first sequence is not a DNA sequence.");
	}
	else
	  if (!is_Protein(SEQ_CHARS(sf1), SEQ_LEN(sf1)))
	    fatal("The first sequence is not a Protein sequence.");
	if (flag_reverse & 01)
	        reverse_inplace(SEQ_CHARS(sf1), SEQ_LEN(sf1));

	sf2 = seq_open(argv[2]);
	if (!seq_read(sf2))
		fatalf("Cannot read sequence from %s.", argv[2]);
	if(bz_flags.t){
	  if (!is_DNA(SEQ_CHARS(sf2), SEQ_LEN(sf2)))
	    fatal("The second sequence is not a DNA sequence.");
	}
	else
	  if (!is_Protein(SEQ_CHARS(sf2), SEQ_LEN(sf2)))
	    fatal("The second sequence is not a Protein sequence.");
	if (flag_reverse & 02)
	        reverse_inplace(SEQ_CHARS(sf2), SEQ_LEN(sf2));


#ifdef TIME
	printf("Init done after %5.1f seconds.\n",difftime(time(NULL),*starttime));
#endif

	print_job_header(ss, &gs, bz_flags.K, bz_flags.L, bz_flags.M, bz_flags.t);

	if(get_cargval('Q', &scorefile) && scorefile[0]=='@')
	  {fp = ckopen(&scorefile[0], "r");}

	do {
	  {
/* Calc. GCcontent */

	    int all1, all2;
	    long i,gc;
	    long seqLength1,seqLength2;

	    for(i=0;i<2;i++){
	      for(gc=0;gc<4;gc++){pi[i][gc]=0;}}

	    seqLength1=bz_flags.seqLength1=SEQ_LEN(sf1);
	    seqLength2=bz_flags.seqLength2=SEQ_LEN(sf2);
	    
	    for(i=gc=0;i<SEQ_LEN(sf1);i++){
	      if (SEQ_CHARS(sf1)[i]=='G'||SEQ_CHARS(sf1)[i]=='C') gc++;
	      switch (SEQ_CHARS(sf1)[i]) {
	      case 'A':
		pi[0][0]++;
		break;
	      case 'C':
		pi[0][1]++;
		break;
	      case 'G':
		pi[0][2]++;
		break;
	      case 'T':
		pi[0][3]++;
		break;
	      default:
		seqLength1--;
	      }
	      
	    }

	    bz_flags.GCcontent1=(double)gc/(seqLength1);
	    for(i=0;i<4;i++)
	      {/*printf("PRE: %f",pi[0][i]);*/pi[0][i]= (double) pi[0][i]/seqLength1; /*printf(" pi %d %f %d %d %f\n",i,pi[0][i],gc,seqLength1,bz_flags.GCcontent1);*/}

	    for(i=gc=0;i<SEQ_LEN(sf2);i++){
	      if (SEQ_CHARS(sf2)[i]=='G'||SEQ_CHARS(sf2)[i]=='C') gc++;
	      switch (SEQ_CHARS(sf2)[i]) {
	      case 'A':
		pi[1][0]++;
		break;
	      case 'C':
		pi[1][1]++;
		break;
	      case 'G':
		pi[1][2]++;
		break;
	      case 'T':
		pi[1][3]++;
		break;
	      default:
		seqLength2--;
	      }
	    }
	    
	    bz_flags.GCcontent2=(double)gc/(seqLength2);
	    for(i=0;i<4;i++)
	      {/*printf("PRE: %f",pi[1][i]);*/ pi[1][i]= (double) pi[1][i]/seqLength2; /*printf(" pi %d %f %d %d %f\n",i,pi[1][i],gc,seqLength2,bz_flags.GCcontent2);*/ }
	    
	  }
	if(get_cargval('Q', &scorefile) && scorefile[0]=='@')
	  {nextscores(fp,ss,&(gs.O),&(gs.E),&(bz_flags.K),&(bz_flags.L));}
	else if(get_cargval('Q', &scorefile) && scorefile[0]=='~')
	  {generate_matrix(ss, &pi[0], &pi[1], (double) atoi(&scorefile[1]),&(gs.O),&(gs.E),&(bz_flags.K),&(bz_flags.L));}
	  
	  if(bz_flags.I||bz_flags.i||bz_flags.P){
	    bt = siteblast_table_new(sf1, sf1->slen>sf2->slen?sf1->slen:sf2->slen, &bz_flags);
	    /*printf("%s:%d\n",__FILE__,__LINE__); dprint*/
	  }
	  else{
	    bt = (bz_flags.T ? blast_1219_new : blast_table_new)(sf1, bz_flags.W);
	  }
#ifdef TIME
	printf("Blast table build after %5.1f seconds.\n",difftime(time(NULL),*starttime));
#endif

	  if (bz_flags.C == 0 || bz_flags.C == 2) {
        	tback = ckalloc(sizeof(TB));
        	tback->space = ckalloc(flag_size*sizeof(uchar));
        	tback->len = flag_size;
	  } else
		tback = NULL;

	  census = new_census(SEQ_LEN(sf1));

	  rf1 = reverse_seq(SEQ_CHARS(sf1), SEQ_LEN(sf1));
	  strand1(sf1, rf1, sf2, bt, gs, ss, sss, &bz_flags, tback, census, 0);

	  if (flag_strand > 0) {
	    seq_revcomp_inplace(sf2);

	    if(bz_flags.i||bz_flags.I||bz_flags.P){
	      if( ((bz_flags.o&1)==1)==((bz_flags.o&8)==8) &&
		  ((bz_flags.o&2)==2)==((bz_flags.o&4)==4) ){
		/*dprint
		  printf("Can do simple search in rev. cmpl. seq!\n");*/
		correct_bt_for_2nd_search(SEQ_LEN(sf2),bt,
					  bz_flags.I?blast_table_type_i:(bz_flags.i?blast_table_type_I:
									 blast_table_type_p));
		strand1(sf1, rf1, sf2, bt, gs, ss, sss, &bz_flags, tback, census, 1);
	      }
	      else{
		/*dprint
		  printf("Cann't do simple search in rev. cmpl. seq. Must use the hard way! \n");*/
		correct_bt_for_2nd_search(0,bt,
					  bz_flags.I?blast_table_type_i:(bz_flags.i?blast_table_type_I:
									 blast_table_type_p));
		strand1(sf1, rf1, sf2, bt, gs, ss, sss, &bz_flags, tback, census, 0);
	      }
	    }
	  }
	  fflush(stdout);
	  ckfree(rf1);rf1=NULL;

	  print_intervals(stdout, census, SEQ_LEN(sf1), bz_flags.M);
	  if (flag_census) print_census(stdout, census, SEQ_LEN(sf1), 0);
	  ckfree(census);census=NULL;
	  if(!(bz_flags.I||bz_flags.i||bz_flags.P)){
	    (bz_flags.T ? blast_1219_free : blast_table_free)(bt);
	  }
	  else{
	    if(bz_flags.I) siteblast_table_i_free(bt);
	    if(bz_flags.i) siteblast_table_I_free(bt);
	    if(bz_flags.P) siteblast_table_p_free(bt);
	  }
	  if (tback) {
	    ckfree(tback->space);
	    ckfree(tback);
	  }
	} while ((seq_read(sf1))&&(seq_read(sf2)));

	print_job_footer();

        seq_close(sf2);
	seq_close(sf1);
	if(get_cargval('Q', &scorefile) && scorefile[0]=='@')
	  {fclose(fp);}

#ifdef TIME
	  printf("All done after %5.1f seconds.\n",difftime(time(NULL),*starttime));
#endif

	return 0;
}
#endif /*ifndef HIDE_bz_main_main*/

