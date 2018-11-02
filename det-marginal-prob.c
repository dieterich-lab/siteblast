/* Based on:
 *
 * Copyright 1996 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   Dept. of Molecular, Cellular, and Developmental Biology
 *   University of Colorado
 *   Campus Box 347
 *   Boulder, CO  80309-0347
 *
 *   hertz@boulder.colorado.edu
 */

/* Round a floating point number to the closest integer.
 * Positive numbers: >= 0.5 is rounded up.
 * Negative numbers: <= -0.5 is rounded down. */
#define ROUND_TO_INT(n) (((n) > 0.0) ? ((int)((n) + 0.5)) : ((int)((n) - 0.5)))

/* Sum numbers given their logarithms.  Subtract the larger logarithm
 * from each logarithm to avoid adding 0's or infinities within the
 * precision of the machine.  Returns the logarithm of the sum. */
#define SUM_LN(n1, n2) (((n1) >= (n2)) ? \
			((n1) + log(exp((n2) - (n1)) + 1.0)) : \
			((n2) + log(exp((n1) - (n2)) + 1.0)))
#define INFINITY 1.0e100     /* A large double. */
#define INT_INF 2147483647   /* The largest 4 byte positive integer. */
#define ERROR 1.0e-6         /* Leeway factor when determining whether
			      * two scores or two column sums are the same. */

#include "bz_all.h"
int A_size=4; /*DNA*/
double **P;    /*background dist.*/

/* Variables for determining the array of marginal probabilities
 * for integer scores. */
double *Marginal_prob;
double Max_int_range;
double Min_score;
double Multiple;
double Max_exact_score;
double Min_exact_score;
char Minus_inf;
int Max_int_score;
int Min_int_score;
double *ArrayA;
double *Ln_P1;


singlepssm_t *pssm;
bz_flags_t *bz_flags;

// #define PRINTDISTRIBUTION
// #define PRINTRANGE 1           /* Print the range of scores and p-values. */
#define T_SIZE A_size      /* The total number of symbols in the matrix. */
#define LN_PROB_SUM 0.0    /* ln(sum of all possible probabilities)
			    * Can be > 0 in the presence of gaps. */

void det_marginal_prob(double **weight_mat, int width);

/* Wrapper */
bool calc_threshold(pssm_t *Pssm, profile_t *profile, bz_flags_t *Bz_flags){
  int i,j, fst;
  double ln_p;
/* Matrix with doubles. */
  double **Matrix;

  double *Marginal_prob_back;
  int Min_int_score_back;
  double *Marginal_prob_profil;
  double power1=.0,power2=.0,pValue1=.0,pValue2=.0;
  bool ret=1;

  pssm=ckalloc(sizeof(singlepssm_t));
  bz_flags=Bz_flags;

  pssm->length=Pssm->length;
  pssm->granularity=Pssm->granularity;

  for(fst=1;fst>=0;fst--){
    pssm->max_remain_score=fst?Pssm->max_remain_score1:Pssm->max_remain_score2;
    pssm->t               =fst?Pssm->t1               :Pssm->t2;
    pssm->rows            =fst?Pssm->rows1            :Pssm->rows2;
    pssm->order           =fst?Pssm->order1           :Pssm->order2;

    P=ckalloc(4*sizeof(double*));
    Matrix=ckalloc(4*sizeof(double*));
    for(j=0;j<4;j++){
      Matrix[j]=ckalloc(Pssm->length*sizeof(double));
      P[j]=ckalloc(Pssm->length*sizeof(double));
      for(i=0;i<Pssm->length;i++){
	Matrix[j][i]=Pssm->granularity*
	  (double)((fst?Pssm->rows1:Pssm->rows2)[i][j]);
      }
    }
  
    Ln_P1=ckalloc(4*sizeof(double));
    Ln_P1[1]=Ln_P1[2]=log(P[1][0]=P[2][0]=
			  (fst?bz_flags->GCcontent1:bz_flags->GCcontent2)/2.);     /*GC*/
    Ln_P1[0]=Ln_P1[3]=log(P[0][0]=P[3][0]=
			  (1-(fst?bz_flags->GCcontent1:bz_flags->GCcontent2))/2.); /*AT*/
    for(j=0;j<4;j++){
      for(i=1;i<Pssm->length;i++){
	P[j][i]=P[j][i-1];
      }
    }

    if(fst) Max_int_range=Pssm->max_remain_score1[0]+Pssm->rows1[0][Pssm->order1[0]&3];
    else Max_int_range=Pssm->max_remain_score2[0]+Pssm->rows2[0][Pssm->order2[0]&3];

    Marginal_prob=NULL;
    Min_int_score=0;
    Min_score=0.;

/* determined later:
 *
 *	ArrayA			Multiple	Minus_inf
 *	Max_int_score		Max_exact_score	Min_exact_score
 *
 */


    if(profile->length!=Pssm->length)
      fatalf("The profile (id:%s) has length %d, but\nthe pssm    (id:%s) has length %d.\n"
	     "Exiting!\t\t(%s:%d)\n",profile->id,profile->length,Pssm->id,Pssm->length,
	     __FILE__,__LINE__);

    det_marginal_prob(Matrix, Pssm->length);
    Marginal_prob_back=Marginal_prob;
    Min_int_score_back=Min_int_score;
    Marginal_prob=NULL;
    Min_int_score=0;
    Min_score=0.;
    for(i=0;i<profile->length;i++){
      double sum;
      for(sum=j=0;j<4;j++) sum+=profile->rows[i][j];
      for(j=0;j<4;j++) P[j][i]=profile->rows[i][j]/sum;
    }
    det_marginal_prob(Matrix, Pssm->length);
    Marginal_prob_profil=Marginal_prob;

    if(bz_flags->pValue!=DEFAULT_pUNDEF){
      ln_p=log(bz_flags->pValue/(fst?bz_flags->seqLength1:bz_flags->seqLength2));
      for(i=Max_int_range-1;i>=0&&Marginal_prob_back[i]<ln_p;i--);
      fst?Pssm->t1:Pssm->t2=((double)i)/Multiple;
      fst?power1:power2=exp(Marginal_prob_profil[i]);
      fst?pValue1:pValue2=exp(Marginal_prob_back[i])*
	(fst?bz_flags->seqLength1:bz_flags->seqLength2);
    }
    else if(bz_flags->power!=DEFAULT_pUNDEF){
      ln_p=log(bz_flags->power);
      for(i=Max_int_range-1;i>=0&&Marginal_prob_profil[i]<ln_p;i--);
      fst?Pssm->t1:Pssm->t2=((double)i)/Multiple;
      fst?power1:power2=exp(Marginal_prob_profil[i]);
      fst?pValue1:pValue2=exp(Marginal_prob_back[i])*
	(fst?bz_flags->seqLength1:bz_flags->seqLength2);
    }
    free(Marginal_prob_back+Min_int_score_back);
    free(Marginal_prob_profil+Min_int_score);
    for(j=0;j<4;j++){
      free(Matrix[j]);
      free(P[j]);
    }
    free(Matrix);
    free(P);
    free(Ln_P1);
  }
  /*checking the limits*/


  if(bz_flags->pValue!=DEFAULT_pUNDEF){
    if(bz_flags->powerLimit!=DEFAULT_pUNDEF&&
       (power1<bz_flags->powerLimit||power2<bz_flags->powerLimit)){
      fprintf(stderr,"The power for Profile \"%s\" is %lf|%lf, lesser than powerLimit=%lf.\n"
	     "\t Ignoring this Profile!\n",Pssm->id,power1,power2,bz_flags->powerLimit);
      Pssm->t1=1000.*(Pssm->max_remain_score1[0]+Pssm->rows1[0][Pssm->order1[0]&3]);
      Pssm->t2=1000.*(Pssm->max_remain_score2[0]+Pssm->rows2[0][Pssm->order2[0]&3]);
      ret=0;
    }
  }
  else if(bz_flags->power!=DEFAULT_pUNDEF){
    if(bz_flags->pValueLimit!=DEFAULT_pUNDEF&&
       (pValue1>bz_flags->pValueLimit||pValue2>bz_flags->pValueLimit)){
      fprintf(stderr,"The pValue for Profile \"%s\" is %lf|%lf, greater than pValueLimit=%lf.\n"
	     "\t Ignoring this Profile!\n",Pssm->id,pValue1,pValue2,bz_flags->pValueLimit);
      Pssm->t1=1000.*(Pssm->max_remain_score1[0]+Pssm->rows1[0][Pssm->order1[0]&3]);
      Pssm->t2=1000.*(Pssm->max_remain_score2[0]+Pssm->rows2[0][Pssm->order2[0]&3]);
      ret=0;
    }
  }
  
#ifdef TRACK_SEARCH
  if(ret){
    /*dprint*/
    printf("%s: \t\t\t(%s:%d)\n",Pssm->id,__FILE__,__LINE__);
    printf("pValue=%lg(%lg)|%lg(%lg), power=%lg|%lg, threshold=%lg|%lg\n",
 	   pValue1,pValue1/bz_flags->seqLength1,
	   pValue2,pValue2/bz_flags->seqLength2,
	   power1,power2,
 	   Pssm->t1,Pssm->t2);
  }
#endif

  free(pssm);
  return ret;
}




/* Functions for determining the exact p-values of the approximate matrix.
 * The following external variables are used:
 *
 *    T_SIZE                Multiple              Max_int_score
 *    P            	    Max_exact_score 	  Min_int_score
 *    Marginal_prob	    Min_exact_score 	  ArrayA
 *    Max_int_range	    Min_pseudo_score	  Ln_P
 *    Min_score    	    Minus_inf
 */


/* Properties of the integral weight matrix. */
static int *STmax_score;   /* Maximum score for each column. */
static int *STmin_score;   /* Minimum (non -INFINITY) score for each column. */
static int STmax_int_range;/* The range: Max_int_score - Min_int_score. */

static double *STarrayA = (double *)NULL;   /* Array for holding
					     * intermediate results. */
static double *STarrayB;   /* Array for holding current column. */
static double *STarrayC;   /* Array for holding intermediate results. */
static int *STarray_index; /* Index to STarrayB[],
			    * the array holding the current column. */


/* Determine the p-values of the approximate matrix.
 * Place p-values in Marginal_prob[];
 * Multiple is the conversion factor from exact scores to integral scores. */
void det_marginal_prob(weight_mat, width)
     double **weight_mat;           /* The weight matrix. */
     int width;                     /* Width of the matrix. */
{
  void determine_variables();
  void appr_mat_space();
  void determine_distribution();
  void determine_marginal_prob();


  /* Make sure the integer range has not been set too large. */
  if (Max_int_range > (double)(INT_INF - width - 3))
    {
      fprintf(stderr, "The integer range of %g is too large,\n",Max_int_range);
      fprintf(stderr, "it must be less than or equal to %d - %d - 3\n",
	      INT_INF, width);
      exit(1);
    }

  /* Reset Marginal_prob[0] to the beginning of the allocated memory. */
  Marginal_prob += Min_int_score;

  /* Allocate space for STmax_score[] and STmin_score[]; determine Minus_inf,
   * Max_exact_score, Min_exact_score, Min_pseudo_score, Multiple,
   * STmax_score[], STmin_score[], Max_int_score, Min_int_score,
   * and STmax_int_range. */
  determine_variables(weight_mat, width);

  /* Allocate space for "STarrayA[]", "STarrayB[]", "STarrayC[]",
   * and "STarray_index[]". */
  appr_mat_space();

  /* Determine the logarithm of the distribution of the integral scores of
   * the "weight_mat[][]" and place in ArrayA[]. */
  determine_distribution(weight_mat, width);

  /* Determine the logarithm of the p-values and place in Marginal_prob[]. */
  determine_marginal_prob();


  /* Free arrays except for STarrayA. */
  ckfree(STmax_score);
  ckfree(STmin_score);
  ckfree(STarrayB);
  ckfree(STarrayC);
  ckfree(STarray_index);
}

/* Allocate space for STmax_score[] and STmin_score[]; determine Minus_inf,
 * Max_exact_score, Min_exact_score, Min_pseudo_score, Multiple, STmax_score[],
 * STmin_score[], Max_int_score, Min_int_score, and STmax_int_range. */
void determine_variables(weight_mat, width)
     double **weight_mat;     /* The weight matrix. */
     int width;               /* Width of the matrix. */
{
  int i, j;
  char max_min_equality;  /* YES: maximum and minimum (non -INFINITY)
			   *      scores are equal. */
  int min_int_score;      /* Minimum integral score greater than -INFINITY. */
  double min_prob_score;  /* Minimum score for estimating probabilities. */

  /* Allocate space for "STmax_score[]" and "STmin_score[]". */
  STmax_score = ckcalloc(width, sizeof(int));
  STmin_score = (int *)ckcalloc(width, sizeof(int));

  /* Determine "STmax_score[]", "STmin_score[]", and Max_int_score. */
  for(i=Max_int_score = min_int_score = 0;i<pssm->length;i++){
    Max_int_score+=(STmax_score[i]=pssm->rows[i][pssm->order[i]&3]);
    min_int_score+=(STmin_score[i]=pssm->rows[i][(pssm->order[i]>>6)&3]);
  }

  Multiple=1./pssm->granularity;
  Min_exact_score=pssm->granularity*(double)min_int_score;
  Max_exact_score=pssm->granularity*(double)Max_int_score;


  /* Determine the minimum score for calculating p-values. */
  if (Min_score >= Min_exact_score) min_prob_score = Min_score;
  else min_prob_score = Min_exact_score;

  /* Make sure the minimum score for calculating p-values is
   * not greater than the maximum possible score. */
  if (min_prob_score > Max_exact_score)
    {
      fprintf(stderr, "The minimum score (%f) for calculating p-values is greater",min_prob_score);
      fprintf(stderr, "\nthan the maximum possible score (%f).\n",Max_exact_score);
      fprintf(stderr, "Min_score=%f, Min_exact_score=%f, min_int_score=%d\n",
	      Min_score, Min_exact_score, min_int_score);
      exit(1);
    }
  /* Adjustment for the special case when the minimum score for
   * calculating p-values is equal to the maximum possible score. */
  else if (min_prob_score == Max_exact_score) max_min_equality = 1;
  else max_min_equality = 0;


  /* Determine "Min_int_score". */
  Min_int_score = ROUND_TO_INT(Multiple * min_prob_score);
  if ((Min_exact_score >= Min_score) && (min_int_score < Min_int_score))
    Min_int_score = min_int_score;
  --Min_int_score;
  if ((max_min_equality) || (Min_int_score >= Max_int_score))
    Min_int_score = Max_int_score - 1;

  /* Determine "STmax_int_range". */
  STmax_int_range = Max_int_score - Min_int_score;
  if (STmax_int_range < 1) fatal("BUG in determine_variables()\n");
}


/* Allocate space for "STarrayA[]", "STarrayB[]", "STarrayC[]",
 * and "STarray_index[]". */
void appr_mat_space()
{
  int i;
  int max_score_1 = STmax_int_range + 1;  /* Maximum size of arrays. */


  STarrayA = (double *)ckrealloc(STarrayA, max_score_1*sizeof(double));

  STarrayB = (double *)ckcalloc(max_score_1, sizeof(double));

  STarrayC = (double *)ckcalloc(max_score_1, sizeof(double));

  STarray_index = (int *)ckcalloc(T_SIZE, sizeof(int));

  for (i = 0; i < max_score_1; ++i)
    {
      STarrayA[i] = -INFINITY;
      STarrayC[i] = -INFINITY;
      STarrayB[i] = 0.0;
    }
}


/* Determine the logarithm of the distribution of the integral scores of
 * the "weight_mat[][]" and place in ArrayA[]. */
void determine_distribution(weight_mat, width)
     double **weight_mat;     /* The weight matrix. */
{
  int i, j;
  double ln_prob_AB;       /* ln(column_probability in STarrayB[] multiplied
			    *    by previous_probability in STarrayA[]) */
  double ln_prob_C;        /* The ln(probability currently in STarrayC[]). */
  double *temp_array;      /* Temporary pointer when exchanging pointers to
			    * "STarrayA" and "STarrayC". */
  int max_score;           /* Maximum score prior to current multiplication. */
  int new_max_score;       /* Maximum score after current multiplication. */
  int min_score;           /* Minimum score prior to current multiplication. */
  int new_min_score;       /* Minimum score after current multiplication. */
  int diff_sub;      /* Number of different substitution scores (<= A_size). */

  /* Pointers to multiplication arrays so maximum index = maximum score. */
  double *arrayA;       /* Pointer to STarrayA: product of previous columns. */
  double *arrayB;       /* Pointer to STarrayB: current column. */
  double *arrayC;       /* Pointer to STarrayC: product of current column
			 *                      with previous columns. */
  /* Minimum index in each of the multiplication arrays. */
  int min_idx_A;        /* arrayA */
  int min_idx_B;        /* arrayB */
  int min_idx_C;        /* arrayC */

  /* The score index for each multiplication array. */
  int idx_A;            /* arrayA */
  int idx_B;            /* arrayB */
  int idx_C;            /* arrayC */


  /* Initialize "arrayA", "max_score", and "min_score" for
   * the first matrix column.*/
  max_score = STmax_score[0];
  min_score = STmin_score[0];
  min_idx_A = -(STmax_int_range - max_score);
  if (min_score > min_idx_A) min_idx_A = min_score - 1;
  arrayA = -min_idx_A + STarrayA;
  for (i = 0; i < T_SIZE; ++i)
    {
      if (weight_mat[i][0] <= -INFINITY) idx_A = min_idx_A;
      else
	{
	  idx_A = ROUND_TO_INT(Multiple * weight_mat[i][0]);
	  if (idx_A < min_idx_A) idx_A = min_idx_A;
	}

      if (arrayA[idx_A] > -INFINITY)
	arrayA[idx_A] = SUM_LN(Ln_P1[i], arrayA[idx_A]);
      else arrayA[idx_A] = Ln_P1[i];
    }

  for (j = 1; j < width; ++j)
    {
      /* Insert the probabilities for the current column into "arrayB[]".
       * Set the "STarray_index[]" to indicate once which scores occur
       * in the current column. */
      min_idx_B = STmax_score[j] - STmax_int_range;
      arrayB = -min_idx_B + STarrayB;
      for (i = 0, diff_sub = 0; i < T_SIZE; ++i)
	{
	  if ((weight_mat[i][j] <= -INFINITY) ||
	     ((idx_B = ROUND_TO_INT(Multiple * weight_mat[i][j])) < min_idx_B))
	    idx_B = min_idx_B;

	  if (arrayB[idx_B] == 0.0)
	    {
	      arrayB[idx_B] = P[i][j];
	      STarray_index[diff_sub++] = idx_B;
	    }
	  else arrayB[idx_B] += P[i][j];
	}
      /* Convert the probabilities in "arrayB[]" to logarithms. */
      for (i = 0; i < diff_sub; ++i)
	arrayB[STarray_index[i]] = log(arrayB[STarray_index[i]]);


      /* Multiply "arrayA[]" by "arrayB[]". */
      new_max_score = max_score + STmax_score[j];
      new_min_score = min_score + STmin_score[j];
      min_idx_C = -(STmax_int_range - new_max_score);
      if (new_min_score > min_idx_C) min_idx_C = new_min_score - 1;
      arrayC = -min_idx_C + STarrayC;
      if (arrayA[min_idx_A] > -INFINITY)
	{
	  arrayC[min_idx_C] = arrayA[min_idx_A] + LN_PROB_SUM;
	  arrayA[min_idx_A] = -INFINITY;
	}
      for (idx_A = max_score; idx_A > min_idx_A; --idx_A)
	{
	  if (arrayA[idx_A] > -INFINITY)
	    {
	      for (i = 0; i < diff_sub; ++i)
		{
		  idx_B = STarray_index[i];

		  if ((idx_C = idx_A + idx_B) < min_idx_C) idx_C = min_idx_C;

		  ln_prob_AB = arrayA[idx_A] + arrayB[idx_B];
		  ln_prob_C = arrayC[idx_C];

		  if (ln_prob_C <= -INFINITY) arrayC[idx_C] = ln_prob_AB;
		  else arrayC[idx_C] = SUM_LN(ln_prob_AB, ln_prob_C);
		}

	      /* Zero the current element of "arrayA[]". */
	      arrayA[idx_A] = -INFINITY;
	    }
	}


      /* Zero "arrayB[]". */
      for (i = 0; i < diff_sub; ++i) arrayB[STarray_index[i]] = 0.0;

      /* Exchange pointers to "STarrayA[]" and "STarrayC[]". */
      temp_array = STarrayA;
      STarrayA = STarrayC;
      STarrayC = temp_array;
      arrayA = arrayC;
      min_idx_A = min_idx_C;

      /* Update "max_score". */
      max_score = new_max_score;
      min_score = new_min_score;
    }

  /* Transfer information on arrayA into external variables.
   * Max_int_score: the maximum index determined in determine_variables().
   * Min_int_score: the minimum index might change from value set
   *                in determine_variables() due to rounding error.
   *      ArrayA[]: holds final ln(probabilities) of integer scores. */
  Min_int_score = min_idx_A;    
  ArrayA = arrayA;
}


/* Determine the logarithm of the p-values and place in Marginal_prob[]. */
void determine_marginal_prob()
{
  int i, j, k;
  int max_score;    /* The greater of the maximum integral score determined
		     * from the maximum exact score and the integral matrix. */
  double delta_ln_p;/* The change in the ln(p-value) between
		     * two observed integral scores. */


  /* Allocate space for "Marginal_prob[]". */
  max_score = ROUND_TO_INT(Multiple * Max_exact_score);
  if (Max_int_score > max_score) max_score = Max_int_score;
  Marginal_prob = (double *)ckrealloc(Marginal_prob, (max_score-Min_int_score+1)*sizeof(double));
  Marginal_prob -= Min_int_score;


  /* Determine the ln(p-value) of the Max_int_score. */
  Marginal_prob[Max_int_score] = ArrayA[Max_int_score];


  /* Copy ln(p-values) for scores greater than Max_int_score. */
  for (i = Max_int_score + 1; i <= max_score; ++i)
    Marginal_prob[i] = Marginal_prob[Max_int_score];


  /* Determine the main part of Marginal_prob[]. */
  for (i = Max_int_score - 1, j = Max_int_score;
       i > Min_int_score; --i)
    {
      if (ArrayA[i] > -INFINITY)
	{
	  Marginal_prob[i] = SUM_LN(ArrayA[i], Marginal_prob[j]);

	  /* Interpolate ln(p-values) for unobserved scores. */
	  delta_ln_p = (Marginal_prob[i] - Marginal_prob[j]) / (double)(j - i);
	  for (k = j--; j > i; k = j--)
	    Marginal_prob[j] = Marginal_prob[k] + delta_ln_p;
	}
    }

  /* Copy ln(p-values) for scores less than the lowest observed score j. */
  for (k = j - 1; k > Min_int_score; --k) Marginal_prob[k] = Marginal_prob[j];

  /* Determine the ln(p-value) for the Min_int_score. */
  Marginal_prob[i] = SUM_LN(ArrayA[i], Marginal_prob[j]);


#ifdef PRINTRANGE
  /* Print minimum score for numerically estimating p-values. */
  printf("           minimum score for calculating p-values: %8.3f\n",
	 (double)(Min_int_score + 1) / Multiple);

  /* Print the maximum and minimum numerically calculated p-values. */
  printf("       maximum ln(numerically calculated p-value): %8.3f\n",
	 Marginal_prob[Min_int_score + 1]);
  printf("       minimum ln(numerically calculated p-value): %8.3f\n",
	 Marginal_prob[Max_int_score]);

  printf("\n");
#endif


#ifdef PRINTDISTRIBUTION
  printf("a priori letter probabilities:\n");
  for (i = 0; i < T_SIZE; ++i) printf("%5.3f\n", P[i]);
  printf("\n");

  printf("Multiple: %.15g\n", Multiple);
  printf("\n");

  printf("integral score | probability | p-value | ln(p-value)\n");
  for (i = max_score; i > Max_int_score; i--)
    printf("%11d:          0 %10g %11f\n",
	   i, exp(Marginal_prob[i]), Marginal_prob[i]);
  for (i = Max_int_score; i >= Min_int_score; i--)
    printf("%11d: %10g %10g %11f\n",
	   i, exp(ArrayA[i]), exp(Marginal_prob[i]), Marginal_prob[i]);
#endif
}

/* Determine and print the cutoff score given the logarithm
 * of the target p-value. */
void cutoff_approx(ln_prob, note)
     double ln_prob;    /* The target p-value. */
     char *note;        /* Note describing the source of the cutoff p-value. */
{
  int i, j, k;
  double ln_prob_score;    /* ln(probability * score) */
  double avg;              /* Average score above the cutoff
			    * [initially holds ln(average)]. */


  /* Determine the cutoff score and the average score. */
  for (i = Max_int_score, j = Max_int_score + 1,
       k = Max_int_score - Min_int_score, avg = -INFINITY;
       (k > 0) && (Marginal_prob[i] <= ln_prob); --i, --k)
    {
      if (ArrayA[i] > -INFINITY)
	{
	  j = i;
	  ln_prob_score = ArrayA[i] + log((double)k);
	  avg = SUM_LN(avg, ln_prob_score);
	}
    }
#ifndef HIGH_CUTOFF
  if (Marginal_prob[j] < ln_prob)
    {
      for ( ; (k > 0) && (ArrayA[i] <= -INFINITY); --i, --k);

      if (k > 0)
	{
	  j = i;
	  ln_prob_score = ArrayA[i] + log((double)k);
	  avg = SUM_LN(avg, ln_prob_score);
	}
    }
#endif
  if (j <= Max_int_score)
    avg = (exp(avg - Marginal_prob[j]) + (double)Min_int_score) / Multiple;
  else avg = 0.0;


  /* Print the target cutoff p-value. */
  printf("ln(cutoff p-value)%s: %8.3f\n", note, ln_prob);

  /* Print the cutoff score and the achieved p-value. */
  if (j > Max_int_score)
    {
      printf("      maximum numerically calculated cutoff score: %8.3f\n",
	     (double)Max_int_score / Multiple);

      printf("       minimum ln(numerically calculated p-value): %8.3f\n",
	     Marginal_prob[Max_int_score]);
    }
  else if ((i == Min_int_score) &&
	   (ArrayA[Min_int_score] > -INFINITY) &&
	   (Marginal_prob[j] < ln_prob))
    {
      /* Information for the minimum score for numerically
       * estimating p-values is now printed automatically above. */
    }
  else
    {
      printf("              numerically calculated cutoff score: %8.3f\n",
	     (double)j / Multiple);

      printf("        ln(numerically calculated cutoff p-value): %8.3f\n",
	     Marginal_prob[j]);
    }

  /* Print the average score above the cutoff. */
  printf("average score above numerically calculated cutoff: %8.3f\n", avg);

  printf("\n");
}


/* Determine and print the cutoff score given the average score of
 * the L-mers above the cutoff. */
double cutoff_approx_avg(avg_score)
     double avg_score;       /* The target average score. */
{
  int i, j, k;
  double ln_int_avg_score;/* ln(target average score converted to "integer") */
  double ln_prob_score;       /* ln(probability * score) */
  double ln_prob_score_sum;   /* ln(SUM ln_prob_score). */
  double ln_int_avg;          /* ln(average score above the current cutoff). */
  double avg;                 /* Average score above the cutoff. */
  double cutoff_score;        /* The cutoff score. */


  /* Convert the targeted average score to its "integral" representation
   * and take its logarithm. */
  if (avg_score * Multiple < (double)(Min_int_score + 1))
    {
      printf("minimum numerically calculated cutoff score: ");
      printf("%8.3f\n", (double)(Min_int_score + 1) / Multiple);
      printf("    targeted average score above the cutoff: %8.3f\n",
	     avg_score);
      printf("\n");

      return(-INFINITY);
    }
  ln_int_avg_score = log(avg_score * Multiple - (double)Min_int_score);

  /* Determine the integral cutoff score and the average score. */
  for (i = Max_int_score, j = Max_int_score + 1,
       k = Max_int_score - Min_int_score,
       ln_int_avg = INFINITY, ln_prob_score_sum = -INFINITY;
       (k > 0) && (ln_int_avg > ln_int_avg_score); --i, --k)
    {
      if (ArrayA[i] > -INFINITY)
	{
	  j = i;
	  ln_prob_score = ArrayA[i] + log((double)k);
	  ln_prob_score_sum = SUM_LN(ln_prob_score_sum, ln_prob_score);

	  ln_int_avg = ln_prob_score_sum - Marginal_prob[j];
	}
    }
  avg = (exp(ln_int_avg) + (double)Min_int_score) / Multiple;

  if ((k == 0) && (ArrayA[Min_int_score] > -INFINITY) &&
      (ln_int_avg > ln_int_avg_score))
    {
      printf("THE MINIMUM SCORE FOR CALCULATING PROBABILITIES\n");
      printf("IS TOO HIGH TO FIND THE TARGETED AVERAGE SCORE.\n");
    }

  /* Determine the cutoff score. */
  cutoff_score = (double)j / Multiple;

  /* Print the target average score, the cutoff score, the p-value,
   * and the calculated average score. */
  printf("  targeted average score above the cutoff: %8.3f\n", avg_score);
  printf("      numerically calculated cutoff score: %8.3f\n", cutoff_score);
  printf("                   ln(calculated p-value): %8.3f\n",
	 Marginal_prob[j]);
  printf("calculated average score above the cutoff: %8.3f\n", avg);
  printf("\n");

  return(cutoff_score);
}

/* Extract the p-value from Marginal_prob[] for the given score.
 * Print the p-value. */
void extract_marginal_prob(score)
     double score;
{
  int int_score = ROUND_TO_INT(Multiple * score);  /* Integral score. */
  double ln_prob;                                  /* ln(p-value) */


  if (Max_int_range > 0.0)
    {
      if ((score >= -INFINITY / 2.0) &&
	  (score >= Min_score - ERROR * fabs(Min_score)) &&
	  (score >= Min_exact_score - ERROR * fabs(Min_exact_score)))
	{
	  if (int_score > Min_int_score) ln_prob = Marginal_prob[int_score];
	  else ln_prob = Marginal_prob[Min_int_score + 1];

	  printf("  ln(p-value)= %7.2f", ln_prob);
	}

      else if (Min_score <= Min_exact_score)
	{
	  ln_prob = 0.0;
	  printf(" ln(p-value)= %5.2f", ln_prob);
	}

#if 0
      else printf("  ln(p-value)> %7.2f", Marginal_prob[Min_int_score + 1]);
#endif
    }
}
