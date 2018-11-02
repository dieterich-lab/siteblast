/* bz_inner.c -- given a list of "outer" alignments, sorted by beginning
   position in seq1, try to interpolate higher-sensitivity "inner alignments"
   between properly-ordered pairs of outer alignments or at ends of such
   chains.  Treat alignments in order, keeping a set of "active" alignments,
   i.e., which have already been treated but whose end position in seq1 is
   within WINDOW bp of the current begin position.  (We expect this set to be
   quite small on average, so a linked list is quite adequate.)  For each
   alignment, A, do the following:

   1. Retire all active alignments that end more than WINDOW bp before the
      current alignment begins.  If an alignment being retired is the end of
      a chain of (1 or more) alignments, look for weak matches to its right.
   2. Look for a current alignment, B, that overlap A, i.e., B's end point is
      within WINDOW diagonals of where A starts, and follows A in one of the
      sequences.  In the abnormal case that B ends after A ends (relative to
      one of the sequences, mark A as "not the right end of a chain".
   3. If no alignments overlap A, look for the alignment B that ends before
      A (in both sequences) and is closest.  If it comes within WINDOW bp in
      both sequences, search for weak alignments between A and B.  If no such
      B exist, think of A as being on the left end of a chain of (1 or more)
      outer alignments, and search for weak alignments to its left.
   When all alignment have been treated that way, retire the remaining active
   alignments.
*/

#include "bz_all.h"

#define CAUTIOUS
//#define DEBUG_INNER
//#define STATS

enum {
  WINDOW  = 20000,	/* max bp between outer alignments to search */
  NEW_W = 7		/* word size for inner alignments */
};

static align_t *inner_list;

// set of alignments that have been checked, but may lie within WINDOW bp
// of the start of a future alignment, relative to sequence #1
typedef struct _active {
	align_t *align;
	int right_end;	// right end of seen-so-far chain
	struct _active *next;
} Active;

static Active *active;

static SEQ *sf1;
static SEQ *sf2;
static gap_scores_t gs;
static ss_t ss;
static ss_t sss;
static int Y;
static int innerK;
static TBack tback;
static bz_flags_t bz_flags;
static connect_t Connect;
static int MSP_Scale;

static void focus(int b1, int e1, int b2, int e2);

// merge two beg1-ordered lists of alignments into a beg1-ordered list
static align_t *merge_align(align_t *a, align_t *b) {
	align_t *ret, *tail;

#ifdef CAUTIOUS
	for (tail = a; tail; tail = tail->next_align)
		if (tail->next_align && tail->next_align->beg1 < tail->beg1)
			fatalf("merge_align: first list out of order at %d",
			  tail->next_align->beg1);
	for (tail = b; tail; tail = tail->next_align)
		if (tail->next_align && tail->next_align->beg1 < tail->beg1)
			fatalf("merge_align: second list out of order at %d",
			  tail->next_align->beg1);
#endif

	if (b == NULL)
		return a;
	if (a == NULL)
		return b;
	if (a->beg1 <= b->beg1) {
		ret = tail = a;
		a = a->next_align;
	} else {
		ret = tail = b;
		b = b->next_align;
	}

	while (a != NULL && b != NULL)
		if (a->beg1 <= b->beg1) {
			tail = tail->next_align = a;
			a = a->next_align;
		} else {
			tail = tail->next_align = b;
			b = b->next_align;
		}
	if (a == NULL)
		tail->next_align = b;
	else
		tail->next_align = a;
	return ret;
}

// put an alignment at the front of the active list
static void activate(align_t *a, int status) {
	Active *c = ckalloc(sizeof(*c));

	c->align = a;
	c->right_end = status;
	c->next = active;
	active = c;
}

// an alignment is being de-activated -- it could it be the outer alignment
// at the end of a chain, with inner alignments to its right
static void retire(Active *c) {
	int a1, a2, b1, b2;

	if (c->right_end) {
		b1 = c->align->end1;
		b2 = c->align->end2;
		a1 = MIN(b1+WINDOW/2, SEQ_LEN(sf1));
		a2 = MIN(b2+WINDOW/2, SEQ_LEN(sf2));
		focus(b1, a1, b2, a2);
	}
}

// create a bogus SEQ item
static SEQ *fakeSEQ(SEQ *sf, int b, int e) {
	SEQ *xf = ckalloc(sizeof(SEQ));
	uchar *s = SEQ_CHARS(sf);
	int i, L;

	xf->fp = NULL;
	xf->flags = xf->count = 0;
	xf->offset = 0;
	xf->maskname = xf->fname = xf->header = NULL;
	xf->from = 1;
	xf->hlen = 0;
	xf->slen = L = e - b + 1;
	xf->seq = ckalloc((L+1)*sizeof(uchar));
	for (i = 0; i < L; ++i)
		xf->seq[i] = s[b+i-1];
	xf->seq[L] = '\0';
	return xf;
}

// perform a high-sensitivity alignment in a specified rectangle
static void focus(int b1, int e1, int b2, int e2) {
	blast_table_t *bt;
	static msp_table_t *mt = 0;
	msp_t *p;
	align_t *a = NULL, *aa;
	SEQ *xf1, *xf2;
	uchar *rf1;
	int f = 0;

#ifdef DEBUG_INNER
	fprintf(stderr, "  focus: (%d,%d) to (%d,%d)\n", b1, b2, e1, e2);
#endif
	/* create SEQ objects for the tiny sequences */
	xf1 = fakeSEQ(sf1, b1, e1);
	xf2 = fakeSEQ(sf2, b2, e2);

	bt = blast_table_new(xf1, NEW_W);
	mt = msp_new_table();
	mt->num = 0;
	mt = blast_search(xf1, xf2, bt, mt, sss, Y, innerK, 0, 0); // P=0,T=0

	/* chain */
	(void)msp_make_chain(mt, bz_flags.G, bz_flags.G, MSP_Scale, Connect);
	for (p = MSP_TAB_FIRST(mt); MSP_TAB_MORE(mt,p); p = MSP_TAB_NEXT(p)) {
		p->filter = 1 - p->cum_score;
		f |= p->filter;
	}

	msp_compress(mt);
	msp_sort_pos1(mt);

	rf1 = reverse_seq(SEQ_CHARS(xf1), SEQ_LEN(xf1));
	if (mt && MSP_TAB_NUM(mt) != 0)
		a = bz_extend_msps(xf1, rf1, xf2, mt, &gs, ss, Y, innerK,
		  tback);
        blast_table_free(bt);
	msp_free_table(mt);
/*
	free(rf1);
	free(xf1->seq);
	free(xf1);
	free(xf2->seq);
	free(xf2);
*/

	/* shift the entries of each alignment */
	for (aa = a; aa != NULL; aa = aa->next_align) {
		aa->seq1 = SEQ_CHARS(sf1);
		aa->seq2 = SEQ_CHARS(sf2);
		aa->beg1 += (b1-1);
		aa->end1 += (b1-1);
		aa->beg2 += (b2-1);
		aa->end2 += (b2-1);
	}
#ifdef DEBUG_INNER
	for (aa = a; aa; aa = aa->next_align)
		fprintf(stderr, "    new alignment (%d,%d) to (%d,%d)\n",
		  aa->beg1, aa->beg2, aa->end1, aa->end2);
#endif

	inner_list = merge_align(a, inner_list);
}

// entry point: interpolate inner alignments in a chain of outer alignments
align_t *inner(align_t *a, SEQ *Sf1, SEQ *Sf2, gap_scores_t Gs, ss_t Ss,
  ss_t Sss, int y, int InnerK, TBack Tback, bz_flags_t Bz_flags,
  connect_t connect, int mSP_Scale) {
	align_t *A, *B;
	Active *c, *c2;
	int i, j, a1, a2, b1, b2, dist_B, d, left_end, overlap, diag_diff;

	sf1 = Sf1;
	sf2 = Sf2;
	gs = Gs;
	for (i = 0; i < NACHARS; ++i)
		for (j = 0; j < NACHARS; ++j) {
			ss[i][j] = Ss[i][j];
			sss[i][j] = Sss[i][j];
		}
	Y = y;
	innerK = InnerK;
	tback = Tback;
	bz_flags = Bz_flags;
	Connect = connect;
	MSP_Scale = mSP_Scale;

	active = NULL;
	inner_list = NULL;

#ifdef CAUTIOUS
	for (B = a; B; B = A)
		if ((A = B->next_align) != NULL && A->beg1 < B->beg1)
			fatal("outer alignments out of order");
#endif

#ifdef DEBUG_INNER
	fprintf(stderr, "outer alignments\n");
	for (A = a; A; A = A->next_align)
		fprintf(stderr, "  (%d,%d) to (%d,%d)\n",
		  A->beg1, A->beg2, A->end1, A->end2);
#endif
	if (a == NULL)
		return NULL;
	for (A = a; A; A = A->next_align) {
		a1 = A->beg1;
		a2 = A->beg2;
#ifdef DEBUG_INNER
		fprintf(stderr, "treating alignment (%d,%d) to (%d,%d)\n",
		  a1, a2, A->end1, A->end2);
#endif

		// Delete any expired alignment from the active list.
		while ((c = active) != NULL && a1 - c->align->end1 > WINDOW) {
			retire(c);
			active = c->next;
			// free(c);
		}
		c = active;
		while (c != NULL && (c2 = c->next) != NULL)
			if (a1 - c2->align->end1 > WINDOW) {
				retire(c2);
				c->next = c2->next;
				// free(c2);
			} else
				c = c2;
#ifdef CAUTIOUS
		for (c = active; c; c = c->next) { 
			b1 = c->align->end1;
			if ( !(a1 <= b1 + WINDOW) )
				fatalf("inner: impossible");
		}
#endif

		// Look for an active alignment that overlaps A.
		overlap = 0;
		for (c = active; c != NULL; c = c2) {
			c2 = c->next;
			B = c->align;
			b1 = B->end1;
			b2 = B->end2;
			diag_diff = (b2 - b1) - (a2 - a1);
			if (diag_diff < 0)
				diag_diff = -diag_diff;
#ifdef DEBUG_INNER
			fprintf(stderr,"  overlap with (%d,%d) to (%d,%d)?",
			  B->beg1, B->beg2, b1, b2);
#endif
			if (diag_diff <= WINDOW && (b1 >= a1 || b2 >= a2)) {
				overlap = 1;
#ifdef DEBUG_INNER
				fprintf(stderr, "  yes\n");
#endif
				if (b1 < A->end1 && b2 < A->end2)
					// B ends properly -- before A ends
					c->right_end = 0;
				else
					break;
			}
#ifdef DEBUG_INNER
			else
				fprintf(stderr, "  no\n");
#endif
		}
#ifdef DEBUG_INNER
		fprintf(stderr, "    overlap = %d\n", overlap);
#endif
		if (overlap) {
			// if c==NULL, we didn't break out of the above loop,
			// i.e., any overlaps are proper, so c->right_end == 1
			activate(A, (c == NULL));
			continue; // don't try to generate inner alignment
				  // on the left side of A
		}

		// Let B be an active alignment that ends at least 0 bp and at
		// most WINDOW bp before A starts (relative to both sequences).
		// Pick B to be closest to A among the candidates.
		B = NULL;
		dist_B = 3*WINDOW;
		left_end = 1;	// A is the first outer alignment in a chain
		for (c = active; c; c = c->next) { 
			b1 = c->align->end1;
			b2 = c->align->end2;
			if (b1 < a1 && b2 < a2 && a2 < b2 + WINDOW) {
				left_end = 0;
				if (c->right_end &&
				    (d = a1 - b1 + a2 - b2) < dist_B) {
					B = c->align;
					dist_B = d;
				}
				c->right_end = 0;
			}
		}
		if (B != NULL) {
			b1 = B->end1;
			b2 = B->end2;
			focus(b1, a1, b2, a2);
		} else if (left_end) {
			// A could be the first outer alignment in a chain,
			// with inner alignments to its left
			b1 = MAX(1, a1-WINDOW/2);
			b2 = MAX(1, a2-WINDOW/2);
			focus(b1, a1, b2, a2);
		}

		activate(A, 1);
	}

	// align in windows after each chain-ending active alignment
	for (c = active; c; c = c->next)
		retire(c);

#ifdef STATS
	a1 = 0;
	for (A = inner_list; A; A = A->next_align)
		a1 += (A->end1 - A->beg1 + 1);
	fprintf(stderr, "  inner alignments added %d bp, ", a1);
#endif

	a = merge_align(a, inner_list);

#ifdef DEBUG_INNER
	fprintf(stderr, "final alignments\n");
	for (A = a; A; A = B) {
		fprintf(stderr, "  (%d,%d) to (%d,%d)\n",
		    A->beg1, A->beg2, A->end1, A->end2);
		if ((B = A->next_align) == NULL)
			break;
		if (A->beg1 > B->beg1)
			fatalf("final intervals out of order");
		if ( MIN(A->end1,B->end1) >= MAX(A->beg1, B->beg1) &&
		     MIN(A->end2,B->end2) >= MAX(A->beg2, B->beg2) )
			fprintf(stderr, "OVERLAP\n");
	}
	exit(0);
#endif

#ifdef STATS
	a1 = 0;
	for (A = a; A; A = A->next_align)
		a1 += (A->end1 - A->beg1 + 1);
	fprintf(stderr, "bringing total to %d bp\n", a1);
#endif

		return a;
}
