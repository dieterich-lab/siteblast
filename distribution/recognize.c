static char rcsid[] = "$Id: recognize.c,v 4.0 2004/08/26 09:49:25 mmichael Exp mmichael $";

#include <string.h>
#include <stdlib.h>
#include "bz_all.h"

void add_ident(Recognizer r, char *pattern, int pattern_id,
	       struct iupac_information *extra_information, char extra_type);
void stop_adding(Recognizer r);

typedef struct goto_node Goto_Node;
typedef struct move_node Move_Node;

typedef struct name_node {
  struct name_node *next; /* points to the next name on the output list */
  int pattern_length;
  int pattern_id;
  struct iupac_information *extra_information;
  char extra_type;
} Name_Node;

struct move_node {
  Move_Node *next;      /* points to the next node on the move list */
  Goto_Node *state;     /* the next state for this character */
  unsigned char c;
};

struct goto_node {
  Name_Node *output;    /* list of words ending in this state */
  Move_Node *moves;     /* list of possible moves */
  Goto_Node *fail;      /* and where to go when no move fits */
  Goto_Node *next;      /* next goto node with same depth */
};

struct recognizer {
  Goto_Node *root[256]; /* might want 128, depending on the character set */
  int max_depth;
  Goto_Node **depths; /* an array, max_depth long, of lists of goto_nodes,
                         created while adding ids, used while building
                         the failure functions */
};

static Goto_Node *goto_lookup(unsigned char c, Goto_Node *g)
{
  Move_Node *m = g->moves;
  while (m && m->c != c)
    m = m->next;
  return m ? m->state : NULL;
}

Recognizer new_recognizer()
{
  Recognizer r = (Recognizer) ckcalloc(1, sizeof(struct recognizer));
  r->max_depth = 10;
  r->depths = (Goto_Node **) ckcalloc(r->max_depth, sizeof(Goto_Node *));
  return r;
}

void add_ident(Recognizer r, char *pattern, int pattern_id,
	       struct iupac_information *extra_information, char extra_type){
  int depth = 2;
  char *p = pattern;
  unsigned char c = *p++;
  Goto_Node *q = r->root[c];
  Name_Node *n;
  char *copy;

  if (!q) 
  {
    q = (Goto_Node *) ckcalloc(1, sizeof(Goto_Node));
    r->root[c] = q;
    q->next = r->depths[1];
    r->depths[1] = q;
  }

  c = *p++;
  while (c) {
    Goto_Node *new = goto_lookup(c, q);
    if (!new)
    {
      Move_Node *new_move = (Move_Node *) ckalloc(sizeof(Move_Node));
      new = (Goto_Node *) ckcalloc(1, sizeof(Goto_Node));
      new_move->state = new;
      new_move->c = c;
      new_move->next = q->moves;
      q->moves = new_move;
      if (depth == r->max_depth)
      {
	int i;
	Goto_Node **new_depths = (Goto_Node **) ckcalloc(2*depth, sizeof(Goto_Node *));
	r->max_depth = 2 * depth;
	for (i=0; i<depth; i++)
	  new_depths[i] = r->depths[i];
	free(r->depths);
	r->depths = new_depths;
      }

      new->next = r->depths[depth];
      r->depths[depth] = new;
    }

    q = new;
    depth++;
    c = *p++;
  }

  if(!q->output){
    n=q->output= (Name_Node *) ckalloc(sizeof(Name_Node));
  }
  else{
    n=q->output;
    while(n->next) n=n->next;
    n->next = (Name_Node *) ckalloc(sizeof(Name_Node));
    n=n->next;
  }
  n->next = NULL;
  n->pattern_length = strlen(pattern);
  n->pattern_id=pattern_id;
  n->extra_information=extra_information;
  n->extra_type=extra_type;
}

void stop_adding(Recognizer r)
{
  int depth;
  for (depth=1; depth<r->max_depth; depth++) {
    Goto_Node *g = r->depths[depth];
    while (g) {
      Move_Node *m = g->moves;
      while (m) {
        unsigned char a = m->c;
        Goto_Node *s = m->state;
        Goto_Node *state = g->fail;
        while (state && !goto_lookup(a, state))
          state = state->fail;
        if (state)
          s->fail = goto_lookup(a, state);
        else
          s->fail = r->root[a];
        if (s->fail) {
          Name_Node *p = s->fail->output;
          while (p) {
            Name_Node *q = (Name_Node *) ckalloc(sizeof(Name_Node));
	    /* depending on memory deallocation 
	       strategy, we may need to copy this */
	    q->pattern_length = p->pattern_length;
	    q->pattern_id = p->pattern_id;
            q->extra_information = p->extra_information;
	    q->extra_type = p->extra_type;
            q->next = s->output;
            s->output = q;
            p = p->next;
          }
        }
        m = m->next;
      }
      g = g->next;
    }
  }
}

void search_for_ident_ext(Recognizer r, char *input, CallbackExt f, 
			  void *closure, 
			  /*param for f*/void *param){
  Goto_Node *state = NULL;
  char *current = input;
  unsigned char c = (unsigned char) *current++;
#ifdef TRACK_SEARCH
  int len=strlen(input);
#endif
  while (c) {
#ifdef TRACK_SEARCH
    /*dprint*/
    if(!((current-input)%(len/50))){
      printf("#");
      fflush(stdout);
    }
#endif

    {
      while (state && !goto_lookup(c, state))
	state = state->fail;
      state = state ? goto_lookup(c, state) : r->root[c];
    }
    {
      if (state) {
	Name_Node *p = state->output;
	while (p) {
	  /*dprint
	  printf("Found \"%*.*s\".\t\t(%s:%d)\n", p->pattern_length, p->pattern_length, 
		 current - p->pattern_length,__FILE__,__LINE__);*/
	  f(closure, p->pattern_length, current - p->pattern_length, 
	    clone_iupac_extra_info(p->extra_information), p->extra_type, p->pattern_id, param);
	  p = p->next;
	}
      }
    }
    
    c = *current++;
  }
#ifdef TRACK_SEARCH
  /*for(c=0;c<55;c++) printf("\b");*/
  printf("\n");
#endif
}
