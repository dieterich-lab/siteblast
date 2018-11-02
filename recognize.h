/*
Based on
Alfred V. Aho and Margaret J. Corasick.
Efficient String Matching: An Aid for Bibliographic Search
Communications of the ACM, Vol. 18., Number 6, Juni 1975, S. 333--340

by
Preston Briggs: found in CTAN web/noweb/src/recognice.*
*/

#ifndef RECOGNICE__H 
#define RECOGNICE__H 

struct iupac_information;

typedef struct recognizer *Recognizer;
typedef void (*CallbackExt) (void *closure, int pattern_length, char *instance, 
			     void *extra_information, char extra_type,
			     int pattern_id, void *paramstruct);
Recognizer new_recognizer();
void add_ident(Recognizer r, char *pattern, int pattern_id,
	       struct iupac_information *extra_information, char extra_type);
void stop_adding(Recognizer r);
void search_for_ident_ext(Recognizer r, char *input, CallbackExt f, 
			  void *closure, void *paramstruct);
#endif


