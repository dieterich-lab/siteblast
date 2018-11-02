#ifndef NIB_H
#define NIB_H
#include "int_types.h"
#include "util.h"

unsigned char *seq_freadnib(FILE *, i32_t , i32_t , i32_t *);
unsigned char *seq_readnib(const char *, i32_t , i32_t , i32_t *);
void seq_fwritenib(FILE *, unsigned const char *, u32_t );
void seq_writenib(char *, unsigned const char *, u32_t );
#endif
