#ifndef __lint
static const char rcsid[] =
"$Id: bz_dna.c,v 1.1 2004/05/05 15:36:24 mmichael Exp $";
#endif

#include "bz_all.h"

bool is_DNA(uchar *s, int len)
{
	int ACGT, i;

	ACGT = 0;
	for (i = 0; i < len; ++i)
		if (strchr("ACGTNXacgtnx", s[i]))
			++ACGT;
	for (i = 0; i < len; ++i)
		if (!is_dchar(s[i]))
			fatalf("Illegal character '%c' in sequence file.\n",
			   s[i]);
	return ((long long)10*ACGT >= (long long)9*len); /* ACGT >= 90% of len */
}
