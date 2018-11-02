#ifndef PIPLIB_INT_TYPES_H
#define PIPLIB_INT_TYPES_H
/* $Id: int_types.h,v 4.0 2004/08/26 09:49:25 mmichael Exp $ */

// Standard C (ISO/IEC 9899:1999 Programming language - C)
//   http://std.dkuug.dk/JTC1/SC22/WG14/www/standards
// specifies stdint.h and inttypes.h so you can ask for
// exactly the right size integers.  Unfortunately, many
// unix vendors don't meet this standard yet.  So we have
// to configure this by hand.

// At the moment, we just need types that are big enough
// rather than exactly enough.

#define I32MAX 2147483647

// Try to guess. 
#include <limits.h>

#if INT_MAX < I32MAX
#error "INT_MAX is too small"
#endif

#if LONG_MAX == I32MAX // ILP32, like x86
typedef int i32_t;
typedef unsigned int u32_t;
typedef long long i64_t;
typedef unsigned long long u64_t;
#elif LONG_MAX > I32MAX // I32 LP64, like sparc v9
typedef int i32_t;
typedef unsigned int u32_t;
typedef long i64_t;
typedef unsigned long u64_t;
#elif INT_MAX > I32MAX  // ILP64, alpha?
typedef int i32_t;
typedef unsigned int u32_t;
typedef long i64_t;
typedef unsigned long u64_t;
#else
#error "cannot guess int sizes"
#endif

#endif
