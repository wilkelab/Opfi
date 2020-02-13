#ifndef MYASSERT_H
#define MYASSERT_H

#ifdef assert
#undef assert
#endif

void assert_(const char *, const char *, unsigned);
#define assert(exp) (void)( (exp) || (assert_(#exp, __FILE__, __LINE__), 0) )

#endif // MYASSERT_H
