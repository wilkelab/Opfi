#include "pilercr.h"

void assert_(const char *Exp, const char *FileName, unsigned LineNr)
	{
	fprintf(stderr, "Assertion failed(%s:%d) %s\n", FileName, LineNr, Exp);
	Quit("Assertion failed(%s:%d) %s\n", FileName, LineNr, Exp);
	}
