#include "pilercr.h"

#if	!WIN32 && !(linux || __linux__)

double GetRAMSize()
	{
	return 1e9;
	}

static unsigned g_uPeakMemUseBytes = (unsigned) 500e6;

unsigned GetPeakMemUseBytes()
	{
	return g_uPeakMemUseBytes;
	}

unsigned GetMemUseBytes()
	{
	return (unsigned) 500e6;
	}

#endif
