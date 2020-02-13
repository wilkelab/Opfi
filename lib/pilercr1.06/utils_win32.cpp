#include "pilercr.h"

#if	WIN32
#include <windows.h>
#include <psapi.h>

double GetRAMSize()
	{
	MEMORYSTATUS MS;
	GlobalMemoryStatus(&MS);
	double m = (double) MS.dwTotalPhys;
	if (m > 1.8e9)
		m = 1.8e9;
	return (unsigned) m;
	}

static unsigned g_uPeakMemUseBytes;

unsigned GetPeakMemUseBytes()
	{
	return g_uPeakMemUseBytes;
	}

unsigned GetMemUseBytes()
	{
	HANDLE hProc = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS PMC;
	BOOL bOk = GetProcessMemoryInfo(hProc, &PMC, sizeof(PMC));
	if (!bOk)
		return 1000000;
	unsigned uBytes = (unsigned) PMC.WorkingSetSize;
	if (uBytes > g_uPeakMemUseBytes)
		g_uPeakMemUseBytes = uBytes;
	return uBytes;
	}

#endif
