#include "pilercr.h"

#if	linux || __linux__
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <fcntl.h>

double GetRAMSize()
	{
	const double DEFAULT_RAM = 1e9;
	static double RAMMB = 0;
	if (RAMMB != 0)
		return RAMMB;

	int fd = open("/proc/meminfo", O_RDONLY);
	if (-1 == fd)
		return DEFAULT_RAM;

	char Buffer[128];
	int n = read(fd, Buffer, sizeof(Buffer) - 1);
	close(fd);
	fd = -1;

	if (n <= 0)
		return DEFAULT_RAM;

	Buffer[n] = 0;
	char *pMem = strstr(Buffer, "Mem: ");
	if (0 == pMem)
		return DEFAULT_RAM;
	int Bytes = atoi(pMem+4);
	double m = Bytes;
	if (m > 1.8e9)
		m = 1.8e9;
	return m;
	}

static unsigned g_uPeakMemUseBytes;

unsigned GetPeakMemUseBytes()
	{
	return g_uPeakMemUseBytes;
	}

unsigned GetMemUseBytes()
	{
	static char statm[64];
	static int PageSize;
	if (0 == statm[0])
		{
		PageSize = sysconf(_SC_PAGESIZE);
		pid_t pid = getpid();
		sprintf(statm, "/proc/%d/statm", (int) pid);
		}

	int fd = open(statm, O_RDONLY);
	if (-1 == fd)
		return 1000000;
	char Buffer[64];
	int n = read(fd, Buffer, sizeof(Buffer) - 1);
	close(fd);
	fd = -1;

	if (n <= 0)
		return 1000000;

	Buffer[n] = 0;
	int Pages = atoi(Buffer);

	unsigned uBytes = Pages*PageSize;
	if (uBytes > g_uPeakMemUseBytes)
		g_uPeakMemUseBytes = uBytes;
	return uBytes;
	}

#endif
