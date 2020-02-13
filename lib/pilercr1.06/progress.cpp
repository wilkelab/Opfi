#include "pilercr.h"
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

static time_t g_StartTime = time(0);
static time_t g_StartTimeStep;
static char *g_Msg = strsave("");
static bool g_NL = false;
bool g_ShowProgress = true;

int GetElapsedSecs()
	{
	return (int) (time(0) - g_StartTime);
	}

void ShowProgress(bool Show)
	{
	g_ShowProgress = Show;
	}

static char *Resources()
	{
	static char Str[64];

	int ElapsedSecs = GetElapsedSecs();
	unsigned MemUseMB = (GetMemUseBytes() + 500000)/1000000;
	unsigned RAMMB = (unsigned) ((GetRAMSize() + 500000)/1000000);
	unsigned MemUsePct = (MemUseMB*100)/RAMMB;

	sprintf(Str, "%4d s  %6d Mb (%3d%%) ", ElapsedSecs, MemUseMB, MemUsePct);
	return Str;
	}

void ProgressStep(int StepIndex, int StepCount)
	{
	if (!g_ShowProgress)
		return;

	double Pct = (double) StepIndex*100.0 / (double) StepCount;

	fprintf(stderr, "\r%s %6.2f%%  %s", Resources(), Pct, g_Msg);
//	fprintf(stderr, "%s  %6.2f%%  %s Step %d of %d\n", Resources(), Pct, g_Msg, StepIndex, StepCount);
	g_NL = false;
	}

void ProgressStart(const char *Format, ...)
	{
	g_StartTimeStep = time(0);
	char Str[4096];
	va_list ArgList;
	va_start(ArgList, Format);
	vsprintf(Str, Format, ArgList);
	if (g_Msg != 0)
		free(g_Msg);
	g_Msg = strsave(Str);

	if (!g_ShowProgress)
		return;

	if (!g_NL)
		fprintf(stderr, "\n");
	ProgressStep(0, 1);
	}

void ProgressDone()
	{
	if (0 == g_Msg)
		return;

	Log("%s %s\n", Resources(), g_Msg);

	if (g_ShowProgress)
		{
		fprintf(stderr, "\r%s %6.2f%%  %s\n", Resources(), 100.0, g_Msg);
		g_NL = true;
		}
	free(g_Msg);
	g_Msg = 0;
	}

void Progress(const char *Format, ...)
	{
	if (!g_ShowProgress)
		return;

	char Str[4096];
	va_list ArgList;
	va_start(ArgList, Format);
	vsprintf(Str, Format, ArgList);
	Log("%s %s\n", Resources(), Str);

	if (g_ShowProgress)
		{
		if (!g_NL)
			fprintf(stderr, "\n");
		fprintf(stderr, "%s %s\n", Resources(), Str);
		g_NL = true;
		}
	}

void ProgressNoRes(const char *Format, ...)
	{
	if (!g_ShowProgress)
		return;

	char Str[4096];
	va_list ArgList;
	va_start(ArgList, Format);
	vsprintf(Str, Format, ArgList);

	Log("%s\n", Str);

	if (g_ShowProgress)
		{
		if (!g_NL)
			fprintf(stderr, "\n");
		fprintf(stderr, "%s\n", Str);
		g_NL = true;
		}
	}

void ProgressExit()
	{
	time_t t = time(0);
	Log("%s", asctime(localtime(&t)));
	Progress("Finished, total elapsed time %ld secs.", GetElapsedSecs());
	}
