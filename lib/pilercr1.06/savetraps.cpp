#include "pilercr.h"

static int GetTargetPos(int QueryPos, int DiagIndex)
	{
	return QueryPos - DiagIndex;
	}

static void SaveTrap(FILE *f, const Trapezoid &Trap)
	{
	int QueryFrom = Trap.bot;
	int QueryTo = Trap.top;
	int TargetFrom = GetTargetPos(Trap.bot, Trap.rgt);
	int TargetTo = GetTargetPos(Trap.top, Trap.lft);
	WriteGFFRecord(f, QueryFrom, QueryTo, TargetFrom, TargetTo, false, "trap", 0, 0);
	}

void SaveTraps(const Trapezoid *Traps)
	{
	const char *FileName = ValueOpt("traps");
	if (FileName == 0)
		return;

	FILE *f = CreateStdioFile(FileName);

	for (const Trapezoid *Trap = Traps; Trap; Trap = Trap->next)
		SaveTrap(f, *Trap);

	fclose(f);
	}

static void LogTrap(const Trapezoid &Trap)
	{
	int TargetFrom = GetTargetPos(Trap.bot, Trap.rgt);
	int TargetTo = GetTargetPos(Trap.top, Trap.lft);

	Log("%10d  %10d  %10d  %10d  %10d  %10d\n",
	  Trap.bot, Trap.top, Trap.lft, Trap.rgt, TargetFrom, TargetTo);
	}

void LogTraps(const Trapezoid *Traps)
	{
	if (!FlagOpt("logtraps"))
		return;

	Log("\n");
	Log("     QFrom         QTo        Left       Right       TFrom         TTo\n");
	Log("==========  ==========  ==========  ==========  ==========  ==========\n");
	for (const Trapezoid *Trap = Traps; Trap; Trap = Trap->next)
		LogTrap(*Trap);
	Log("\n");
	}
