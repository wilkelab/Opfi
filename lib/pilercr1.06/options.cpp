#include "pilercr.h"

struct VALUE_OPT
	{
	const char *m_pstrName;
	const char *m_pstrValue;
	};

struct FLAG_OPT
	{
	const char *m_pstrName;
	bool m_bSet;
	};

static VALUE_OPT ValueOpts[] =
	{
	"in",					0,
	"minhitlength",			0,
	"minid",				0,
	"maxmem",				0,
	"log",					0,
	"loga",					0,
	"out",					0,
	"hits",					0,
	"discardedhits",		0,
	"piles",				0,
	"filterout",			0,
	"diameter",				0,
	"traps",				0,
	"fhits",				0,
	"seq",					0,
	"cmp",					0,
	"ref",					0,
	"minrepeat",			0,
	"maxrepeat",			0,
	"minspacer",			0,
	"maxspacer",			0,
	"minrepeatratio",		0,
	"minspacerratio",		0,
	"minarray",				0,
	"mincons",				0,
	};
static int ValueOptCount = sizeof(ValueOpts)/sizeof(ValueOpts[0]);

static FLAG_OPT FlagOpts[] =
	{
	"quiet",				0,
	"version",				0,
	"help",					0,
	"options",				0,
	"loghits",				0,
	"logtraps",				0,
	"showhits",				0,
	"loghits",				0,
	"logimages",			0,
	"logalns",				0,
	"noinfo",				0,
	"segv",					0,
	"trimseqs",				0,
	};
static int FlagOptCount = sizeof(FlagOpts)/sizeof(FlagOpts[0]);

void CommandLineError(const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);

	char Str[4096];
	vsprintf(Str, Format, ArgList);
	Usage();
	fprintf(stderr, "\n*** Invalid command line ***\n");
	fprintf(stderr, "%s\n", Str);
	exit(1);
	}

static bool TestSetFlagOpt(const char *Arg)
	{
	for (int i = 0; i < FlagOptCount; ++i)
		if (!stricmp(Arg, FlagOpts[i].m_pstrName))
			{
			FlagOpts[i].m_bSet = true;
			return true;
			}
	return false;
	}

static bool TestSetValueOpt(const char *Arg, const char *Value)
	{
	for (int i = 0; i < ValueOptCount; ++i)
		if (!stricmp(Arg, ValueOpts[i].m_pstrName))
			{
			if (0 == Value)
				CommandLineError("Option -%s must have value\n", Arg);
			ValueOpts[i].m_pstrValue = strsave(Value);
			return true;
			}
	return false;
	}

bool FlagOpt(const char *Name)
	{
	for (int i = 0; i < FlagOptCount; ++i)
		if (!stricmp(Name, FlagOpts[i].m_pstrName))
			return FlagOpts[i].m_bSet;
	Quit("FlagOpt(%s) invalid", Name);
	return false;
	}

const char *ValueOpt(const char *Name)
	{
	for (int i = 0; i < ValueOptCount; ++i)
		if (!stricmp(Name, ValueOpts[i].m_pstrName))
			return ValueOpts[i].m_pstrValue;
	Quit("ValueOpt(%s) invalid", Name);
	return 0;
	}

const char *RequiredValueOpt(const char *Name)
	{
	const char *s = ValueOpt(Name);
	if (0 == s)
		CommandLineError("Required option -%s not specified\n", Name);
	return s;
	}

void ProcessArgVect(int argc, char *argv[])
	{
	if (argc == 0)
		{
		Usage();
		exit(0);
		}

	if (argc == 1)
		{
		if (strcmp("-help", argv[0]) == 0 || strcmp("--help", argv[0]) == 0)
			{
			Usage();
			exit(0);
			}
		}

	for (int iArgIndex = 0; iArgIndex < argc; )
		{
		const char *Arg = argv[iArgIndex];
		if (Arg[0] != '-')
			Quit("Command-line option \"%s\" must start with '-'\n", Arg);
		const char *ArgName = Arg + 1;
		if (TestSetFlagOpt(ArgName))
			{
			++iArgIndex;
			continue;
			}
		
		char *Value = 0;
		if (iArgIndex < argc - 1)
			Value = argv[iArgIndex+1];
		if (TestSetValueOpt(ArgName, Value))
			{
			iArgIndex += 2;
			continue;
			}
		CommandLineError("Invalid command line option \"%s\"\n", ArgName);
		}
	}

void IntOpt(const char *Name, int *ptrValue)
	{
	const char *strValue = ValueOpt(Name);
	if (strValue != 0)
		*ptrValue = atoi(strValue);
	}

void FloatOpt(const char *Name, double *ptrValue)
	{
	const char *strValue = ValueOpt(Name);
	if (strValue != 0)
		*ptrValue = atof(strValue);
	}

void Options()
	{
	fprintf(stderr,
"\n"
"\n"
"Basic options:\n"
"   -in <filename>         Sequence file to analyze (FASTA format).\n"
"   -out <filename>        Report file name (plain text).\n"
"   -seq <filename>        Save consensus sequences to this FASTA file.\n"
"   -trimseqs              Eliminate similar seqs from -seq file.\n"
"   -noinfo                Don't write help to report file.\n"
"   -quiet                 Don't write progress messages to stderr.\n"
"\n"
"Criteria for CRISPR detection, defaults in parentheses:\n"
"   -minarray <N>          Must be at least <n> repeats in array (3).\n"
"   -mincons <F>           Minimum conservation (0.9).\n"
"                            At least N repeats must have identity\n"
"                            >= F with the consensus sequence.\n"
"                            Value is in range 0 .. 1.0.\n"
"                            It is recommended to use a value < 1.0\n"
"                            because using 1.0 may suppress true\n"
"                            arrays due to boundary misidentification.\n"
"   -minrepeat <L>         Minimum repeat length (16).\n"
"   -maxrepeat <L>         Maximum repeat length (64).\n"
"   -minspacer <L>         Minimum spacer length (8).\n"
"   -maxspacer <L>         Maximum spacer length (64).\n"
"   -minrepeatratio <R>    Minimum repeat ratio (0.9).\n"
"   -minspacerratio <R>    Minimum spacer ratio (0.75).\n"
"                            'Ratios' are defined as minlength / maxlength,\n"
"                            thus a value close to 1.0 requires lengths to\n"
"                            be similar, 1.0 means identical lengths.\n"
"                            Spacer lengths sometimes vary significantly, so\n"
"                            the default ratio is smaller. As with -mincons,\n"
"                            using 1.0 is not recommended.\n"
"\n"
"Parameters for creating local alignments:\n"
"   -minhitlength <L>      Minimum alignment length (16).\n"
"   -minid <F>             Minimum identity (0.94).\n"
"\n"
"\n"
);
    }
