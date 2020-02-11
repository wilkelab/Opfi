#include "pilercr.h"

void Credits()
	{
	static bool Displayed = false;
	if (Displayed)
		return;

	fprintf(stderr, "\n" PILERCR_LONG_VERSION "\n");
	fprintf(stderr, "http://www.drive5.com/piler\n");
	fprintf(stderr, "Written by Robert C. Edgar\n");
	fprintf(stderr, "This software is donated to the public domain.\n");
	fprintf(stderr, "Please visit web site for requested citation.\n\n");
	Displayed = true;
	}

void Usage()
	{
	Credits();
	fprintf(stderr, "\n");
	fprintf(stderr, "Basic usage:\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "   pilercr -in <sequence_file> -out <report_file>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Sequence file is in FASTA format.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "For advanced options, type pilercr -options.\n");
	fprintf(stderr, "\n");
	}
