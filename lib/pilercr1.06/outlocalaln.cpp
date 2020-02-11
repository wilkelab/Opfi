#include "pilercr.h"

void OutLocalAln(const char *A, unsigned LA, const char *B, unsigned LB,
  unsigned StartA, unsigned StartB, const char *Path_)
	{
	Log("A=%s\n", A);
	Log("B=%s\n", B);
	Log("Path=%s\n", Path_);

	unsigned i = StartA;
	const char *Path = Path_;
	while (char c = *Path++)
		{
		switch (c)
			{
		case 'M':
		case 'D':
			assert(i < LA);
			Log("%c", A[i++]);
			break;

		case 'I':
			Log("-");
			break;
		default:
			assert(false);
			}
		}
	Log("\n");

	unsigned j = StartB;
	Path = Path_;
	while (char c = *Path++)
		{
		switch (c)
			{
		case 'M':
		case 'I':
			assert(j < LB);
			Log("%c", B[j++]);
			break;

		case 'D':
			Log("-");
			break;
		default:
			assert(false);
			}
		}
	Log("\n");
	Log("\n");
	}
