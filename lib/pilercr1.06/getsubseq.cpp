#include "pilercr.h"

void GetSubseq(int Lo, int Hi, Seq &s)
	{
	s.clear();
	for (int i = Lo; i <= Hi; ++i)
		s.push_back(g_SeqQ[i]);
	}
