#include "pilercr.h"
#include <algorithm>

static bool LocalAcceptHit(const DPHit &Hit)
	{
	int Length = GetHitLength(Hit);
	if (Length < g_DraftMinRepeatLength || Length > g_DraftMaxRepeatLength)
		return false;

	//int StartSpace = std::min(Hit.aHi, Hit.bHi);
	//int EndSpace = std::max(Hit.aLo, Hit.bLo);
	//int Space = EndSpace - StartSpace + 1;
	//if (Space < g_DraftMinSpacerLength || Space > g_DraftMaxSpacerLength)
	//	return false;

	if (GlobalPosToContigIndex(Hit.aLo) != GlobalPosToContigIndex(Hit.bHi))
		return false;

	return true;
	}

bool AcceptHit(const DPHit &Hit)
	{
	bool Accepted = LocalAcceptHit(Hit);
	if (Accepted)
		++g_AcceptedHitCount;
	else
		++g_NotAcceptedHitCount;
	return Accepted;
	}
