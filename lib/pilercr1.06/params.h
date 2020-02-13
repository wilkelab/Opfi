#ifndef PARAMS_H

#include "types.h"

extern int g_DraftMinRepeatLength;
extern int g_DraftMaxRepeatLength;
extern int g_DraftMinSpacerLength;
extern int g_DraftMaxSpacerLength;
extern int g_DraftMinArraySize;

// For determing linkage:
extern double g_DraftMinHitRatio;
extern double g_DraftMinSpacerRatio;
extern int g_SpacerPadding;

extern double g_DraftMinCons;
extern double g_DraftMinColCons;
extern int g_MinOneDiff;

extern int g_MaxMSASeqs;
extern int g_PileExpandMargin;

extern int g_HitPadding;
extern int g_FlankSize;
extern int g_MinDiag;
extern int g_MaxDPHitLen;

extern double g_ClusterMaxDist;
extern double g_ColonThresh;
extern double g_MinPalindromeFractId;

extern TERMGAPS g_TermGaps;

extern bool g_ShowHits;

extern bool g_LogHits;
extern bool g_LogLinks;
extern bool g_LogImages;
extern bool g_LogAlns;
extern bool g_FlushLog;
extern bool g_NoInfo;

extern int g_MinRepeatLength;
extern int g_MaxRepeatLength;
extern int g_MinSpacerLength;
extern int g_MaxSpacerLength;
extern int g_MinArraySize;
extern double g_MinRepeatRatio;
extern double g_MinSpacerRatio;
extern double g_MinCons;

#endif	// PARAMS_H
