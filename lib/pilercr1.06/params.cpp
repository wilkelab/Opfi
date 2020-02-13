#include "params.h"

int g_DraftMinRepeatLength = 16;
int g_DraftMaxRepeatLength = 128;
int g_DraftMinSpacerLength = 8;
int g_DraftMaxSpacerLength = 256;
int g_DraftMinArraySize = 3;
double g_DraftMinHitRatio = 0.6;
double g_DraftMinSpacerRatio = 0.6;
double g_DraftMinCons = 0.95;
double g_DraftMinColCons = 0.90;

int g_MinRepeatLength = 16;
int g_MaxRepeatLength = 64;
int g_MinSpacerLength = 8;
int g_MaxSpacerLength = 64;
int g_MinArraySize = 3;
double g_MinRepeatRatio = 0.9;
double g_MinSpacerRatio = 0.75;
double g_MinCons = 0.90;

int g_MaxMSASeqs = 4096;
int g_HitPadding = 8;
int g_FlankSize = 10;
int g_MinOneDiff = 3;
int g_PileExpandMargin = 8;
int g_MinDiag = 8;
int g_MaxDPHitLen = 512;
int g_SpacerPadding = 16;

double g_ClusterMaxDist = 0.5;
double g_ColonThresh = 0.95;

TERMGAPS g_TermGaps = TERMGAPS_Ext;

bool g_ShowHits = false;
bool g_LogHits = false;
bool g_LogImages = false;
bool g_LogAlns = false;
bool g_FlushLog = false;
bool g_NoInfo = false;
