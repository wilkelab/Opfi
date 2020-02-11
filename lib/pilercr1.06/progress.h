#ifndef PROGRESS_H
#define PROGRESS_H

void ShowProgress(bool Show);
void Progress(const char *Format, ...);
void ProgressNoRes(const char *Format, ...);
void ProgressStart(const char *Format, ...);
void ProgressStep(int StepIndex, int StepCount);
void ProgressDone();
void ProgressExit();

extern bool g_ShowProgress;

#endif // PROGRESS_H
