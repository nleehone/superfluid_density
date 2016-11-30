#ifndef LOGGING_H
#define LOGGING_H

#include <time.h>
#include <stdio.h>

#define DEBUG 0

// See http://stackoverflow.com/questions/1644868/c-define-macro-for-debug-printing for detailed reasons for
// the form of this debug statement
#define debug_print(fmt, ...) \
              do { if (DEBUG) fprintf(fpLog, fmt, __VA_ARGS__); } while (0)
#define DTTMFMT "%Y-%m-%d %H:%M:%S "
#define DTTMSZ 21


extern char timebuff[DTTMSZ];
extern FILE* fpLog; // Shared log file

void openLog(char const* logFilename);
char* getTime();
void closeLog();

// Error handler for gsl errors
void error_handler(const char* reason, const char* file, int line, int gsl_errno);

#endif
