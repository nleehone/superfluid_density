#include "logging.h"

char timebuff[DTTMSZ];
FILE* fpLog = NULL; // Shared log file

void openLog(char const* logFilename){
  if(!fpLog)
    fpLog = fopen(logFilename, "a");
}

char* getTime(){
  time_t t = time(0);
  strftime(timebuff, DTTMSZ, DTTMFMT, localtime(&t));
  return timebuff;
}

void closeLog(){
  if(fpLog){
    fclose(fpLog);
    fpLog = NULL;
  }
}

void error_handler(const char* reason, const char* file, int line, int gsl_errno){
  fprintf(fpLog, "gsl: %s:%d: errno:%d: %s\n", file, line, gsl_errno, reason);
  fflush(fpLog);
}
