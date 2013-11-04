#include <stdio.h>
#include <sys/time.h>
char timedate (ival,len)
char *ival;
int len;
{ 
int ierr, gettimeofday();
char *ctime();
struct timeval Tp;
struct timezone Tzp;
     if ((ierr=gettimeofday(&Tp,&Tzp))!=0)printf("Err gettimeofday!\n");
     strcpy(ival,ctime(&Tp.tv_sec));
     ival[24] = ival[25] = '\b';
}
