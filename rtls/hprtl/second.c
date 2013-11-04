#include <sys/times.h>
#include <unistd.h>                     /* POSIX P1003.1 */

/*
 * Retourne le temps CPU (en secondes) depuis le debut du programme
 * Note : le temps CPU est le cumul des temps systeme et usager
 */

void second (ph)
double *ph ;
{
    struct tms buf ;
    long hertz ;

    times (&buf) ;
    hertz = sysconf (_SC_CLK_TCK) ;

    *ph = (double) (buf.tms_utime + buf.tms_stime) / hertz ;
}

