#include <time.h>

/*
 * Retourne la date sous forme d'une chaine de 26 caracteres
 */

void timedate (pc)
char *pc ;
{
    time_t seconds_since_epoch ;

    time (&seconds_since_epoch) ;
    strcpy (pc, ctime (&seconds_since_epoch)) ;
}

