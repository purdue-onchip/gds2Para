// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: stoupper.C 83 2010-09-24 20:14:35Z schumack $ 
#include <kvstypes.h>
#include <kvsstring_c.h>
#include <ctype.h>
/***** STOUPPER **************************************************************/
/* uses toupper repeatedly on a string   -kvs                                */
/*****************************************************************************/
char* stoupper(char* string, char* upstring)
{
    while (*string != '\0') 
    {
        if (islower(*string))  *upstring++ = toupper(*string++);
        else                   *upstring++ = *string++;
    }
    *upstring = '\0';       /* add null back on for end of string */
    return(upstring);
}

