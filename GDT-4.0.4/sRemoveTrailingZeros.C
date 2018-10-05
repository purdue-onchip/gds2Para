// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: sRemoveTrailingZeros.C 81 2010-09-24 20:14:33Z schumack $ 
#include <string.h>
#include <kvstypes.h>
#include <kvsstring_c.h>
#include <gcc4.h>

/***** sRemoveTrailingZeros *********************************************/
/* Remove zeros at the end of a real number string                 -kvs */
/************************************************************************/
char* sRemoveTrailingZeros(char* inString, char* outString)
{
    int i, done;
    done = 0;
    inString[LENGTHSSTRING] = '\0'; //safety
    strcpy(outString,inString);

    for(i=strlen(inString); (i>0) && (! done); i--) 
    {
        if (inString[i] == '.')
        {
            outString[i] = '\0';
            done = 1;
            break;
        }
        else if ((inString[i] != '0') && (inString[i] != '\0'))
        {
            done = 1;
            break;
        }
        else
        {
            outString[i] = '\0';
        }
    }
    return(outString);
}

