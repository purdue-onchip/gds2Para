// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: sRemoveSpaces.C 69 2010-09-24 20:14:21Z schumack $ 
 
//#include <string>
#include <string.h>
#include <kvstypes.h>
#include <kvsstring_c.h>
#include <gcc4.h>

/***** SREMOVESPACES ********************************************************/
/* Removes all occurances of a spaces & tabs in a string             -kvs   */
/****************************************************************************/
char* sRemoveSpaces(char* parent, char* child)
{
    int i, j, length;

    length = strlen(parent);
    for(i=0, j=0; i<length; i++) 
    {
        if (!(parent[i] == ' ' || parent[i] == '\t')) 
        {
            child[j] = parent[i];
            j++;
        }
    }
    child[j] = '\0';
    return(child);
}

