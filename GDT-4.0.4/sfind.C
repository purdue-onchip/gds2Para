// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: sfind.C 72 2010-09-24 20:14:24Z schumack $  
#include <kvstypes.h>
#include <kvsstring_c.h>
#include <ctype.h>
#include <string.h>
#include <gcc4.h>

/***** SFIND FUNCTION ********************************************************/
/* MODIFIED FROM MAW's search_string                                         */
/* trys to find a string in a line (string). Returns position in line that   */
/* string occurs or ERR if not found    -kvs                                 */
/* If 'I'or'i' option is used, case will be IGNORED when searching for       */
/* the string.                                                               */
/* If 'F' option is used - returns field in which string occurs.             */
/* If 'f' option is used - returns field in which string occurs &&           */
/* case IGNORED                                                              */
/*                                                                           */
/* note: 0 is the 1st position in the string.                                */
/*     ie.  if line == "hello there"                                         */
/*          sfind("The", line, 'i') == 6                                     */
/*     whereas sfind("bb", line) == ERR    if ERR is #defined as -1  -kvs    */
/*****************************************************************************/
int sfind(char* string, char* line, char option)
{
    int          line_len = strlen(line),
                 str_len,
                 i,
                 numwds,
                 ptr;
    Boolean      inword;
    char         c;

    switch (option) 
    {
      case 'i': case 'I':
        for (i=0; i < line_len; i++) 
        {
            if (line[i] == string[0] || line[i] == toupper(string[0]) || line[i] == tolower(string[0]) ) 
            {
                if (match_string(&line[i], string, 'm'))  return(i);
            }
        }
        break;
      case 'f':
        for (i=0; i < line_len; i++) 
        {
            if (line[i] == string[0] || line[i] == toupper(string[0]) || line[i] == tolower(string[0]) ) 
            {
                if (match_string(&line[i], string, 'm')) 
                {
                    for (numwds=ptr=0,inword=FALSE;(c = string[ptr]) != '\n' && c != '\0' && ptr<i;ptr++) 
                    {
                        if (c == ' ' || c == '\t')
                            inword = FALSE;
                        else if (inword == FALSE) 
                        {
                            inword = TRUE;
                            ++numwds;
                        }
                    }
                    return(numwds);
                }
            }
        }
        break;
      case 'F':  /* case sensitive */
        str_len = strlen(string);
        for (i=0; i < line_len; i++) 
        {
            if (line[i] == string[0]) 
            {
                if (strncmp(&line[i], string, str_len) == 0) 
                {
                    for (numwds=ptr=0,inword=FALSE;(c = string[ptr]) != '\n' && c != '\0' && ptr<i;ptr++) 
                    {
                        if (c == ' ' || c == '\t')
                            inword = FALSE;
                        else if (inword == FALSE) 
                        {
                            inword = TRUE;
                            ++numwds;
                        }
                    }
                    return(numwds);
                }
            }
        }
        break;
       default: 
        str_len = strlen(string);
        for (i=0; i < line_len; i++) 
        {
            if (line[i] == string[0]) 
            {
                if (strncmp(&line[i], string, str_len) == 0)  return(i);
            }
        } 
        break;
    }
    return(ERR); /* not found */
}

