// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: get_field.C 77 2010-09-24 20:14:29Z schumack $ 
#include <kvstypes.h>
#include <kvsstring_c.h>
#include <ctype.h>

/***** GET FIELD *************************************************************/
/* Returns position of string passed in which ends the field starting at the */
/* position passed in with the integer variable. The field string is updated */
/* to the new field.                                                         */
/* ie.                                                                       */
/*      for (n = numfields(string), i=0; n>0; n--) {                         */
/*          i = get_field(string, name, i);                                  */
/*          printf("%s\n", name);                                            */
/*      }                                                                    */
/* or                                                                        */
/*      i=0                                                                  */
/*      while ((i=get_field(string, name, i)) != ERR) printf("%s\n", name);  */
/*                                                         -KVS              */
/*****************************************************************************/
int get_field(char* string, char* field_string, int position)
{ 
    int     i, j, k;

    for(i=position; (string[i] != '\0') &&   (isspace(string[i])); i++) ;  /* skip white space */
    for(j=i, k=0;  ((string[j] != '\0') && (!(isspace(string[j])))); j++, k++) field_string[k] = string[j];
    field_string[k] = '\0'; 

    if(string[i] == '\0') return(ERR);
    else return(j);
}

