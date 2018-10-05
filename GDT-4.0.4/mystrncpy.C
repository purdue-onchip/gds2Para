// Copyright 2014 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: mystrncpy.C 96 2013-05-22 20:50:29Z schumack $
/*****    mystrncpy  ********************************************************/
/* "safe" strcpy -kvs                                                       */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;
#include <gcc4.h>

char* mystrncpy(char *dest, char *src, size_t n)
{
    src[n] = '\0';
    strcpy(dest, src);
    return(dest);
}

