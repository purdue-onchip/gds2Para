// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: kvstypes.h 79 2010-09-24 20:14:31Z schumack $ 
#ifndef __kvstypes__
#define __kvstypes__
#include <stdio.h>

static const char* kvstypes_hwhat = "@(#) $Id: kvstypes.h 79 2010-09-24 20:14:31Z schumack $ $Revision: 79 $ $Date: 2010-09-24 15:14:31 -0500 (Fri, 24 Sep 2010) $";

#define TRUE    1
#define FALSE   0
#define YES     1
#define NO      0
#define ERROR   (-1)
#define ERR     (-1)
#define READ    0
#define WRITE   1

#define LENGTHLSTRING 512
#define LENGTHLSTRINGINIT LENGTHLSTRING + 1

#define LENGTHSSTRING 80
#define LENGTHSSTRINGINIT LENGTHSSTRING + 1

typedef short   Boolean;
typedef char    stringS[LENGTHSSTRINGINIT];
typedef char    stringL[LENGTHLSTRINGINIT];
typedef double* coordArray;

#define SEEK_CUR 1

#define BUFSIZE     100
#define SIZE_FOREST 37

#define strEqual(s1, s2)   ( ! strcmp(s1, s2) )

typedef struct sort_linest 
{
    char*  strng;
    struct sort_linest* linkr;
} sort_linest;

typedef struct synonymForestMember 
{
    char*    name;
    char*    synonym;
    struct synonymForestMember* linkl;
    struct synonymForestMember* linkr;
}  synonymForestMember;

typedef struct GPO {  
    char*       name;
    char*       synonym;
    Boolean     Switch;
    int         Int;
    double      Double;
    struct GPO* headSynonym;
    struct GPO* linkl;
    struct GPO* linkr;
} GPO;

typedef GPO **GPOforest;
typedef GPO *GPOtree;

#endif

