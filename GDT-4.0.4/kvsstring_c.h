// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: kvsstring_c.h 78 2010-09-24 20:14:30Z schumack $ 
#ifndef _kvsstring
#define _kvsstring
#include <kvstypes.h>

extern void    copy_array(float *a1, float *a2, int n);
extern char*   expandPosIntRange(char *input, char *output);
extern char*   fillCharArray(char *string, char character, int arr_length);
extern int     get_delim_field(char *string, char *field_string, char character, int position);
extern int     get_field(char *string, char *field_string, int position);
extern char*   get_field_number(char *string, char *field_string, int fieldNumber);
extern int     get_field_pos(char *string, int field, int position);
extern int     getargs(char *args[][80]);
extern void    getbody(char *input, char *output);
extern int     gsub(char *match, char *replacement, char *string);
extern int     gsubi(char *match, char *replacement, char *string);
extern char*   hexstopadedbins(char *string, char *binarystring, int numberofbits);
extern int     isdecint(char *string);
extern Boolean isfloat(char *string);
extern Boolean ishexint(char *string);
extern int     itobs(int i, char *bitstring);
extern void    legalstring(char *sinput, char *soutput, int maxlength, char fillerc);
extern int     lexidex(int y, int *array);
extern Boolean match_string(char *master, char *string, char option);
extern int     numfields(char *string);
extern char*   padCharArray(char *string, char character, int start, int arr_length);
extern char*   revstring(char *string, char *reversed);
extern char*   sGetWordPart(char *string, char *word, int position, int startOrEnd);
extern char*   sRemoveLeadingSpaces(char *parent, char *child);
extern char*   sRemoveSpaces(char *parent, char *child);
extern char*   sRemoveTrailingSpaces(char *parent, char *child);
extern char*   sRemoveWhiteSpace(char* parent, char* child);
extern char*   sRemoveTrailingZeros(char* parent, char* child);
extern char*   scolumn(char *string, char *control);
extern void    sdelcr(char *string);
extern char*   sdup(char *string);
extern int     sexpcmp(char *s1, char *s2);
extern int     sfind(char *string, char *line, char option);
extern int     sfindlast(char *string, char *line, char option);
extern char**  sget_argvs(char *string);
extern int     sgetword(char *string, char *word, int lim);
extern int     shexcmp(char *s1, char *s2);
extern char*   sintonly(char *parent, char *child);
extern int     sislower(char *string);
extern int     sisupper(char *string);
extern char    slastNonSpacec(char *string);
extern char*   smakealnum(char *parent, char *child);
extern int     snumc(char *string, char c);
extern int     snumcmp(char *s1, char *s2);
extern void    sortfile(char *options);
extern int     soverlap(char *string1, char *string2);
extern int     spatrep(char *string, char *pattern, char *replacement, int num_to_replace, Boolean ignore_case);
extern int     spatrepcol(char *string, char *pattern, char *replacement, int start_column, int num_to_replace, Boolean ignore_case, int number_replaced);
extern int     sposc(char *string, char c);
extern int     sposlastc(char *string, char c);
extern int     sposnextc(char *string, char c, int previousPosition);
extern char*   sprefix(char *prefix_string, char *string);
extern char*   sremovec(char *parent, char *child, char character);
extern char*   srepc(char *parent, char *child, char character, char replacement);
extern char*   sreppunct(char *parent, char *child, char replacement);
extern char*   sspacemin(char *parent, char *child);
extern char*   stolower(char *string, char *lowstring);
extern char*   stoupper(char *string, char *upstring);
extern void    str_qsort(struct sort_linest **v, int left, int right, int (*comparison)(), int start_field, int end_field);
extern char*   strip_version(char *input, char *output);
extern int     sub(char *match, char *replacement, char *string);
extern int     subi(char *match, char *replacement, char *string);
extern char*   substring(char *parent, char *child, int start, int end);
extern void    swap(void *v[], int i, int j);
extern int     wildInLine(char *wildToFind, char *line, Boolean ignore_case);
extern int     wild_match(char *string, char *wildstring, int ignore_case);
extern void    write_sortfile_help(void);

#endif

