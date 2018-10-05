/*
 ******************************************************************************
 * Copyright 2003-2014 by Ken Schumack
 *
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author:
 *           Ken Schumack (Schumack@cpan.org)
 *
 * Suggestions for improvements may be sent to Ken Schumack.
 ******************************************************************************
*/
// @(#) $Id: gdt2gds.C 96 2014-12-04 22:03:27Z schumack $
#include <stdio.h>
#include <math.h>
#include <gdsStream.h>
#include <kvsstring_c.h>
//#include <iostream.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
//#include <ostream.h>
#include <ostream>
using namespace std;
#include <gcc4.h>
#include <sourceForge.h>

#define LINELENGTH 1000000 //gdt lines can be this long
#define ERRORSTRING "ERROR:  "
#define NUMPROPS 30 //can handle this many props on an element

extern char* mystrncpy(char *dest, char *src, size_t n);
extern char* sRemoveTrailingZeros(char* parent, char* child);
extern void sStdWhiteSpace(char* parent, int isText);
extern void rmLeadingAndTrailingSpace(char* parent);
extern char* sRemoveWhiteSpace(char* parent, char* child);
extern char* sUnescapedString(char* parent, char* child);
extern void print_help();
extern int getUnescapedPropString(char* parent, char* child, int position);
extern int getUnescapedTextString(char* parent, char* child, int position);
extern int getField(char* string, char* field_string, int position);
extern int sIsPosInt(char* string);

//////////////////////////////////////////////////////////////////////////////
// M A I N
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    int     i,j,k,stringLen,len,pos=0,layer,vertIndex,lineCnt,
            firstQuotePos=0,
            col,row,
            precision = 6,
            inPrecision = 0,
            dataType,pathType,fontType,textType,boxType,nodeType,ref,gds2Version,
            inText,inXy,inCr,
            propIndex,propNumArray[NUMPROPS],
            mYear=0,mMon=0,mDay=0,mHour=0,mMin=0,mSec=0,
            aYear=0,aMon=0,aDay=0,aHour=0,aMin=0,aSec=0;
    double  width=0.0,
            bgnExtn=0.0,
            endExtn=0.0,
            mag=1.0,
            angle=0.0,
            dbUnits=1e-9,
            uUnits=0.001,
            xd1=0.0,
            xd2=0.0,
            xd3=0.0,
            yd1=0.0,
            yd2=0.0,
            yd3=0.0,
            xArray[MAXPAIRSXY], yArray[MAXPAIRSXY];
    stringL inputFile, outputFile, tmpFile;
    Boolean haveInputName  = FALSE,
            haveOutputName = FALSE,
            inInFileName   = FALSE,
            inOutFileName  = FALSE,
            inLib          = FALSE,
            debug          = FALSE,
            done           = FALSE,
            grabName       = FALSE,
            sawXy          = FALSE,
            inStr          = FALSE,
            useStdout      = FALSE;
    char    field[MAXPAIRSXY],
            field2[MAXPAIRSXY],
            text[MAXPAIRSXY],
            cellName[LENGTHLSTRINGINIT],
            libName[LENGTHLSTRINGINIT],
            textJust[3],
            propValueArray[NUMPROPS][LENGTHLSTRING],
            lineWithSpace[LINELENGTH],
            line[LINELENGTH];
    FILE    *fpIn=0;

    // next line uses '@ ( # )' so the unix/linux 'what' command works on the executable
    cerr << "# ** gdt2gds @(#) Version " << sourceForgeVersion << " $Id: gdt2gds.C 96 2014-10-28 23:03:27Z schumack $ ** #" << endl; //using subversion props
    textJust[2] = '\0';
    strcpy(inputFile, "");
    strcpy(outputFile, "");
    for(i=1; i<argc; i++)
    {
        if((argv[i])[0] == '-')
        {
            if ((inInFileName) && (! strcmp("-", argv[i])))
            {
                haveInputName = 1;
            }
            else if ((inOutFileName) && (! strcmp("-", argv[i])))
            {
                haveOutputName = TRUE;
                outputFile[0] = '\0';
                useStdout = TRUE;
            }
            else if (match_string("-debug", argv[i], 'm'))
            {
                debug = TRUE;
            }
            else if (match_string("-help", argv[i], 'm'))
            {
                print_help();
                exit(0);
            }
            else if (match_string("-infilename", argv[i], 'm'))
            {
                inInFileName = TRUE;
            }
            else if (match_string("-outfilename", argv[i], 'm'))
            {
                inOutFileName = TRUE;
            }
            else if (match_string("-precision", argv[i], 'm'))
            {
                inPrecision = 1;
            }
            else if (match_string("-version", argv[i], 'm') || match_string("-swversion", argv[i], 'm'))
            {
                exit(0);   // version already printed */
            }
        }
        else
        {
            if (inInFileName)
            {
                haveInputName = TRUE;
                inInFileName = FALSE;
                strcpy(inputFile, argv[i]);
            }
            else if (inOutFileName)
            {
                haveOutputName = TRUE;
                inOutFileName = FALSE;
                strcpy(outputFile, argv[i]);
            }
            else if (inPrecision)
            {
                inPrecision = 0;
                precision = atoi(argv[i]);
            }
            else if (! haveInputName)
            {
                haveInputName = TRUE;
                strcpy(inputFile, argv[i]);
            }
            else if (! haveOutputName)
            {
                haveOutputName = TRUE;
                strcpy(outputFile, argv[i]);
            }
        }
    }

    if (! haveInputName)
    {
        fpIn = stdin;
    }
    else if((fpIn = fopen(inputFile,"r")) == NULL)
    {
        fprintf(stderr,"%sunable to read file %s.\n", ERRORSTRING, inputFile);
        exit(1);
    }

    if (! haveOutputName)
    {
        useStdout = TRUE;
        haveOutputName = TRUE;
    }

    if ((! useStdout) && strEqual(inputFile,outputFile))
    {
        sprintf(outputFile,"%s.gds",outputFile);
        fprintf(stderr,"Warning ** you gave the input filename for the output file.  Will use %s instead\n", outputFile);
    }
    if (debug) printf("DEBUG: ON\n");
    // ************* end of command line and menu stuff ***************************
    /*
    # Key: <required> [optional]
    # File format:
    # gds2{<ver>
    # m=<modificationTimeStamp> a=<accessTimeStamp>
    # lib '<libName>' <userUnits> <dataUnits>
    # <cellDefine>
    # }
    # - - - - -
    # cellDefine is one of more of:
    # cell {c=<creationTimeStamp> m=<modificationTimeStamp> '<cellName>'
    # <cellStuff>
    # }
    # - - - - - - - - -
    # cellStuff is one or more of:
    # boundary:
    # b{<layer> [dt<dataType>] xy(<xyList>)}
    #
    # path:
    # p{<layer> [dt<dataType>] [pt<pathType>] [w<width>] [bx<real>] [ex<real>] xy(<xyList>)}
    #
    # text:
    # t{<layer> [tt<textType>] [f<fontType>] [<textJust>] [pt<pathType>] [fx] [m<magification>] [a<angle>] xy(<xyList>) <'text'> }
    #
    # sref:
    # s{<'cellName'> [fx] [a<angle>] xy(<xyList>)}
    #
    # aref:
    # a{<'cellName'> [fx] [a<angle>] cr(<columns> <rows>) xy(<xyList>)}
    #
    # box:
    # x{<layer> [bt<boxType>] xy(<xyList>)}
    #
    # node:
    # n{<layer> [nt<nodeType>] xy(<xyList>)}
    #
    # - - - - - - - - -
    # property : pr{<propAttr> <'propValue'>}
    # <text> : ASCII String
    # textJust: one of mc bl br bc mr ml tr tc tl
    #            Where m is middle b is bottom t is top
    #                  c is center l is left r is right
    # <propAttr> : a 2 byte (small) integer
    # <propValue> : ASCII String
    # # as first character on a line is a comment
    */
    if (fgets(line, LINELENGTH, fpIn)) // gds2{600
    {
        sRemoveWhiteSpace(line,line);
        stringLen = strlen(line);
        if (strncmp("gds2",line,4))
        {
            fprintf(stderr,"%sDoes not appear to be a GDT file. Did not find 'gds2' as 1st characters on 1st line.\n",ERRORSTRING);
            exit(4);
        }
        inLib = TRUE;
        for(i=0, j=0; i<stringLen; i++)
        {
            line[j] = line [i];
            j++;
            if (line[i] == '{') j = 0;
        }
        line[j] = '\0';
        i = atoi(line);
        if ((i>0) && (i<=600))
        {
            gds2Version = i;
        }
        else
        {
            fprintf(stderr,"%sInvalid (%s) gds2 version\n",ERRORSTRING,line);
            exit(44);
        }
    }
    else
    {
        fprintf(stderr,"%sUnable to read line one of %s\n",ERRORSTRING,inputFile);
        exit(5);
    }
    //////////////////////////////////////////////////////////////////////

    if (fgets(line, LINELENGTH, fpIn)) // m=2004-05-11 14:02:34 a=2004-05-11 14:02:34
    {
        sRemoveWhiteSpace(line,line); //m=2004-05-1114:02:34a=2004-05-1114:02:34
        stringLen = strlen(line);
        if (strncmp("m=",line,2))
        {
            fprintf(stderr,"%sDoes not appear to be a valid GDT file. Did not find 'm=' as 1st characters on 2nd line.\n",ERRORSTRING);
            exit(4);
        }
        len = strlen(line);
        for(i=0, j=0; i<len; i++) //remove all non integers
        {
            if (isdigit(line[i]))
            {
                line[j] = line[i];
                j++;
            }
        }
        line[j] = '\0';
        if (strlen(line) == 28)
        {
            strncpy(field,&line[0],4);
            field[4] = '\0';
            mYear = atoi(field);
            strncpy(field,&line[4],2);
            field[2] = '\0';
            mMon = atoi(field);
            strncpy(field,&line[6],2);
            field[2] = '\0';
            mDay = atoi(field);
            strncpy(field,&line[8],2);
            field[2] = '\0';
            mHour = atoi(field);
            strncpy(field,&line[10],2);
            field[2] = '\0';
            mMin = atoi(field);
            strncpy(field,&line[12],2);
            field[2] = '\0';
            mSec = atoi(field);
            strncpy(field,&line[14],4);
            field[4] = '\0';
            aYear = atoi(field);
            strncpy(field,&line[18],2);
            field[2] = '\0';
            aMon = atoi(field);
            strncpy(field,&line[20],2);
            field[2] = '\0';
            aDay = atoi(field);
            strncpy(field,&line[22],2);
            field[2] = '\0';
            aHour = atoi(field);
            strncpy(field,&line[24],2);
            field[2] = '\0';
            aMin = atoi(field);
            strncpy(field,&line[26],2);
            field[2] = '\0';
            aSec = atoi(field);
            if (mYear > 1900) mYear -= 1900;
            if (aYear > 1900) aYear -= 1900;
        }
        else
        {
            fprintf(stderr,"%sDoes not appear to be a valid GDT file. Did not find valid 'm=' 2nd line.\n",ERRORSTRING);
            exit(5);
        }
    }
    else
    {
        fprintf(stderr,"%sUnable to read line one of %s\n",ERRORSTRING,inputFile);
        exit(5);
    }
    //////////////////////////////////////////////////////////////////////

    GDSFILE gds2file(outputFile, WRITE);
    lineCnt=2;
    while (fgets(line, LINELENGTH, fpIn) != NULL)
    {
        mystrncpy(lineWithSpace,line,LINELENGTH);
        lineCnt++;
        rmLeadingAndTrailingSpace(line);
        stringLen = strlen(line);
        propIndex = -1;
        if (line[0] == 'b') //b{12 xy(90.825 28.9 96.725 28.9 96.725 31.65 90.825 31.65)} --- boundary
        {
            sawXy = FALSE;
            layer = 0;
            dataType = 0;
            for(done=FALSE, i=0; (! done) && (i <= stringLen); i++)
            {
                if (line[i] == '{' )
                {
                    done = TRUE;
                }
                line[i] = ' ';
            }
            if (done)
            {
                sStdWhiteSpace(line,0);
            }
            else
            {
                fprintf(stderr,"%sBad boundary '%s' lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            pos = 0;
            pos = getField(line, field, pos); // layer
            if (sIsPosInt(field))
            {
                layer = atoi(field);
            }
            else
            {
                fprintf(stderr,"%sBad boundary '%s' bad layer. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            dataType = 0;
            vertIndex = 0;
            inXy = FALSE;
            while ((pos = getField(line, field, pos)) > 0)
            {
                if ((field[0] == ')') || (field[0] == '}'))
                {
                    inXy = FALSE;
                }
                if (inXy)
                {
                    if (vertIndex < 1)
                    {
                        xd1 = atof(field);
                        xd2 = xd1;
                    }
                    else
                    {
                        xd2 = atof(field);
                    }
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        if (vertIndex < 1)
                        {
                            yd1 = atof(field);
                            yd2 = yd1;
                        }
                        else
                        {
                            yd2 = atof(field);
                        }
                    }
                    else
                    {
                        fprintf(stderr,"%sBad boundary '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                    xArray[vertIndex] = xd2;
                    yArray[vertIndex] = yd2;
                    vertIndex++;
                }
                else if (! strncmp("dt",field,2))
                {
                    strcpy(field2,&field[2]);
                    dataType = atoi(field2);
                }
                else if (field[0] == '}') // end of boundary
                {
                    if ((xArray[0] != xArray[(vertIndex - 1)]) || (yArray[0] != yArray[(vertIndex - 1)])) //boundary closure
                    {
                        xArray[vertIndex] = xd1;
                        yArray[vertIndex] = yd1;
                        vertIndex++;
                    }
                    if (propIndex >= 0)
                    {
                        gds2file.putBndDbl(layer, dataType, xArray, yArray, vertIndex, propIndex, propNumArray, propValueArray, uUnits );
                    }
                    else
                    {
                        gds2file.putBndDbl(layer, dataType, xArray, yArray, vertIndex, uUnits);
                    }
                }
                else if (
                    (! strncmp("xy(",field,3)) &&
                    (! sawXy)
                )
                {
                    inXy = TRUE;
                    sawXy = TRUE;
                }
                else if (
                          (! strncmp("pr{",field,3)) ||
                          (! strncmp("pr(",field,3))
                )
                {
                    propIndex++;
                    i = atoi(&field[3]);
                    propNumArray[propIndex] = i; //PROPATTR
                    propValueArray[propIndex][0] = '\0';
                    pos = getUnescapedPropString(line,field,++pos);
                    strcpy(propValueArray[propIndex],field); //PROPVALUE
                }
                else if ((field[0] != ')') && (field[0] != '}'))
                {
                    fprintf(stderr,"%sBad boundary '%s'. Field='%s' LineNum=%d\n",ERRORSTRING,line,field,lineCnt);
                    exit(5);
                }
            }
        }
        else if (line[0] == 'p') //p{16 pt2 xy(-0.425 -0.1 -0.425 -0.425 1.525 -0.425 1.525 0.425 -0.425 0.425 -0.425 0.1)} -- path
        { //p{<layer> [dataType] [pt<pathType>] [w<width>] [bx<real>] [ex<real>] xy(<xyList>)}
            layer = 0;
            for(done=FALSE, i=0; (! done) && (i <= stringLen); i++)
            {
                if (line[i] == '{' )
                {
                    done = TRUE;
                }
                line[i] = ' ';
            }
            if (done)
            {
                sStdWhiteSpace(line,0);
            }
            else
            {
                fprintf(stderr,"%sBad path '%s' lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            pos = 0;
            pos = getField(line, field, pos); // layer
            if (sIsPosInt(field))
            {
                layer = atoi(field);
            }
            else
            {
                fprintf(stderr,"%sBad path '%s' bad layer. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            dataType = 0;
            pathType = 0;
            vertIndex = 0;
            inXy = FALSE;
            width = 0.0;
            bgnExtn = 0.0;
            endExtn = 0.0;
            while ((pos = getField(line, field, pos)) > 0)
            {
                if ((field[0] == ')') || (field[0] == '}'))
                {
                    inXy = FALSE;
                }

                if (inXy)
                {
                    xd2 = atof(field);
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        yd2 = atof(field);
                    }
                    else
                    {
                        fprintf(stderr,"%sBad boundary '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                    xArray[vertIndex] = xd2;
                    yArray[vertIndex] = yd2;
                    vertIndex++;
                }
                else if (! strncmp("dt",field,2))
                {
                    strcpy(field2,&field[2]);
                    dataType = atoi(field2);
                }
                else if (! strncmp("pt",field,2))
                {
                    strcpy(field2,&field[2]);
                    pathType = atoi(field2);
                }
                else if (! strncmp("w",field,1))
                {
                    strcpy(field2,&field[1]);
                    width = atof(field2);
                }
                else if (! strncmp("bx",field,2))
                {
                    strcpy(field2,&field[2]);
                    bgnExtn = atof(field2);
                }
                else if (! strncmp("ex",field,2))
                {
                    strcpy(field2,&field[2]);
                    endExtn = atof(field2);
                }
                else if (field[0] == '}') // end of path
                {
                    if (propIndex >= 0)
                    {
                        gds2file.putPathDbl(layer, dataType, pathType, width, bgnExtn, endExtn, xArray, yArray, vertIndex, propIndex, propNumArray, propValueArray, uUnits);
                    }
                    else
                    {
                        gds2file.putPathDbl(layer, dataType, pathType, width, bgnExtn, endExtn, xArray, yArray, vertIndex, uUnits);
                    }
                }
                else if (! strncmp("xy(",field,3)) inXy = TRUE;
                else if (
                  (! strncmp("pr{",field,3)) ||
                  (! strncmp("pr(",field,3))
                )
                {
                    propIndex++;
                    i = atoi(&field[3]);
                    propNumArray[propIndex] = i; //PROPATTR
                    propValueArray[propIndex][0] = '\0';
                    pos = getUnescapedPropString(line,field,++pos);
                    strcpy(propValueArray[propIndex],field);
                }
                else if ((field[0] != ')') && (field[0] != '}'))
                {
                    fprintf(stderr,"%sBad path '%s'. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                    exit(5);
                }
            }
        }
        else if (line[0] == 's') // s{<'cellName'> [fx] [a<angle>] xy(<xyList>)}
        { //s{'PECLCMLINDIFFJPLS33__P_DR_DV' fx a180 xy(5.62 0)}
            inXy = FALSE;
            for(done=FALSE, i=0; (! done) && (i <= stringLen); i++)
            {
                if (line[i] == '{' )
                {
                    done = TRUE;
                }
                line[i] = ' ';
            }
            if (done)
            {
                sStdWhiteSpace(line,0);
            }
            else
            {
                fprintf(stderr,"%sBad sref '%s' lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            pos = 0;
            ref = 0;
            mag = 1.0;
            angle = 0.0;
            pos = getField(line, field, pos); // layer
            if ((field[0] == '\'') && (field[strlen(field) - 1] == '\''))
            {
                strcpy(cellName,&field[1]);
                cellName[strlen(cellName) - 1] = '\0'; //remove the last '
            }
            else
            {
                fprintf(stderr,"%sBad sref name '%s'. lineNum=%d\n",ERRORSTRING,field,lineCnt);
                exit(5);
            }
            while ((pos = getField(line, field, pos)) > 0)
            {
                if (field[0] == ')')
                {
                    inXy = FALSE;
                }
                else if (field[0] == '}') // end of sref
                {
                    if (propIndex >= 0)
                    {
                        gds2file.putSref(cellName,ref,mag,angle,xd1,yd1,propIndex,propNumArray,propValueArray,uUnits);
                    }
                    else
                    {
                        gds2file.putSref(cellName,ref,mag,angle,xd1,yd1,uUnits);
                    }
                }
                else if (inXy)
                {
                    xd1 = atof(field);
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        yd1 = atof(field);
                        inXy = FALSE;
                        continue;
                    }
                    else
                    {
                        fprintf(stderr,"%sBad sref '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                }
                else if (! strncmp("xy(",field,3)) inXy = TRUE;
                else if (! strncmp("fx",field,2))
                {
                    ref = 1;
                }
                else if (! strncmp("a",field,1))
                {
                    strcpy(field2,&field[1]);
                    angle = atof(field2);
                }
                else if (! strncmp("m",field,1))
                {
                    strcpy(field2,&field[1]);
                    mag = atof(field2);
                }
                else if (
                          (! strncmp("pr{",field,3)) ||
                          (! strncmp("pr(",field,3))
                )
                {
                    propIndex++;
                    i = atoi(&field[3]);
                    propNumArray[propIndex] = i; //PROPATTR
                    propValueArray[propIndex][0] = '\0';
                    pos = getUnescapedPropString(line,field,++pos);
                    strcpy(propValueArray[propIndex],field); //PROPVALUE
                }
            }
        }
        else if (line[0] == 'a') // a{<'cellName'> [fx] [a<angle>] cr(<columns> <rows>) xy(<xyList>)}
        { //a{'CD00' fx cr(22 1) xy(99.49 -5.095 108.51 -5.095 99.49 -5.505)}
            inXy = FALSE;
            for(done=FALSE, i=0; (! done) && (i <= stringLen); i++)
            {
                if (line[i] == '{' )
                {
                    done = TRUE;
                }
                line[i] = ' ';
            }
            if (done)
            {
                sStdWhiteSpace(line,0);
            }
            else
            {
                fprintf(stderr,"%sBad sref '%s' lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            pos = 0;
            ref = 0;
            mag = 1.0;
            angle = 0.0;
            col = 1;
            row = 1;
            pos = getField(line, field, pos); // layer
            if ((field[0] == '\'') && (field[strlen(field) - 1] == '\''))
            {
                field[LENGTHLSTRING] = '\0'; //safety
                strcpy(cellName,&field[1]);
                cellName[strlen(cellName) - 1] = '\0'; //remove the last '
            }
            else
            {
                fprintf(stderr,"%sBad sref name '%s'. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            while ((pos = getField(line, field, pos)) > 0)
            {
                if (field[0] == ')')
                {
                    inXy = FALSE;
                    inCr = FALSE;
                }
                else if (field[0] == '}') // end of aref
                {
                    if (propIndex >= 0)
                    {
                        gds2file.putAref(cellName,ref,mag,angle,col,row,xd1,yd1,xd2,yd2,xd3,yd3,propIndex,propNumArray,propValueArray,uUnits);
                    }
                    else
                    {
                        gds2file.putAref(cellName,ref,mag,angle,col,row,xd1,yd1,xd2,yd2,xd3,yd3,uUnits);
                    }
                }
                else if (inXy)
                {
                    xd1 = atof(field);
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        yd1 = atof(field);
                        if ((pos = getField(line, field, pos)) > 0)
                        {
                            xd2 = atof(field);
                            if ((pos = getField(line, field, pos)) > 0)
                            {
                                yd2 = atof(field);
                                if ((pos = getField(line, field, pos)) > 0)
                                {
                                    xd3 = atof(field);
                                    if ((pos = getField(line, field, pos)) > 0)
                                    {
                                        yd3 = atof(field);
                                    }
                                    else
                                    {
                                        fprintf(stderr,"%sBad aref '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                                        exit(5);
                                    }
                                }
                                else
                                {
                                    fprintf(stderr,"%sBad aref '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                                    exit(5);
                                }
                            }
                            else
                            {
                                fprintf(stderr,"%sBad aref '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                                exit(5);
                            }
                        }
                        else
                        {
                            fprintf(stderr,"%sBad aref '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                            exit(5);
                        }
                        inXy = FALSE;
                        continue;
                    }
                    else
                    {
                        fprintf(stderr,"%sBad aref '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                }
                else if (inCr)
                {
                    col = atoi(field);
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        row = atoi(field);
                    }
                    else
                    {
                        fprintf(stderr,"%sBad aref '%s' bad cr() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                }
                else if (! strncmp("xy(",field,3)) inXy = TRUE;
                else if (! strncmp("cr(",field,3)) inCr = TRUE;
                else if (! strncmp("fx",field,2))
                {
                    ref = 1;
                }
                else if (! strncmp("a",field,1))
                {
                    strcpy(field2,&field[1]);
                    angle = atof(field2);
                }
                else if (! strncmp("m",field,1))
                {
                    strcpy(field2,&field[1]);
                    mag = atof(field2);
                }
                else if (
                          (! strncmp("pr{",field,3)) ||
                          (! strncmp("pr(",field,3))
                )
                {
                    propIndex++;
                    i = atoi(&field[3]);
                    propNumArray[propIndex] = i; //PROPATTR
                    propValueArray[propIndex][0] = '\0';
                    pos = getUnescapedPropString(line,field,++pos);
                    strcpy(propValueArray[propIndex],field); //PROPVALUE
                }
            }
        }
        else if (line[0] == 't') // t{<layer> [textType] [f<fontType>] <textJust> [pt<pathType>] [m<magification>] [a<angle>] xy(<xyList>) <'text'> }
        {
            layer = 0;
            for(done=FALSE, i=0; (! done) && (i <= stringLen); i++)
            {
                if (line[i] == '{' )
                {
                    done = TRUE;
                }
                line[i] = ' ';
            }
            if (done)
            {
                sStdWhiteSpace(line,1);
            }
            else
            {
                fprintf(stderr,"%sBad text '%s' lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            pos = 0;
            pos = getField(line, field, pos); // layer
            if (sIsPosInt(field))
            {
                layer = atoi(field);
            }
            else
            {
                fprintf(stderr,"%sBad text '%s' bad layer. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            fontType = 0;
            textType = 0;
            pathType = 0;
            mag = 1.0;
            ref = 0;
            angle=0.0;
            sprintf(textJust,"tl");
            inXy = FALSE;
            inText = FALSE;
            text[0] = '\0';
            width = 0.0;
            while ((pos = getField(line, field, pos)) > 0)
            {
                if (field[0] == ')')
                {
                    inXy = FALSE;
                } //{
                else if (field[0] == '}') // end of text
                {
                    inText = FALSE; //reset
                    for(i=0; i<=strlen(text); i++)
                    {
                        if (text[i] == '\r')
                        {
                            text[i] = '\n';
                        }
                    }
                    if (propIndex >= 0)
                    {
                        gds2file.putText(layer, textType, fontType, textJust, pathType, width, ref, mag, angle, xd2, yd2, text, propIndex, propNumArray, propValueArray, uUnits);
                    }
                    else
                    {
                        gds2file.putText(layer, textType, fontType, textJust, pathType, width, ref, mag, angle, xd2, yd2, text, uUnits);
                    }
                }
                else if (inXy)
                {
                    xd2 = atof(field);
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        yd2 = atof(field);
                        inXy = FALSE;
                        continue;
                    }
                    else
                    {
                        fprintf(stderr,"%sBad text '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                }
                else if (
                          (! strncmp("pr{",field,3)) ||
                          (! strncmp("pr(",field,3))
                )
                {
                    propIndex++;
                    i = atoi(&field[3]);
                    propNumArray[propIndex] = i; //PROPATTR
                    propValueArray[propIndex][0] = '\0';
                    pos = getUnescapedPropString(line,field,++pos);
                    strcpy(propValueArray[propIndex],field); //PROPVALUE
                }
                else if (field[0] == '\'') inText = TRUE;
                else if (! inText)
                {
                    if (! strncmp("xy(",field,3)) inXy = TRUE;
                    else if (! strncmp("pt",field,2))
                    {
                        strcpy(field2,&field[2]);
                        pathType = atoi(field2);
                    }
                    else if (! strncmp("fx",field,2))
                    {
                        ref = 1;
                    }
                    else if (! strncmp("f",field,1))
                    {
                        strcpy(field2,&field[1]);
                        fontType = atoi(field2);
                    }
                    else if (
                            (! strncmp("mc",field,2)) ||
                            (! strncmp("bl",field,2)) ||
                            (! strncmp("br",field,2)) ||
                            (! strncmp("bc",field,2)) ||
                            (! strncmp("mr",field,2)) ||
                            (! strncmp("ml",field,2)) ||
                            (! strncmp("tr",field,2)) ||
                            (! strncmp("tc",field,2)) ||
                            (! strncmp("tl",field,2))
                        )
                    {
                        strncpy(textJust,field,2);
                    }
                    else if (! strncmp("tt",field,2))
                    {
                        strcpy(field2,&field[2]);
                        textType = atoi(field2);
                    }
                    else if (! strncmp("a",field,1))
                    {
                        strcpy(field2,&field[1]);
                        angle = atof(field2);
                    }
                    else if (! strncmp("m",field,1))
                    {
                        strcpy(field2,&field[1]);
                        mag = atof(field2);
                    }
                    else
                    {
                        fprintf(stderr,"%sBad text \"%s\" lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        if (debug) printf("DEBUG %d: field=%s\n",__LINE__,field);
                        exit(5);
                    }
                }
                if (inText)
                {
                    firstQuotePos = 1;
                    while (lineWithSpace[firstQuotePos] != '\'') { ++firstQuotePos; }
                    if (debug) printf("\nDEBUG %d: %s pos=%d firstQuotePos=%d",__LINE__,line,pos,firstQuotePos);
                    pos = getUnescapedTextString(line,text,firstQuotePos);
                    if (debug) printf(" text=%s\n",text);
                    inText = FALSE;
                }
            }
        }
        else if (! strncmp("cell",line,4)) //cell{c=1998-08-17 14:31:10 m=1998-08-17 14:33:47 'CC1804_R1'   #}
        {
            inStr = TRUE;
            grabName = FALSE;
            j=0;
            k=0;
            for(i=3;i<strlen(line);i++)
            {
                if (grabName) cellName[k++] = line[i];
                if (isdigit(line[i])) field2[j++] = line[i];
                if (line[i] == '\'') grabName = TRUE;
            }
            cellName[k - 1] = '\0';
            field2[j] = '\0';
            if (strlen(line) >= 28)
            {
                strncpy(field,&field2[0],4);
                field[4] = '\0';
                mYear = atoi(field);
                strncpy(field,&field2[4],2);
                field[2] = '\0';
                mMon = atoi(field);
                strncpy(field,&field2[6],2);
                field[2] = '\0';
                mDay = atoi(field);
                strncpy(field,&field2[8],2);
                field[2] = '\0';
                mHour = atoi(field);
                strncpy(field,&field2[10],2);
                field[2] = '\0';
                mMin = atoi(field);
                strncpy(field,&field2[12],2);
                field[2] = '\0';
                mSec = atoi(field);
                strncpy(field,&field2[14],4);
                field[4] = '\0';
                aYear = atoi(field);
                strncpy(field,&field2[18],2);
                field[2] = '\0';
                aMon = atoi(field);
                strncpy(field,&field2[20],2);
                field[2] = '\0';
                aDay = atoi(field);
                strncpy(field,&field2[22],2);
                field[2] = '\0';
                aHour = atoi(field);
                strncpy(field,&field2[24],2);
                field[2] = '\0';
                aMin = atoi(field);
                strncpy(field,&field2[26],2);
                field[2] = '\0';
                aSec = atoi(field);
            }
            else
            {
                fprintf(stderr,"%sDoes not appear to be a valid GDT file. Did not find valid 'cell{' line.\n",ERRORSTRING);
                exit(5);
            }
            if (mYear > 1900) mYear -= 1900;
            if (aYear > 1900) aYear -= 1900;
            gds2file.beginStr(cellName,mYear,mMon,mDay,mHour,mMin,mSec,aYear,aMon,aDay,aHour,aMin,aSec);
        } //{
        else if (line[0] == 'x') //x{12 xy(90.825 28.9 96.725 28.9 96.725 31.65 90.825 31.65)} --- box
        {
            sawXy = FALSE;
            layer = 0;
            boxType = 0;
            for(done=FALSE, i=0; (! done) && (i <= stringLen); i++)
            {
                if (line[i] == '{' )
                {
                    done = TRUE;
                }
                line[i] = ' ';
            }
            if (done)
            {
                sStdWhiteSpace(line,0);
            }
            else
            {
                fprintf(stderr,"%sBad box '%s' lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            pos = 0;
            pos = getField(line, field, pos); // layer
            if (sIsPosInt(field))
            {
                layer = atoi(field);
            }
            else
            {
                fprintf(stderr,"%sBad box '%s' bad layer. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            boxType = 0;
            vertIndex = 0;
            inXy = FALSE;
            while ((pos = getField(line, field, pos)) > 0)
            {
                if ((field[0] == ')') || (field[0] == '}'))
                {
                    inXy = FALSE;
                }
                if (inXy)
                {
                    if (vertIndex < 1)
                    {
                        xd1 = atof(field);
                        xd2 = xd1;
                    }
                    else
                    {
                        xd2 = atof(field);
                    }
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        if (vertIndex < 1)
                        {
                            yd1 = atof(field);
                            yd2 = yd1;
                        }
                        else
                        {
                            yd2 = atof(field);
                        }
                    }
                    else
                    {
                        fprintf(stderr,"%sBad box '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                    xArray[vertIndex] = xd2;
                    yArray[vertIndex] = yd2;
                    vertIndex++;
                }
                else if (! strncmp("dt",field,2))
                {
                    strcpy(field2,&field[2]);
                    boxType = atoi(field2);
                }
                else if (field[0] == '}') // end of box
                {
                    if ((xArray[0] != xArray[(vertIndex - 1)]) || (yArray[0] != yArray[(vertIndex - 1)])) //boundary closure
                    {
                        xArray[vertIndex] = xd1;
                        yArray[vertIndex] = yd1;
                        vertIndex++;
                    }
                    if (propIndex >= 0)
                    {
                        gds2file.putBndDbl(layer, boxType, xArray, yArray, vertIndex, propIndex, propNumArray, propValueArray, uUnits );
                    }
                    else
                    {
                        gds2file.putBndDbl(layer, boxType, xArray, yArray, vertIndex, uUnits);
                    }
                }
                else if (
                    (! strncmp("xy(",field,3)) &&
                    (! sawXy)
                )
                {
                    inXy = TRUE;
                    sawXy = TRUE;
                }
                else if (
                          (! strncmp("pr{",field,3)) ||
                          (! strncmp("pr(",field,3))
                )
                {
                    propIndex++;
                    i = atoi(&field[3]);
                    propNumArray[propIndex] = i; //PROPATTR
                    propValueArray[propIndex][0] = '\0';
                    pos = getUnescapedPropString(line,field,++pos);
                    strcpy(propValueArray[propIndex],field); //PROPVALUE
                }
                else if ((field[0] != ')') && (field[0] != '}'))
                {
                    fprintf(stderr,"%sBad box '%s'. Field='%s' LineNum=%d\n",ERRORSTRING,line,field,lineCnt);
                    exit(5);
                }
            }
        }
        else if (line[0] == 'n') //n{12 xy(90.825 28.9 96.725 28.9 96.725 31.65 90.825 31.65)} --- node
        {
            sawXy = FALSE;
            layer = 0;
            nodeType = 0;
            for(done=FALSE, i=0; (! done) && (i <= stringLen); i++)
            {
                if (line[i] == '{' )
                {
                    done = TRUE;
                }
                line[i] = ' ';
            }
            if (done)
            {
                sStdWhiteSpace(line,0);
            }
            else
            {
                fprintf(stderr,"%sBad node '%s' lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            pos = 0;
            pos = getField(line, field, pos); // layer
            if (sIsPosInt(field))
            {
                layer = atoi(field);
            }
            else
            {
                fprintf(stderr,"%sBad node '%s' bad layer. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            nodeType = 0;
            vertIndex = 0;
            inXy = FALSE;
            while ((pos = getField(line, field, pos)) > 0)
            {
                if ((field[0] == ')') || (field[0] == '}'))
                {
                    inXy = FALSE;
                }
                if (inXy)
                {
                    if (vertIndex < 1)
                    {
                        xd1 = atof(field);
                        xd2 = xd1;
                    }
                    else
                    {
                        xd2 = atof(field);
                    }
                    if ((pos = getField(line, field, pos)) > 0)
                    {
                        if (vertIndex < 1)
                        {
                            yd1 = atof(field);
                            yd2 = yd1;
                        }
                        else
                        {
                            yd2 = atof(field);
                        }
                    }
                    else
                    {
                        fprintf(stderr,"%sBad node '%s' bad xy() section. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                        exit(5);
                    }
                    xArray[vertIndex] = xd2;
                    yArray[vertIndex] = yd2;
                    vertIndex++;
                }
                else if (! strncmp("dt",field,2))
                {
                    strcpy(field2,&field[2]);
                    nodeType = atoi(field2);
                }
                else if (field[0] == '}') // end of node
                {
                    if ((xArray[0] != xArray[(vertIndex - 1)]) || (yArray[0] != yArray[(vertIndex - 1)])) //boundary closure
                    {
                        xArray[vertIndex] = xd1;
                        yArray[vertIndex] = yd1;
                        vertIndex++;
                    }
                    if (propIndex >= 0)
                    {
                        gds2file.putBndDbl(layer, nodeType, xArray, yArray, vertIndex, propIndex, propNumArray, propValueArray, uUnits );
                    }
                    else
                    {
                        gds2file.putBndDbl(layer, nodeType, xArray, yArray, vertIndex, uUnits);
                    }
                }
                else if (
                    (! strncmp("xy(",field,3)) &&
                    (! sawXy)
                )
                {
                    inXy = TRUE;
                    sawXy = TRUE;
                }
                else if (
                          (! strncmp("pr{",field,3)) ||
                          (! strncmp("pr(",field,3))
                )
                {
                    propIndex++;
                    i = atoi(&field[3]);
                    propNumArray[propIndex] = i; //PROPATTR
                    propValueArray[propIndex][0] = '\0';
                    pos = getUnescapedPropString(line,field,++pos);
                    strcpy(propValueArray[propIndex],field); //PROPVALUE
                }
                else if ((field[0] != ')') && (field[0] != '}'))
                {
                    fprintf(stderr,"%sBad node '%s'. Field='%s' LineNum=%d\n",ERRORSTRING,line,field,lineCnt);
                    exit(5);
                }
            }
        }
        else if ((strlen(line)==1) && (line[0] == '}'))
        {
            if (inStr)
            {
                inStr = FALSE;
                gds2file.endStr();
            }
            else
            {
                inLib = FALSE;
                gds2file.endLib();
            }
        }
        else if (! strncmp("lib ",line,4)) //lib 'cc1804' 0.001 1e-09
        {
            xd2 = atof(field);
            pos = 3; //"lib '"
            if ((pos = getField(line, field, pos)) > 0)
            {
                field[LENGTHLSTRING] = '\0'; //safety
                strcpy(libName,&field[1]);
                libName[strlen(libName) - 1] = '\0';
            }
            else
            {
                fprintf(stderr,"%sBad text '%s' bad lib line. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            if ((pos = getField(line, field, pos)) > 0)
            {
                uUnits = atof(field);
            }
            else
            {
                fprintf(stderr,"%sBad text '%s' bad lib line. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            if ((pos = getField(line, field, pos)) > 0)
            {
                dbUnits = atof(field);
            }
            else
            {
                fprintf(stderr,"%sBad text '%s' bad lib line. lineNum=%d\n",ERRORSTRING,line,lineCnt);
                exit(5);
            }
            gds2file.initLib(libName,uUnits,dbUnits,mYear,mMon,mDay,mHour,mMin,mSec,aYear,aMon,aDay,aHour,aMin,aSec,gds2Version);
        }
        else if (line[0] != '#')
        {
            fprintf(stderr,"%sBad/unknown line \"%s\"\n",ERRORSTRING,line);
            //exit(5);
        }
    }
    //////////////////////////////////////////////////////////////////////

    if ((strlen(line)>1) || (strncmp("}",line,1)))
    {
        fprintf(stderr,"%sBad last line of %s - should have been '}'\n",ERRORSTRING,inputFile);
        exit(6);
    }
    if (inLib)
    {
        fprintf(stderr,"%sDid not find '}' end of lib closure\n",ERRORSTRING);
        exit(7);
    }
    gds2file.clstrm();
    return(0);
}  // end of main

/***** GET FIELD *************************************************************/
/* Returns position of string passed in which ends the field starting at the */
/* position passed in with the integer variable. The field string is updated */
/* to the new field.                                                         */
/* ie.                                                                       */
/*      for (n = numfields(string), i=0; n>0; n--) {                         */
/*          i = getField(string, name, i);                                   */
/*          printf("%s\n", name);                                            */
/*      }                                                                    */
/* or                                                                        */
/*      i=0                                                                  */
/*      while ((i=getField(string, name, i)) != ERR) printf("%s\n", name);   */
/*                                                         -KVS              */
/*****************************************************************************/
int getField(char* string, char* field_string, int position)
{
    int     i, j, k;

    for(i=position; (string[i] != '\0') &&   (isspace(string[i])); i++) ;  /* skip white space */
    for(j=i, k=0;  ((string[j] != '\0') && (!(isspace(string[j])))); j++, k++) field_string[k] = string[j];
    field_string[k] = '\0';

    if(string[i] == '\0') return(ERR);
    else return(j);
}

/*****************************************************************************/
int sIsPosInt(char* string)
{
    int i,len;
    len = strlen(string);
    for(i=0;i<len;i++)
    {
        if (! isdigit(string[i])) return(0);
    }
    return(1);
}

/***** rmLeadingAndTrailingSpace ********************************************/
/* Removes all occurances of white space at beginning of a string    -kvs   */
/****************************************************************************/
void rmLeadingAndTrailingSpace(char* parent)
{
    int childIndex = 0,
        parentIndex = 0,
        lastChar = 0,
        length = 0;
    char c;
    static char child[LINELENGTH];

    if (parent[0] != '#')
    {
        length = strlen(parent);
        lastChar = length;
        for(childIndex=0; childIndex<=lastChar; childIndex++)
        {
            child[childIndex] = '\0';
        }
        childIndex = 0;
        //find 1st non white space character position
        while (isspace(parent[parentIndex])) parentIndex++;
        //find last non white space character position
        while (isspace(parent[lastChar]) || (parent[lastChar] == '\0')) lastChar--;

        for(childIndex=0; parentIndex<=lastChar; parentIndex++)
        {
            c = parent[parentIndex];
            child[childIndex++] = c;
        }
        child[childIndex] = '\0';
        strcpy(parent,child);
    }
}

/***** sStdWhiteSpace *******************************************************/
/* Removes all occurances of white space fore and aft in a string    -kvs   */
/* and standardize to one space between fields                              */
/****************************************************************************/
void sStdWhiteSpace(char* parent, int isText)
{
    int childIndex = 0,
        parentIndex = 0,
        lastChar = 0,
        inSpace = 0,
        length = 0,
        okToModify = 1;
    char c;
    static char child[LINELENGTH];

    if (parent[0] != '#')
    {
        length = strlen(parent);
        lastChar = length;
        for(childIndex=0; childIndex<=lastChar; childIndex++)
        {
            child[childIndex] = '\0';
        }
        childIndex = 0;
        //find 1st non white space character position
        while (isspace(parent[parentIndex])) parentIndex++;
        //find last non white space character position
        while (isspace(parent[lastChar]) || (parent[lastChar] == '\0')) lastChar--;

        for(childIndex=0; parentIndex<=lastChar; parentIndex++)
        {
            c = parent[parentIndex];
            if (isspace(c)) inSpace++;
            else            inSpace = 0;
            if (! inSpace)
            {
                if (isText && (c == '\''))
                {
                    if (okToModify) // the first ' of a text string
                    {
                        okToModify = 0;
                    }
                    else // the last ' of a text string
                    {
                        okToModify = 1;
                    }
                }
                // add space before ) or } if not there already
                if (
                    okToModify &&
                    ((c == ')') || (c == '}')) &&
                    (parentIndex > 0) &&
                    (parent[(parentIndex - 1)] != ' ')
                )
                {
                    child[childIndex++] = ' ';
                }

                child[childIndex] = c;

                // add space after '(' if not there already
                if ((c == '(') && (parent[(parentIndex + 1)] != ' '))
                {
                    child[++childIndex] = ' ';
                }
                childIndex++;
            }
            else if (inSpace == 1) //first space only
            {
                child[childIndex++] = ' ';
                inSpace++;
            }
        }
        child[childIndex] = '\0';
        strcpy(parent,child);
    }
}

/***** getUnescapedPropString *****************************************************/
/* remove enclosing 's and escape characters                         -kvs   */
/****************************************************************************/
int getUnescapedPropString(char* parent, char* child, int position)
{
    int i, j;
    Boolean done = FALSE;

    for(i=(position + 1), j=0; (! done); i++)
    {
        if (parent[i] == '\0')
        {
            done = TRUE;
        }
        else if ((parent[i] == '\'') && (parent[(i - 1)] != '\\') && (parent[(i + 1)] == ' ') && (parent[(i + 2)] == '}'))
        {
            done = TRUE;
            i+=2;
        }
        else
        {
            if (parent[i] == '\\') i++;
            child[j++] = parent[i];
        }
    }
    child[j] = '\0';
    return(i);
}

/***** getUnescapedTextString ***********************************************/
/* remove enclosing 's and escape characters                         -kvs   */
/****************************************************************************/
int getUnescapedTextString(char* parent, char* child, int position)
{
    int i, j;
    Boolean done = FALSE;

    for(i=(position + 1), j=0; (! done); i++)
    {
        if (parent[i] == '\0')
        {
            done = TRUE;
        }
        else if ((parent[i] == '\'') &&
                (parent[(i - 1)] != '\\') &&
                ((parent[(i + 1)] == ' ') || (parent[(i + 1)] == '}'))
        )
        {
            done = TRUE;
            i++;
        }
        else
        {
            if (parent[i] == '\\') i++;
            child[j++] = parent[i];
        }
    }
    child[j] = '\0';
    return(i);
}

/***** sUnescapedString *****************************************************/
/* remove enclosing 's and escape characters                         -kvs   */
/****************************************************************************/
char* sUnescapedString(char* parent, char* child)
{
    int i, j, k;

    k = strlen(parent);
    for(i=1, j=0; i<k; i++)
    {
        if (parent[i] == '\\') i++;
        child[j] = parent[i];
    }
    child[j] = '\0';
    return(child);
}

// ****************************************************************************
// * print_help()
// ****************************************************************************
void print_help()
{
    cout << endl;
    cout << "gdt2gds provides you with a tool to convert a textual gdt file to a binary gds2 file " << endl;
    cout << endl;
    cout << "Usage:" << endl;
    cout << "    gdt2gds [inputFile] [outputFile] [OPTIONS]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << " Note: options are case-insensitive and can be shortened as long they remain unique." << endl;
    cout << endl;
    cout << "  -help" << endl;
    cout << "    print this and exit" << endl;
    cout << endl;
    cout << "  -inFilename [inputFile]" << endl;
    cout << "    default is STDIN" << endl;
    cout << "    if inputFile == '-' read from STDIN" << endl;
    cout << endl;
    cout << "  -outFilename [outputFile]  : print this" << endl;
    cout << "    if outputFile == '-' print to STDOUT (screen)" << endl;
    cout << endl;
    cout << "  -precision [integer]" << endl;
    cout << "    Currenly a NOOP. Only included to be symmetric with gds2gdt" << endl;
    cout << endl;
    cout << "  -version" << endl;
    cout << "    print version of program and quit" << endl;
    cout << endl;
    cout << "Examples:" << endl;
    cout << "  gdt2gds test.gdt test.gds" << endl;
    cout << endl;
    cout << "  cat test.gdt | gdt2gds | gzip > test.gds.gz" << endl;
    cout << endl;
}

/*

__END__
use pod2html to make web page help
use pod2man to make man page

=pod

=head1 NAME

gdt2gds - tool to convert textual gdt format to binary gds2

=head1 SYNOPSIS

gdt2gds [options] [inputFile] [outputFile]

=head1 OPTIONS

Note: options are case-insensitive and can be shortened as long they remain unique.

  -help
    print this and exit

  -inFilename <fileName>
    default is STDIN
    also filename of '-' means STDIN

  -outFilename <fileName>
    default is STDOUT
    also filename of '-' means STDOUT

  -precision <integer>
    Currenly a NOOP. Only included to be symmetric with gds2gdt

  -version
    print version of program and quit

=head1 EXAMPLES

gdt2gds test.gdt test.gds

cat test.gdt | gdt2gds | gzip > test.gds.gz

=cut

*/

