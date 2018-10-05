/*
 ******************************************************************************
 * Copyright 2004-2014 by Ken Schumack
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
// @(#) $Id: gds2gdt.C 95 2014-12-08 17:16:36Z schumack $
//#include <stdio.h>
//#include <iostream.h>
//#include <iostream.h>
//#include <sstream>
//#include <stdlib.h>
//#include <ostream.h>
//#include <math.h>
//#include <kvstypes.h>
#include <kvsstring_c.h>
#include <gdsStream.h>
#include <iostream>
#include <gcc4.h>
using namespace std;
#include <sourceForge.h>

extern void print_help();

//////////////////////////////////////////////////////////////////////////////
// M A I N
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    int     i,
            rectyp,
            length,
            debug=0,
            extn,
            width,
            precision = 0,
            haveInputName=0,
            haveOutputName=0,
            inInFileName=0,
            inOutFileName=0,
            inPrecision=0,
            inProperty=0,
            inBoundary=0;
    short   Ref, columns, rows, layer, dataType, textType, propAttr,
            year, month, day, hour, min, sec,
            version_num;
    double  userUnits=0.001,           //1e-3 is common default
            dataBaseUnits=0.000000001, //1e-9 is common default
            Ang=0.0,
            epsilon=(userUnits/2000.0), //use to "fix" floating point error
            Mag=0.0;
    char    oneChar,
            precisionString[128] = "0.001",
            bigString[204800] = "",
            string512[513] = ""; // 512 + room for '\0'
    stringL inputFile, outputFile, tmpFile, strname, tmpString1, tmpString2;
    stringS xstring, ystring, extnstring, wstring;
    FILE    *fpOut;
    long    xcoord=0,
            ycoord=0;

    // next line uses '@ ( # )' so the unix/linux 'what' command works on the executable
    cerr << "# ** gds2gdt @(#) Version " << sourceForgeVersion << " $Id: gds2gdt.C 95 2014-12-08 17:16:36Z schumack $ ** #" << endl; //using subversion props

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
                haveOutputName = 1;
                fpOut = stdout;
            }
            else if (match_string("-debug", argv[i], 'm'))
            {
                debug = 1;
            }
            else if (match_string("-help", argv[i], 'm'))
            {
                print_help();
                exit(0);
            }
            else if (match_string("-infilename", argv[i], 'm'))
            {
                inInFileName = 1;
            }
            else if (match_string("-precision", argv[i], 'm'))
            {
                inPrecision = 1;
            }
            else if (match_string("-outfilename", argv[i], 'm'))
            {
                inOutFileName = 1;
            }
            else if (match_string("-version", argv[i], 'm') || match_string("-swversion", argv[i], 'm'))
            {
                exit(0);   // version already printed */
            }
            else
            {
                cerr << "ERROR **** unknown option " << argv[i] << endl;
                exit(2);
            }
        }
        else
        {
            if (inInFileName)
            {
                haveInputName = 1;
                inInFileName = 0;
                strcpy(inputFile, argv[i]);
            }
            else if (inOutFileName)
            {
                haveOutputName = 1;
                inOutFileName = 0;
                strcpy(outputFile, argv[i]);
            }
            else if (inPrecision)
            {
                inPrecision = 0;
                precision = atoi(argv[i]);
            }
            else if (! haveInputName)
            {
                haveInputName = 1;
                strcpy(inputFile, argv[i]);
            }
            else if (! haveOutputName)
            {
                haveOutputName = 1;
                strcpy(outputFile, argv[i]);
            }
        }
    }

    GDSFILE gds2file(inputFile, 0);
    if (! haveOutputName)
    {
        fpOut = stdout;
    }
    else if (fpOut != stdout)
    {
        if (strEqual(gds2file.fileName(),outputFile))
        {
            sprintf(outputFile,"%s.gdt",outputFile);
            cerr << "Warning ** you gave the input filename for the output file.  Will use " << outputFile << "instead" << endl;

        }
        if((fpOut = fopen(outputFile,"w")) == NULL)
        {
            cerr << "ERROR **** unable to create file " << outputFile << endl;
            exit(1);
        }
    }
    #ifdef DEBUG
    setbuf(fpOut, NULL);
    #endif
    // ************* end of command line and menu stuff ***************************

    gds2file.rdstrm();  // header
    version_num = gds2file.getI16();
    fprintf(fpOut, "gds2{%d\n",version_num);
    gds2file.rdstrm();  // bgnlib
        year  = gds2file.getI16(0);
        if (year < 999) year += 1900;
        month = gds2file.getI16(2);
        day   = gds2file.getI16(4);
        hour  = gds2file.getI16(6);
        min   = gds2file.getI16(8);
        sec   = gds2file.getI16(10);
    fprintf(fpOut, "m=%d-%02d-%02d %02d:%02d:%02d",year,month,day,hour,min,sec);
        year  = gds2file.getI16(12);
        if (year < 999) year += 1900;
        month = gds2file.getI16(14);
        day   = gds2file.getI16(16);
        hour  = gds2file.getI16(18);
        min   = gds2file.getI16(20);
        sec   = gds2file.getI16(22);
    fprintf(fpOut, " a=%d-%02d-%02d %02d:%02d:%02d\n",year,month,day,hour,min,sec);
    while (! gds2file.eof())
    {
        gds2file.rdstrm();
        rectyp = gds2file.rectyp();
        if (debug) cerr << "DEBUG line:" << __LINE__ << "rectyp=" << rectyp << endl;
        if ((rectyp < 0) || (rectyp > 59))
        {
            fprintf(stderr, "ERROR: invalid record type:%d found in gds2 file. Note: May have just overflowed due to a super long record\n", rectyp);
            exit(1);
        }
        if (rectyp == LIBNAME)
        {
            gds2file.libName((char*) gds2file.record());
            fprintf(fpOut, "lib '%s'",gds2file.record());
            strcpy(strname, gds2file.record());
        }
        else if (rectyp == BGNSTR)
        {
            year  = gds2file.getI16(0);
                if (year < 999) year += 1900;
                month = gds2file.getI16(2);
                day   = gds2file.getI16(4);
                hour  = gds2file.getI16(6);
                min   = gds2file.getI16(8);
                sec   = gds2file.getI16(10);
            fprintf(fpOut, "cell{c=%d-%02d-%02d %02d:%02d:%02d",year,month,day,hour,min,sec);
                year  = gds2file.getI16(12);
                if (year < 999) year += 1900;
                month = gds2file.getI16(14);
                day   = gds2file.getI16(16);
                hour  = gds2file.getI16(18);
                min   = gds2file.getI16(20);
                sec   = gds2file.getI16(22);
            fprintf(fpOut, " m=%d-%02d-%02d %02d:%02d:%02d",year,month,day,hour,min,sec);
        }
        else if (rectyp == UNITS)
        {
            userUnits     = gds2file.getDbl();   // Calma default is 1.0e-3
            epsilon = (userUnits/2000.0);  //reset
            dataBaseUnits = gds2file.getDbl(8);  // Calma default is 1.0e-9
            if (precision <= 0) // then not set to positive integer on command line
            {
                sprintf(tmpString1,"%0.8f",userUnits);
                sRemoveTrailingZeros(tmpString1,precisionString);
                precision = strlen(precisionString) - 2; // "0.001" -> 3, "0.05" -> 2
                if (debug) cerr << "DEBUG line:" << __LINE__ << " precisionString=" << precisionString << " precision=" << precision << endl;
            }
            if (debug) cerr << "DEBUG line:" << __LINE__ << " userUnits=" << userUnits << " dataBaseUnits=" << dataBaseUnits << endl;
            fprintf(fpOut, " %g %g\n",userUnits,dataBaseUnits);
            fprintf(fpOut, "# lines above need to stay as is (may be read by other tools)\n");
            fprintf(fpOut, "# File created by gds2gdt  http://sourceforge.net/projects/gds2/  version:%s \n", sourceForgeVersion);
            fprintf(fpOut, "# Key: <required> [optional]\n");
            fprintf(fpOut, "# File format:\n");
            fprintf(fpOut, "# gds2{<ver>\n");
            fprintf(fpOut, "# m=<modificationTimeStamp> a=<accessTimeStamp>\n");
            fprintf(fpOut, "# lib '<libName>' <userUnits> <dataUnits>\n");
            fprintf(fpOut, "# <cellDefine>\n");
            fprintf(fpOut, "# }\n");
            fprintf(fpOut, "# - - - - -\n");
            fprintf(fpOut, "# cellDefine is one of more of:\n");
            fprintf(fpOut, "# cell {c=<creationTimeStamp> m=<modificationTimeStamp> '<cellName>'\n");
            fprintf(fpOut, "# <cellStuff>*\n");
            fprintf(fpOut, "# }\n");
            fprintf(fpOut, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
            fprintf(fpOut, "## <cellStuff>\n");
            fprintf(fpOut, "# cellStuff is one or more of:\n");
            fprintf(fpOut, "# boundary:\n");
            fprintf(fpOut, "# b{<layer> [dt<dataType>] xy(<xyList>) [property]*}\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# path:\n");
            fprintf(fpOut, "# p{<layer> [dt<dataType>] [pt<pathType>] [w<real>] [bx<real>] [ex<real>] xy(<xyList>) [property]*}\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# text:\n");
            fprintf(fpOut, "# t{<layer> [tt<textType>] [f<fontType>] [<textJust>] [pt<pathType>] [fx] [w<real>] [m<magification>] [a<angle>] xy(<xyList>) <'text'> [property]*}\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# sref:\n");
            fprintf(fpOut, "# s{<'cellName'> [fx] [a<angle>] xy(<xyList>) [property]*}\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# aref:\n");
            fprintf(fpOut, "# a{<'cellName'> [fx] [a<angle>] cr(<columns> <rows>) xy(<xyList>) [property]*}\n");
            fprintf(fpOut, "#   aref xyList: 1st coord: origin, 2nd coord: X of col * xSpacing + origin, 3rd coord: Y of row * ySpacing + origin\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# box:\n");
            fprintf(fpOut, "# x{<layer> [xt<boxType>] xy(<xyList>) [property]*}\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# node:\n");
            fprintf(fpOut, "# n{<layer> [nt<boxType>] xy(<xyList>) [property]*}\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
            fprintf(fpOut, "# property : pr{<propAttr> <'propValue'>}\n");
            fprintf(fpOut, "# <text> : ASCII String\n");
            fprintf(fpOut, "# <textJust> : two letter combination of bmt (bottom,middle,top) and rcl (right,center,left) e.g. bl (default is tl)\n");
            fprintf(fpOut, "# <propAttr> : a 2 byte (small) integer\n");
            fprintf(fpOut, "# <propValue> : ASCII String\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# x==origin\n");
            fprintf(fpOut, "#     _____\n");
            fprintf(fpOut, "#    |    a0      x_________        a180    | x             a270\n");
            fprintf(fpOut, "#    |--                |   |               |        |\n");
            fprintf(fpOut, "#    |                      |             --|        |___|______\n");
            fprintf(fpOut, "#  x |                  a90            _____|                   x\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "#                   fx a90             _____          __________x\n");
            fprintf(fpOut, "#  x |    fx                |               |        |   |\n");
            fprintf(fpOut, "#    |            ______|___|             --|        |\n");
            fprintf(fpOut, "#    |--          x              fx a180    |            fx a270\n");
            fprintf(fpOut, "#    |_____                                 | x\n");
            fprintf(fpOut, "#\n");
            fprintf(fpOut, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
            fprintf(fpOut, "# # as first character on a line is a comment\n");
        }
        else if (rectyp == STRNAME)
        {
            strcpy(strname, gds2file.record());
            fprintf(fpOut, " '%s'\n",strname);
            //cell{c=1998-08-17 14:31:10 m=1998-08-17 14:33:47 'CC1804_R1'  //...}
        }
        else if (rectyp == BOUNDARY)
        {
            inBoundary = 1;
            fprintf(fpOut, "b{");
        }
        else if (rectyp == PATH)
        {
            fprintf(fpOut, "p{");
        }
        else if ((rectyp == ENDEL) || (rectyp == ENDSTR) || (rectyp == ENDLIB))
        {
            inBoundary = 0;
            fprintf(fpOut, "}\n");
        }
        else if (rectyp == COLROW)
        {
            columns = gds2file.getI16();
            rows = gds2file.getI16(2);
            fprintf(fpOut," cr(%d %d)",columns,rows);
        }
        else if (rectyp == PATHTYPE)
        {
            Ref = gds2file.getI16();
            if (Ref != 0)
            {
                fprintf(fpOut," pt%d",Ref);
            }
        }
        else if (rectyp == STRANS)
        {
            Ref = gds2file.getI16();
            if (Ref & 0x8000)
            {
                fprintf(fpOut," fx");
            }
        }
        else if (rectyp == PRESENTATION)
        {
            Ref = gds2file.getI16();
            //font number
            if ((Ref & 0x30) == 0x30)
            {
                fprintf(fpOut," f3");
            }
            else if ((Ref & 0x20) == 0x20)
            {
                fprintf(fpOut," f2");
            }
            else if ((Ref & 0x10) == 0x10)
            {
                fprintf(fpOut," f1");
            }
            //else sprintf(tmpString1," f0");

            tmpString1[0] = '\0';
            // bottom middle top
            if ((Ref & 0x8) == 0x8)
            {
                sprintf(tmpString1,"b");
            }
            else if ((Ref & 0x4) == 0x4)
            {
                sprintf(tmpString1,"m");
            }
            else
            {
                sprintf(tmpString1,"t");
            }
            // right center left
            if ((Ref & 0x2) == 0x2)
            {
                strcat(tmpString1,"r");
            }
            else if ((Ref & 0x1) == 0x1)
            {
                strcat(tmpString1,"c");
            }
            else
            {
                strcat(tmpString1,"l");
            }
            if (strncmp("tl",tmpString1,2)) fprintf(fpOut," %s",tmpString1); //save space - default
        }
        else if (rectyp == TEXT)
        {
            fprintf(fpOut, "t{");
        }
        else if (rectyp == SREF)
        {
            fprintf(fpOut, "s{");
        }
        else if (rectyp == AREF)
        {
            fprintf(fpOut, "a{");
        }
        else if (rectyp == SNAME)
        {
            fprintf(fpOut, "'%s'", gds2file.record());
        }
        else if ((rectyp == STRING) || (rectyp == PROPVALUE))
        {
            fprintf(fpOut, " '");
            strcpy(bigString,gds2file.record());
            strncpy(string512,bigString,512);
            for (i=0; i<strlen(string512); i++)
            {
                oneChar = string512[i];
                if (i==512)
                {
                    string512[512] = '\0';
                    cerr << "WARNING **** STRING may be longer than 512 characters - truncated to" << string512 << argv[i] << endl;
                    break;
                }
                if (oneChar == '\r')
                {
                    //fprintf(fpOut, "%c", '\¤'); //do nothing
                }
                else if (oneChar == '\n')
                {
                    fprintf(fpOut, "%c", '\r'); // will print out as ^M char
                }
                else if (oneChar == '\'')
                {
                    fprintf(fpOut, "%c%c", '\\','\''); // escape
                }
                else if (oneChar == '\\')
                {
                    fprintf(fpOut, "%c%c", '\\','\\'); // escape
                }
                else
                {
                    fprintf(fpOut, "%c", oneChar);
                }
            }
            fprintf(fpOut, "'");
            if (inProperty) fprintf(fpOut, "}");
            inProperty=0;
        }
        else if (rectyp == XY)
        {
            length = gds2file.length();
            if (debug) cerr << "DEBUG line:" << __LINE__ << " XY length=" << length << endl;
            if ((length > 65524) || (length < 0)) // 0xffff - 4 for header / 8 = 8191 full 8 byte pairs -- 0=>overflow
            {
                fprintf(stderr, "ERROR: unsupport XY length (%d - more than 8191 points) found in gds2 file.\n",length);
                exit(1);
            }
            if (inBoundary) length -= 8; //remove closure -- last 2 points
            fprintf(fpOut, " xy(");
            int skip=0;
            for(i=0; i<length; i+=8)
            {
                if (i == MAXPAIRSXY)
                {
                    cerr << "ERROR: element has > " << MAXPAIRSXY << " coordinates. GDT output lines xy() list is truncated." << endl;
                    skip=1;
                }
                xcoord = gds2file.getI32(i);
                ycoord = gds2file.getI32(i+4);
                if (xcoord < 0.0)
                {
                    sprintf(xstring,"%0.*f",precision,(xcoord * userUnits) - epsilon);
                }
                else
                {
                    sprintf(xstring,"%0.*f",precision,(xcoord * userUnits) + epsilon);
                }
                if (ycoord < 0.0)
                {
                    sprintf(ystring,"%0.*f",precision,(ycoord * userUnits) - epsilon);
                }
                else
                {
                    sprintf(ystring,"%0.*f",precision,(ycoord * userUnits) + epsilon);
                }
                if (! skip) fprintf(fpOut,"%s%s %s",(i?" ":""),sRemoveTrailingZeros(xstring,tmpString1),sRemoveTrailingZeros(ystring,tmpString2));
            }
            if (debug) cerr << "DEBUG line:" << __LINE__ << " LAST XY i=" << i << " xstring=" << xstring << " ystring=" << ystring << endl;
            inBoundary = 0; //done
            fprintf(fpOut, ")");
        }
        else if (rectyp == LAYER)
        {
            layer = gds2file.getI16();
            fprintf(fpOut,"%d",layer);
        }
        else if (rectyp == WIDTH)
        {
            width = gds2file.getI32();
            if (width != 0)
            {
                sprintf(wstring,"%0.5f",((double)width * userUnits) + epsilon);
                fprintf(fpOut," w%s",sRemoveTrailingZeros(wstring,tmpString1));
            }
        }
        else if (rectyp == DATATYPE)
        {
            dataType = gds2file.getI16();
            if (dataType != 0) fprintf(fpOut," dt%d",dataType);
        }
        else if (rectyp == TEXTTYPE)
        {
            textType = gds2file.getI16();
            if (textType != 0) fprintf(fpOut," tt%d",textType);
        }
        else if (rectyp == ANGLE)
        {
            Ang = gds2file.getDbl();
            if (Ang != 0.0)
            {
                sprintf(xstring,"%0.5f",Ang);
                fprintf(fpOut," a%s",sRemoveTrailingZeros(xstring,tmpString1));
            }
        } // end ANGLE
        else if (rectyp == MAG)
        {
            Mag = gds2file.getDbl();
            if (Mag != 1.0)
            {
                sprintf(xstring,"%0.5f",Mag);
                fprintf(fpOut," m%s",sRemoveTrailingZeros(xstring,tmpString1));
            }
        } // end MAG
        else if (rectyp == BGNEXTN)
        {
            extn = gds2file.getI32();
            sprintf(extnstring,"%0.5f",((double)extn * userUnits) + epsilon);
            fprintf(fpOut," bx%s",sRemoveTrailingZeros(extnstring,tmpString1));
        }
        else if (rectyp == ENDEXTN)
        {
            extn = gds2file.getI32();
            sprintf(extnstring,"%0.5f",((double)extn * userUnits) + epsilon);
            fprintf(fpOut," ex%s",sRemoveTrailingZeros(extnstring,tmpString1));
        }
        else if (rectyp == PROPATTR)
        {
            propAttr = gds2file.getI16();
            fprintf(fpOut," pr{%d",propAttr);
            inProperty=1;
        }
        else if (rectyp == NODE)
        {
            fprintf(fpOut, "n{");
        }
        else if (rectyp == NODETYPE)
        {
            Ref = gds2file.getI16();
            if (Ref != 0)
            {
                fprintf(fpOut," nt%d",Ref);
            }
        }
        else if (rectyp == BOX)
        {
            fprintf(fpOut, "n{");
        }
        else if (rectyp == BOXTYPE)
        {
            Ref = gds2file.getI16();
            if (Ref != 0)
            {
                fprintf(fpOut," xt%d",Ref);
            }
        }
    }
    return(0);
}  // end of main

// ****************************************************************************
// * print_help()
// ****************************************************************************
void print_help()
{
    cout << endl;
    cout << "gds2gdt provides you with a tool to convert a binary gds2 files to a textual" << endl;
    cout << "gdt file." << endl;
    cout << endl;
    cout << "Usage:" << endl;
    cout << "    gds2gdt [inputFile] [outputFile] [OPTIONS]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << " Note: options are case-insensitive and can be shortened as long they remain unique." << endl;
    cout << endl;
    cout << "  -help" << endl;
    cout << "    Print this and exit" << endl;
    cout << endl;
    cout << "  -inFilename [inputFile]" << endl;
    cout << "    Default is STDIN" << endl;
    cout << "    If inputFile == '-' read from STDIN" << endl;
    cout << endl;
    cout << "  -outFilename [outputFile]  : print this" << endl;
    cout << "    If outputFile == '-' print to STDOUT (screen)" << endl;
    cout << endl;
    cout << "  -precision [integer]" << endl;
    cout << "    Default is taken from GDS2 files userUnits. I.E. 0.001 gives 3" << endl;
    cout << "    Print xy floats rounded this many places right of decimal" << endl;
    cout << endl;
    cout << "  -version" << endl;
    cout << "    Print version of program and quit" << endl;
    cout << endl;
    cout << "Examples:" << endl;
    cout << "  gds2gdt test.gds test.gdt" << endl;
    cout << endl;
    cout << "  zcat test.gds.gz | gds2gdt | grep ^cell" << endl;
    cout << endl;
}

/*

__END__
use pod2html to make web page help
use pod2man to make man page

=pod

=head1 NAME

gds2gdt - tool to convert binary gds2 to textual gdt format

=head1 SYNOPSIS

gds2gdt [options] [inputFile] [outputFile]

=head1 OPTIONS

Note: options are case-insensitive and can be shortened as long they remain unique.

  -help
    Print this and exit

  -inFilename <fileName>
    Default is STDIN
    Also filename of '-' means STDIN

  -outFilename <fileName>
    Default is STDOUT
    Also filename of '-' means STDOUT

  -precision <integer>
    Default is taken from GDS2 files userUnits. I.E. 0.001 gives 3
    Print xy floats rounded this many places right of decimal

  -version
    Print version of program and quit

=head1 EXAMPLES

gds2gdt test.gds test.gdt

zcat test.gds.gz | gds2gdt | grep ^cell

=cut
*/

