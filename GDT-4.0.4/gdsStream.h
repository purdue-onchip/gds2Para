// Copyright 1995-2014 by Ken Schumack (Schumack@cpan.org)
// @(#) $Id: gdsStream.h 95 2014-12-08 17:16:36Z schumack $
#ifndef _gdsStream__
#define _gdsStream__

#include <stdio.h>
#include <malloc.h>
#include <stdio.h>
#include <time.h>
#include <kvstypes.h>

#define NUMGDSLAYERS 1024
#define MAXPAIRSXY 8191     //  int(0xffff / 8) -> 0xffff - 4 byte record header / 4 bytes per int
//G_epsilon: to take care of floating point representation problems
#define G_epsilon 0.00005

static const char* gdsStream_hwhat = "@(#) $Id: gdsStream.h 95 2014-12-08 17:16:36Z schumack $ $Revision: 95 $ $Date: 2014-12-08 11:16:36 -0600 (Mon, 08 Dec 2014) $";

class GDSFILE
{
protected:
    char*   LibName;
    char*   FileName;
    char*   CurrentStrName;
    time_t  Time_val;
    char    Record[204800];     //100x 2048 for super large boundaries ...
    int     Eof;                // are we at the End Of File?
    int     EndOfLib;           // are we at the End Of Lib?
    int     Length;             //current record length
    int     Rectyp;             //current record type
    int     Dattyp;             //current data type
    char    Buffer[204800];
    int     Fd;                 // current file descriptor
    int     Writtn;             // 1 = write, 0 = read
    int     Ptr;                // current pointer in buffer
    short   Glayers[NUMGDSLAYERS];
    short   Tlayers[NUMGDSLAYERS];
    short   LayerDataTypes[NUMGDSLAYERS][NUMGDSLAYERS];
    short   LayerTextTypes[NUMGDSLAYERS][NUMGDSLAYERS];

    void    endEl(); // end of element
public:
    GDSFILE(char* fileName, int readOrWrite); // constructor to open (create if WRITE) stream file (calls opstrm())
    void    opstrm();                 // open stream file
    int     rdstrm();                 // read Record from stream file
    void    wrstrm();                 // write Record to stream file
            // write external record to stream file, Rectyp,Dattyp,Length taken from
            // another stream file (used when changing XY value for instance)
    void    wrstrm(char record[204800], GDSFILE* gdsfile);
            // write external record to stream file, Rectyp,Dattyp,Length taken from args
    void    wrstrm(char record[204800], int Rectyp, int Dattyp, int Length);
    void    cpstrm(GDSFILE* gdsfile); // cp Record from one GDSFILE to another
    void    cpend(GDSFILE* gdsfile);  // cp remaining bytes after ENDLIB from one GDSFILE to another
    void    clstrm();                 // close stream file

    void    initLib(char *library, double dbu_uu, double dbu_m, int myear, int mmon, int mmday, int mhour, int mmin, int msec, int ayear, int amon, int amday, int ahour, int amin, int asec, int version);
    void    initLib(char *library, double dbu_uu, double dbu_m, int myear, int mmon, int mmday, int mhour, int mmin, int msec, int ayear, int amon, int amday, int ahour, int amin, int asec);
    void    initLib(char* library, double dbu_uu, double dbu_m);
    void    initLib(char* library);   //uses dbu_uu==1.0e-3 and dbu_m==1.0e-9
    void    endLib();                 //write end of library record
    char*   libName();                //get stored LibName
    void    libName(char* name);      //store LibName
    void    copyRecord(char* copy);   //copy current record to "copy"

    void    beginStr(char* str_name); //open a new structure
    void    beginStr(char *str_name, int myear, int mmon, int mmday, int mhour, int mmin, int msec, int ayear, int amon, int amday, int ahour, int amin, int asec);

    void    endStr();                 //close a structure

    int     roundInt(int input, int grid);
    // place a structure reference
    void    putSref(char* sname,unsigned short ref,double mag,double angle,double x_coord,double y_coord,int propIndex,int propNumArray[],char propValueArray[][LENGTHLSTRING],double dbu_uu);
    void    putSref(char* sname, unsigned short ref, double mag, double angle, double x_coord, double y_coord, double dbu_uu);
    void    putSref(char* sname, unsigned short ref, double mag, double angle, double x_coord, double y_coord); //default uu==0.001
    // place an array reference
    void    putAref(char*  sname,unsigned short ref,double mag,double angle,short  col,short  row,double x1,double y1,double x2,double y2,double x3,double y3,int propIndex,int propNumArray[],char propValueArray[][LENGTHLSTRING],double dbu_uu);
    void    putAref(char* sname, unsigned short ref, double mag, double angle, short col, short row, double x1, double y1, double x2, double y2, double x3, double y3, double dbu_uu);
    void    putAref(char* sname, unsigned short ref, double mag, double angle, short col, short row, double x1, double y1, double x2, double y2, double x3, double y3); //default uu==0.001

    // place a rectangle
    void    putRt(int layer, int datatyp, double minX, double minY, double maxX, double maxY, double dbu_uu);
    void    putRt(int layer, int datatyp, double minX, double minY, double maxX, double maxY); //default uu==0.001
    // place a text
    void    putText(unsigned short layer, unsigned short ref, double mag, double angle, double x, double y, char* txt, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING], double dbu_uu);
    void    putText(unsigned short layer, unsigned short ref, double mag, double angle, double x, double y, char* txt, double dbu_uu);
    void    putText(unsigned short layer, unsigned short ref, double mag, double angle, double x, double y, char* txt); //default uu==0.001
    void    putText(unsigned short layer, unsigned short textType, unsigned short fontType, char* textJust, unsigned short pathType, double width, unsigned short ref, double mag, double angle, double x, double y, char* txt); //default uu==0.001
    void    putText(unsigned short layer, unsigned short textType, unsigned short fontType, char* textJust, unsigned short pathType, double width, unsigned short ref, double mag, double angle, double x, double y, char* txt, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING], double dbu_uu);
    void    putText(unsigned short layer, unsigned short textType, unsigned short fontType, char* textJust, unsigned short pathType, double width, unsigned short ref, double mag, double angle, double x, double y, char* txt, double dbu_uu);

    // place a boundary using arrays of doubles and props
    int     putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING], double dbu_uu);
    int     putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING]); //default uu==0.001
    // place a boundary using arrays of doubles
    int     putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert, double dbu_uu);
    int     putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert); //default uu==0.001
    // place a boundary using arrays of ints
    int     putBndInt(int layer, int datatyp, int xArray[], int yArray[], int nVert);
    // place a path using arrays of ints
    int     putPathInt(int layer, int datatyp, int width, int xArray[], int yArray[], int nVert);
    // place a path using arrays of doubles
    int     putPathDbl(int layer, int datatyp, int pathtyp, double width, double bgnextn, double endextn, double xArray[], double yArray[], int nVert, double dbu_uu);
    int     putPathDbl(int layer, int dataTyp, int pathTyp, double width, double bgnExtn, double endExtn, double xArray[], double yArray[], int nVert, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING], double uUnits);
    int     putPathDbl(int layer, int datatyp, int pathtyp, double width, double bgnextn, double endextn, double xArray[], double yArray[], int nVert); //default uu==0.001

    int     endoflib(); // are you at END OF LIB?
    int     eof();      // are you at EOF?
    int     rectyp();   // retrive current rectype
    void    rectyp(int);// set current rectype
    int     dattyp();   // retrive current dataType
    void    dattyp(int);// set current dataType
    int     length();   // retrive length of current record
    void    length(int);// set length of current record
    char*   record();   // return current record
    char*   fileName(); // return FileName of stream file.

    void    foundGraphicsLayer(short layerFound);
    void    foundLayerDatatype(short layerFound, short dataTypeFound);
    int     gLayer(short layerNum);    // does this graphics layer exist in stream file?
    int     layerDataType(short layer, short dataType); // does this graphics layer / datatype combo exist?

    void    foundTextLayer(short layerFound);
    void    foundLayerTexttype(short layerFound, short textTypeFound);
    int     tLayer(short layerNum);   // does this text layer exist in stream file?
    int     layerTextType(short layer, short textType); // does this text layer / datatype combo exist?;

    double  getDbl();                       // get double
    double  getDbl(int offset);             // get double at specified offset
    void    putDbl(double dbl, int offset); // put double at specified offset

    int     getI16();                               // get 16 bit integer
    int     getI16(int offset);                     // get 16 bit integer at specified offset
    void    putI16(unsigned short i16, int offset); // put 16 bit integer at specified offset

    int     getI32();                     // get 32 bit integer
    int     getI32(int offset);           // get 32 bit integer at specified offset
    void    putI32(int i32, int offset);  // put 32 bit integer at specified offset

    void    copy(char src_rec[], char dst_rec[], int  num);
    int     iround(long int number, int places);
} //end class GDSFILE

;
//extern foobar(void);
/*************************************************************************************
 * CALMA STREAM FORMAT
 *
The stream format file is composed of variable length records. The mininum
length record is 4 bytes. The 1st 2 btyes of a record contain a count (in 8 bit
bytes) of the total record length.  The 3rd byte of the header is the record
type. The 4th byte describes the type of data contained w/in the record. The
5th through last bytes are data.

If the output file is a mag tape, then the records of the library are written
out in 2048-byte physical blocks. Records may overlap block boundaries.

A null word consists of 2 consecutive zero bytes. Use null words to fill the
space between:
    o the last record of a library and the end of its block
    o the last record of a tape in a mult-reel stream file.

DATA TYPE        VALUE  RECORD
---------        -----  -----------------------
no data present     0   4bytes long
Bit Array           1   2bytes long
2byte Signed Int    2  SMMMMMMM MMMMMMMM  -> S - sign ;  M - magnitude.
                       Twos complement format, with the most significant byte first.
                       I.e.
                        0x0000 = 1
                        0x0020 = 2
                        0x0089 = 137
                        0xffff = -1
                        0xfffe = -2
                        0xff77 = -137

4byte Signed Int    3  SMMMMMMM MMMMMMMM MMMMMMMM MMMMMMMM
8byte Real          5  SEEEEEEE MMMMMMMM MMMMMMMM MMMMMMMM E-expon in excess-64 representation
                       MMMMMMMM MMMMMMMM MMMMMMMM MMMMMMMM
                       Mantissa -> pos fraction >=1/16 and <1 bit 8 = 1/2, 9=1/4 etc...
                        The first bit is the sign (1 = negative), the next 7 bits
                        are the exponent, you have to subtract 64 from this number to get the real
                        value. The next three bytes are the mantissa, divide by 2^24 to get the
                        denominator.
                        value = (mantissa/(2^24)) * (16^(exponent-64))
string              6  odd length strings must be padded w/ null character and byte count++

*****************************************************************************/

/****************************************************************************
 * CALMA STREAM SYNTAX
 *
 <STREAM FORMAT>::=     HEADER BGNLIB [LIBDIRSIZE] [SRFNAME] [LIBSECR]
                        LIBNAME [REFLIBS] [FONTS] [ATTRTABLE] [GENERATIONS]
                        [<FormatType>] UNITS {<structure>}* ENDLIB

 <FormatType>::=        FORMAT | FORMAT {MASK}+ ENDMASKS

 <structure>::=         BGNSTR STRNAME [STRCLASS] {<element>}* ENDSTR

 <element>::=           {<boundary> | <path> | <SREF> | <AREF> | <text> |
                         <node> | <box>} {<property>}* ENDEL

 <boundary>::=          BOUNDARY [ELFLAGS] [PLEX] LAYER DATATYPE XY

 <path>::=              PATH [ELFLAGS] [PLEX] LAYER DATATYPE [PATHTYPE]
                        [WIDTH] XY

 <SREF>::=             SREF [ELFLAGS] [PLEX] SNAME [<strans>] XY

 <AREF>::=             AREF [ELFLAGS] [PLEX] SNAME [<strans>] COLROW XY

 <text>::=             TEXT [ELFLAGS] [PLEX] LAYER <textbody>

 <textbody>::=         TEXTTYPE [PRESENTATION] [PATHTYPE] [WIDTH] [<strans>] XY
                       STRING

 <strans>::=           STRANS [MAG] [ANGLE]

 <node>::=             NODE [ELFLAGS] [PLEX] LAYER NODETYPE XY

 <box>::=              BOX [ELFLAGS] [PLEX] LAYER BOXTYPE XY

 <property>::=         PROPATTR PROPVALUE

****************************************************************************/

// CALMA STREAM RECORD DATATYPES
#define  NO_DATA       0
#define  BIT_ARRAY     1
#define  INTEGER_2     2
#define  INTEGER_4     3
#define  REAL_4        4
#define  REAL_8        5
#define  ACSII_STRING  6

// CALMA STREAM RECORD TYPES
#define  HEADER         0  // 2-byte Signed Integer
#define  BGNLIB         1  // 2-byte Signed Integer
#define  LIBNAME        2  // ASCII String
#define  UNITS          3  // 8-byte Real
#define  ENDLIB         4  // no data present
#define  BGNSTR         5  // 2-byte Signed Integer
#define  STRNAME        6  // ASCII String
#define  ENDSTR         7  // no data present
#define  BOUNDARY       8  // no data p/resent
#define  PATH           9  // no data present
#define  SREF          10  // no data present
#define  AREF          11  // no data present
#define  TEXT          12  // no data present
#define  LAYER         13  // 2-byte Signed Integer
#define  DATATYPE      14  // 2-byte Signed Integer
#define  WIDTH         15  // 4-byte Signed Integer
#define  XY            16  // 2-byte Signed Integer
#define  ENDEL         17  // no data present
#define  SNAME         18  // ASCII String
#define  COLROW        19  // 2-byte Signed Integer
#define  TEXTNODE      20  // no data present
#define  NODE          21  // no data present
#define  TEXTTYPE      22  // 2-byte Signed Integer
#define  PRESENTATION  23  // Bit Array
#define  SPACING       24  // discontinued
#define  STRING        25  // ASCII String
#define  STRANS        26  // Bit Array
#define  MAG           27  // 8-byte Real
#define  ANGLE         28  // 8-byte Real
#define  UINTEGER      29  // UNKNOWN User int, used only in V2.0
#define  USTRING       30  // UNKNOWN User string, used only in V2.0
#define  REFLIBS       31  // ASCII String
#define  FONTS         32  // ASCII String
#define  PATHTYPE      33  // 2-byte Signed Integer
#define  GENERATIONS   34  // 2-byte Signed Integer
#define  ATTRTABLE     35  // ASCII String
#define  STYPTABLE     36  // ASCII String "Unreleased feature"
#define  STRTYPE       37  // 2-byte Signed Integer "Unreleased feature"
#define  EFLAGS        38  // BIT_ARRAY  Flags for template and exterior data.  bits 15 to 0, l to r 0=template, 1=external data, others unused
#define  ELKEY         39  // INTEGER_4  "Unreleased feature"
#define  LINKTYPE      40  // UNKNOWN    "Unreleased feature"
#define  LINKKEYS      41  // UNKNOWN    "Unreleased feature"
#define  NODETYPE      42  // INTEGER_2  Nodetype specification. On GDSII this could be 0 to 63, LTL allows 0 to 255. Of course a 2 byte integer allows up to 65535...
#define  PROPATTR      43  // INTEGER_2  Property number.
#define  PROPVALUE     44  // STRING     Property value. On GDSII, 128 characters max, unless an SREF, AREF, or NODE, which may have 512 characters.
#define  BOX           45  // NO_DATA    The beginning of a BOX element. An unfilled boundary. Not used for IC polygon.
#define  BOXTYPE       46  // INTEGER_2  Boxtype specification.
#define  PLEX          47  // INTEGER_4  Plex number and plexhead flag. The least significant bit of the most significant byte is the plexhead flag.
#define  BGNEXTN       48  // INTEGER_4  Path extension beginning for pathtype 4 in CustomPlus. In database units, may be negative.
#define  ENDEXTN       49  // INTEGER_4  Path extension end for pathtype 4 in CustomPlus. In database units, may be negative.
#define  TAPENUM       50  // INTEGER_2  Tape number for multi-reel stream file.
#define  TAPECODE      51  // INTEGER_2  Tape code to verify that the reel is from the proper set. 12 bytes that are supposed to form a unique tape code.
#define  STRCLASS      52  // BIT_ARRAY  Calma use only.
#define  RESERVED      53  // INTEGER_4  Used to be NUMTYPES per GDSII Stream Format Manual, v6.0.
#define  FORMAT        54  // INTEGER_2  Archive or Filtered flag.  0: Archive 1: filtered
#define  MASK          55  // STRING     Only in filtered streams. Layers and datatypes used for mask in a filtered stream file. A string giving ranges of layers and datatypes separated by a semicolon. There may be more than one mask in a stream file.
#define  ENDMASKS      56  // NO_DATA    The end of mask descriptions.
#define  LIBDIRSIZE    57  // INTEGER_2  Number of pages in library director, a GDSII thing, it seems to have only been used when Calma INFORM was creating a new library.
#define  SRFNAME       58  // STRING     Sticks rule file name.
#define  LIBSECUR      59  // INTEGER_2  Access control list stuff for ancient CalmaDOS. INFORM used this when creating a new library. Had 1 to 32 entries with group numbers, user numbers and access rights.

#endif

//end

