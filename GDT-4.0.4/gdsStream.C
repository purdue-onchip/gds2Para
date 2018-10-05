// Copyright 1995-2014 by Ken Schumack (Schumack@cpan.org)
// STREAM stuff
// @(#) $Id: gdsStream.C 95 2014-12-08 17:16:36Z schumack $
#include <time.h>
#include <sstream>
#include <sys/file.h>
#include <fcntl.h>
#include <gdsStream.h>
#include <unistd.h>
#include <math.h>
//#include <iostream.h>
#include <iostream>
using namespace std;

#include <sstream>
#include <cstdlib>
#include <gcc4.h>
//#include <string.h>

extern char* mystrncpy(char *dest, char *src, size_t n);

static const char* gdsStream_Cwhat = "@(#) $Id: gdsStream.C 95 2014-12-08 17:16:36Z schumack $ $Revision: 95 $ $Date: 2014-12-08 11:16:36 -0600 (Mon, 08 Dec 2014) $";

GDSFILE::GDSFILE(char* fileName, int readOrWrite)
{
    int   i,j;
    char* tmpName;

    // remake FileName if found with one of these std extensions ".gds2, .gdsii, .sf, .gds"
    tmpName = new char[2048];
    mystrncpy(tmpName, fileName, 2048);
    if (! strcmp(fileName,""))
    {
        Fd = fileno(stdin);
        //setbuf(stdin, NULL);
        //setbuf(stdin, Buffer);
    }
    else
    {
        Fd = open(tmpName, O_RDONLY, 0777);
    }

    if (Fd == -1)
    {
        strcpy(tmpName, fileName);
        strcat(tmpName, ".gds2");
        Fd = open(tmpName, O_RDONLY, 0777);
    }

    if (Fd == -1)
    {
        strcpy(tmpName, fileName);
        strcat(tmpName, ".gdsii");
        Fd = open(tmpName, O_RDONLY, 0777);
    }

    if (Fd == -1)
    {
        strcpy(tmpName, fileName);
        strcat(tmpName, ".sf");
        Fd = open(tmpName, O_RDONLY, 0777);
    }

    if (Fd == -1)
    {
        strcpy(tmpName, fileName);
        strcat(tmpName, ".gds");
        Fd = open(tmpName, O_RDONLY, 0777);
    }

    FileName = new char[1024];
    mystrncpy(FileName, fileName, 1024);
    // let the opstrm handle the error if needed

    if (Fd != fileno(stdin)) close(Fd);

    Eof      = FALSE;
    EndOfLib = FALSE;
    Writtn   = readOrWrite;
    Ptr      = 0;
    for (i=0; i<NUMGDSLAYERS; i++)
    {
        Tlayers[i] = Glayers[i] = FALSE;
        for (j=0; j<NUMGDSLAYERS; j++)
        {
            LayerDataTypes[i][j] = LayerTextTypes[i][j] = FALSE;
        }
    }
    opstrm();
}


// OPSTRM - opens STREAM file for read/write
void GDSFILE::opstrm()
{
    if (Writtn == WRITE)
    {
        if (! strcmp(FileName,""))
        {
            Fd = fileno(stdout);
            #ifdef DEBUG
            setbuf(stdout, NULL);
            #endif
        }
        else
        {
            Fd = creat(FileName, 0777);
            if(Fd == -1)
            {
                cout << "ERROR ***** Unable to create file \"" << FileName << "\". Exiting..." << endl;
                exit(1);
            }
        }
        Ptr = 0;
    }
    else
    {
        if (! strcmp(FileName,""))
        {
            Fd = fileno(stdin);
            setbuf(stdin, NULL);
            //Fd = dup2(fileno(stdin),9);
        }
        else
        {
            Fd = open(FileName, O_RDONLY, 0777);
            if(Fd == -1)
            {
                cout << "ERROR ***** Unable to read file \"" << FileName << "\". Exiting..." << endl;
                exit(1);
            }
        }
        Ptr = 204800;
    }
}


// RDSTRM - reads next record from STREAM file
int GDSFILE::rdstrm()
{
    int     remain,         // remaining size of tape block
            amount_read,
            amount_read_total = 0;

    Length = 0;
    while (Length == 0)
    {
        if (Ptr > 204799)
        {
            if ((amount_read = read(Fd, Buffer, 204800)) <= 0)
            {
                Eof = TRUE;
                if (amount_read < 0)
                {
                    fprintf(stderr,"ERROR **** problem reading Fd:%d\n",Fd);

                }
                return 0;
            }
            amount_read_total += amount_read;
            Ptr = 0;
        }
        Length = (int)((unsigned short)((Buffer[Ptr]) << 8) | ((unsigned short)(Buffer[Ptr+1]) & 0xff));
        Ptr += 2;
    }

    if (Ptr > 204799)
    {
        if ((amount_read = read(Fd, Buffer, 204800)) <= 0)
        {
            Eof = TRUE;
            if (amount_read < 0)
            {
                fprintf(stderr,"ERROR **** problem reading Fd:%d\n",Fd);

            }
            return 0;
        }
        amount_read_total += amount_read;
        Ptr = 0;
    }
    Rectyp = (int)Buffer[Ptr];
    Dattyp = (int)Buffer[Ptr + 1];
    if (Rectyp==ENDLIB)
    {
        EndOfLib = TRUE;
        Eof = TRUE;
    }
    Ptr += 2;
    Length -= 4;
    if(Length >= 0)
    {
        remain = 204800 - Ptr;
        if(remain <= Length)
        {
            copy(&Buffer[Ptr], Record, remain);
            if ((amount_read = read(Fd, Buffer, 204800)) <= 0)
            {
                Eof = TRUE;
                if (amount_read < 0)
                {
                    fprintf(stderr,"ERROR **** problem reading Fd:%d\n",Fd);

                }
                return 0;
            }
            amount_read_total += amount_read;
            Ptr = 0;
        }
        if(remain < Length)
        {
            copy(Buffer, &Record[remain], (Length - remain));
            Ptr = Length - remain;
        }
        if(remain > Length)
        {
            copy(&Buffer[Ptr], Record, Length);
            Ptr += Length;
        }
    }
    Record[Length] = 0;
    return amount_read_total;
}


// CPSTRM -
void GDSFILE::cpstrm(GDSFILE* gds2file)
{
    wrstrm(gds2file->record(), gds2file);
}


// CPEND -
void GDSFILE::cpend(GDSFILE* gds2file)
{
    int errorCnt = 0;

    for(int i = Ptr; i < 204800; i++) Buffer[i] = gds2file->Buffer[i];
    if (write(Fd, Buffer, 204800) < 0) errorCnt++;
    if (Fd != fileno(stdout)) close(Fd);
}


// WRSTRM - writes private Record to STREAM file
void GDSFILE::wrstrm()
{
    int     remain;   // remaining size of "tape" block
    int     len;      // length of data + rectyp + dattyp
    int     i;
    int     errorCnt = 0;

    len = Length + 4;
    Buffer[Ptr] = (char)((len >> 8) & 0xff);

    if((len & 0x80) == 0) Buffer[Ptr+1] = (char)(len & 0xff);
    else                  Buffer[Ptr+1] = (char)(len & 0xff) | 0x80;

    Ptr += 2;
    if(Ptr > 204799)
    {
        if (write(Fd, Buffer, 204800) < 0) errorCnt++;
        Ptr = 0;
    }

    Buffer[Ptr]     = (char)Rectyp;
    Buffer[Ptr + 1] = (char)Dattyp;
    Ptr += 2;
    if(Ptr > 204799)
    {
        if (write(Fd, Buffer, 204800) < 0) errorCnt++;
        Ptr = 0;
    }

    if(Length >= 0)
    {
        remain = 204800 - Ptr;
        if(remain <= Length)
        {
            copy(Record, &Buffer[Ptr], remain);
            if (write(Fd, Buffer, 204800) < 0) errorCnt++;
            Ptr = 0;
            for(i = 0; i < 204800; Buffer[i] = 0, i++);
        }

        if(remain < Length)
        {
            copy(&Record[remain], Buffer, (Length - remain));
            Ptr = Length - remain;
        }

        if(remain > Length)
        {
            copy(Record, &Buffer[Ptr], Length);
            Ptr += Length;
        }
    }
}


// WRSTRM - writes external Record to STREAM file using Length,Rectyp, and Dattyp of existing record
void GDSFILE::wrstrm(char record[204800], GDSFILE* gds2file)
{
    int remain;   // remaining size of "tape" block
    int len;      // length of data + rectyp + dattyp
    int i;
    int errorCnt = 0;

    Length = gds2file -> Length;
    Rectyp = gds2file -> Rectyp;
    Dattyp = gds2file -> Dattyp;

    len = Length + 4;
    Buffer[Ptr] = (char)((len >> 8) & 0xff);

    if((len & 0x80) == 0) Buffer[Ptr+1] = (char)(len & 0xff);
    else                  Buffer[Ptr+1] = (char)(len & 0xff) | 0x80;

    Ptr += 2;
    if(Ptr > 204799)
    {
        if (write(Fd, Buffer, 204800) < 0) errorCnt++;
        Ptr = 0;
    }

    Buffer[Ptr]     = (char)Rectyp;
    Buffer[Ptr + 1] = (char)Dattyp;
    Ptr += 2;
    if(Ptr > 204799)
    {
        if (write(Fd, Buffer, 204800) < 0) errorCnt++;
        Ptr = 0;
    }

    if(Length >= 0)
    {
        remain = 204800 - Ptr;
        if(remain <= Length)
        {
            copy(record, &Buffer[Ptr], remain);
            if (write(Fd, Buffer, 204800) < 0) errorCnt++;
            Ptr = 0;
            for(i = 0; i < 204800; Buffer[i] = 0, i++);
        }

        if(remain < Length)
        {
            copy(&record[remain], Buffer, (Length - remain));
            Ptr = Length - remain;
        }

        if(remain > Length)
        {
            copy(record, &Buffer[Ptr], Length);
            Ptr += Length;
        }
    }
}


// WRSTRM - writes external Record to STREAM file
void GDSFILE::wrstrm(char record[204800], int rectyp, int dattyp, int length)
{
    int remain;   // remaining size of "tape" block
    int len;      // length of data + rectyp + dattyp
    int i;
    int errorCnt = 0;

    len = length + 4;
    Buffer[Ptr] = (char)((len >> 8) & 0xff);

    if((len & 0x80) == 0) Buffer[Ptr+1] = (char)(len & 0xff);
    else                  Buffer[Ptr+1] = (char)(len & 0xff) | 0x80;

    Ptr += 2;
    if(Ptr > 204799)
    {
        if (write(Fd, Buffer, 204800) < 0) errorCnt++;
        Ptr = 0;
    }

    Buffer[Ptr]     = (char)rectyp;
    Buffer[Ptr + 1] = (char)dattyp;
    Ptr += 2;
    if(Ptr > 204799)
    {
        if (write(Fd, Buffer, 204800) < 0) errorCnt++;
        Ptr = 0;
    }

    if(length >= 0)
    {
        remain = 204800 - Ptr;
        if(remain <= length)
        {
            copy(record, &Buffer[Ptr], remain);
            if (write(Fd, Buffer, 204800) < 0) errorCnt++;
            Ptr = 0;
            for(i = 0; i < 204800; Buffer[i] = 0, i++);
        }

        if(remain < length)
        {
            copy(&record[remain], Buffer, (length - remain));
            Ptr = length - remain;
        }

        if(remain > length)
        {
            copy(record, &Buffer[Ptr], length);
            Ptr += length;
        }
    }
}


// CLSTRM - closes STREAM file
void GDSFILE::clstrm()
{
    int i;              // loop counter for zero fill
    int errorCnt = 0;

    if(Writtn == WRITE) // pad w/ zeros
    {
        for(i = Ptr; i < 204800; i++) Buffer[i] = 0;

        if (write(Fd, Buffer, 204800) < 0) errorCnt++;
    }

    if (Fd != fileno(stdout)) close(Fd);
}


// writes an ENDEL to stream file
void GDSFILE::endEl()
{
    Length = 0;
    Rectyp = ENDEL;
    Dattyp = NO_DATA;
    wrstrm();
}


// WRITES A GDSII ENDLIB RECORD
void GDSFILE::endLib()
{
    // WRITE ENDLIB
    Length = 0;
    Rectyp = ENDLIB;
    Dattyp = NO_DATA;
    wrstrm();
}


// WRITES A GDSII ENDSTR RECORD
void GDSFILE::endStr()
{
    // WRITE ENDSTR
    Length = 0;
    Rectyp = ENDSTR;
    Dattyp = NO_DATA;
    wrstrm();
}


// GET_DBL - returns 64 bit real from GDS STREAM data representation
double GDSFILE::getDbl()
{
    int     i,
            byte,
            negative,
            expon;
    double  mant,
            dbl;

    byte = (int)(*Record & 0xff);
    if (byte > 127)
    {
        negative = TRUE;
        expon = byte - 192;
    }
    else
    {
        negative = FALSE;
        expon = byte - 64;
    }

    mant = 0.0;
    for (i = 1; i <= 7; i++)
    {
        byte = (int)(*(Record + i) & 0xff);
        mant = mant + ((double)byte) / pow((double)256, (double)i);
    }
    dbl = mant * pow((double)16, (double)expon);
    if (negative) dbl = -dbl;

    return dbl;
}


// GET_DBL - returns 64 bit real from GDS STREAM data representation
double GDSFILE::getDbl(int offset)
{
    int     i,
            byte,
            negative,
            expon;
    double  mant,
            dbl;

    byte = (int)(*(Record + offset) & 0xff);
    if (byte > 127)
    {
        negative = TRUE;
        expon = byte - 192;
    }
    else
    {
        negative = FALSE;
        expon = byte - 64;
    }

    mant = 0.0;
    for (i = 1; i <= 7; i++)
    {
        byte = (int)(*(Record + offset + i) & 0xff);
        mant = mant + ((double)byte) / pow((double)256, (double)i);
    }
    dbl = mant * pow((double)16, (double)expon);
    if (negative) dbl = -dbl;

    return dbl;
}


// GET_I16 - returns 16 bit integer from GDS STREAM data representation
int GDSFILE::getI16(int offset)
{
    int     negative;
    int     byte;
    int     int16;

    byte = (int)(*(Record + offset) & 0xff);

    if(byte > 127) negative = TRUE;
    else           negative = FALSE;

    int16 = byte;
    byte = (int)(*(Record + offset + 1) & 0xff);
    int16 = int16 * 256 + byte;

    if (negative) int16 = int16 - 65536;

    return(int16);
}


// GET_I16 - returns 16 bit integer from GDS STREAM data representation
int GDSFILE::getI16()
{
    int     negative;
    int     byte;
    int     int16;

    byte = (int)(*Record & 0xff);

    if(byte > 127) negative = TRUE;
    else           negative = FALSE;

    int16 = byte;
    byte = (int)(*(Record + 1) & 0xff);
    int16 = int16 * 256 + byte;

    if (negative) int16 = int16 - 65536;

    return(int16);
}


// GET_I32 - returns 32 bit integer from GDS STREAM data representation
int GDSFILE::getI32()
{
    int     i,
            negative,
            byte,
            int32;

    byte = (int)(*Record & 0xff);
    if(byte > 127)
    {
        byte = byte - 255;
        negative = TRUE;
    }
    else  negative = FALSE;

    int32 = byte;
    for(i = 1; i <= 3; i++)
        {
        byte = (int)(*(Record + i) & 0xff);
        if(negative) byte = byte - 255;
        int32 = int32 * 256 + byte;
    }

    if (negative) int32 = int32 - 1;

    return int32;
}


// GET_I32 - returns 32 bit integer from GDS STREAM data representation
int GDSFILE::getI32(int offset)
{
    int     i,
            negative,
            byte,
            int32;

    byte = (int)(*(Record + offset) & 0xff);
    if(byte > 127)
    {
        byte = byte - 255;
        negative = TRUE;
    }
    else  negative = FALSE;

    int32 = byte;
    for(i = 1; i <= 3; i++)
        {
        byte = (int)(*(Record + offset + i) & 0xff);
        if(negative) byte = byte - 255;
        int32 = int32 * 256 + byte;
    }

    if (negative) int32 = int32 - 1;

    return int32;
}


// 0 || 1 depending on whether we are at EOF
int GDSFILE::eof()
{
    return Eof;
}


// 0 || 1 depending on whether we are at EOF
int GDSFILE::endoflib()
{
    return EndOfLib;
}


// return Length of Record
int GDSFILE::length()
{
    return Length;
}


// set Length of Record
void GDSFILE::length(int len)
{
    if (len < 4)
    {
        cout << "ERROR:: Program attempted to set invalid Length" << endl;
        Length = 4;
    }
    else Length = len;
}


// return Dattyp of Record
int GDSFILE::dattyp()
{
    return Dattyp;
}


// set Dattyp of Record
void GDSFILE::dattyp(int dattype)
{
    if ((dattype > 6) || (dattype < 0))
    {
        cout << "ERROR:: Program attempted to set invalid Dattyp" << endl;
        Dattyp = NO_DATA;
    }
    else Dattyp = dattype;
}


// return Rectyp of Record
int GDSFILE::rectyp()
{
    return Rectyp;
}


// set Rectyp of Record
void GDSFILE::rectyp(int rectype)
{
    if ((rectype > 59) || (rectype < 0))
    {
        cout << "ERROR:: Program attempted to set invalid Rectyp" << endl;
        Rectyp = HEADER;
    }
    else Rectyp = rectype;
}


// return GDS Record
char* GDSFILE::record()
{
    return Record;
}


// store library name
void GDSFILE::libName(char* name)
{
    LibName = new char[strlen(name) + 1];
    strcpy(LibName, name);
}


// get stored library name
char* GDSFILE::libName()
{
    return LibName;
}


// get stored stream file name
char* GDSFILE::fileName()
{
    return FileName;
}


// use to save the fact that you found a text layer in the gds file
void GDSFILE::foundTextLayer(short layerNum)
{
    if (layerNum < NUMGDSLAYERS)
    {
        Tlayers[layerNum] = TRUE;
    }
    else
    {
        cout << "ERROR **** Found graphics layer " << layerNum << " in structure " << CurrentStrName << endl;
    }
}


// use to save the fact that you found a graphics layer in the gds file
void GDSFILE::foundGraphicsLayer(short layerNum)
{
    if (layerNum < NUMGDSLAYERS)
    {
        Glayers[layerNum] = TRUE;
    }
    else
    {
        cout << "ERROR **** Found graphics layer " << layerNum << " in structure " << CurrentStrName << endl;
    }
}


// use to save the fact that you found a datatype in the gds file
void GDSFILE::foundLayerDatatype(short layerNum, short dataTypeNum)
{
    if ((layerNum < NUMGDSLAYERS) && (dataTypeNum < NUMGDSLAYERS))
    {
        LayerDataTypes[layerNum][dataTypeNum] = TRUE;
    }
    else
    {
        cout << "ERROR **** Found graphics layer " << layerNum << " with datatype " << dataTypeNum << " in structure " << CurrentStrName << endl;
    }
}


// use to save the fact that you found a datatype in the gds file
void GDSFILE::foundLayerTexttype(short layerNum, short textTypeNum)
{
    if ((layerNum < NUMGDSLAYERS) && (textTypeNum < NUMGDSLAYERS))
    {
        LayerTextTypes[layerNum][textTypeNum] = TRUE;
    }
    else
    {
        cout << "ERROR **** Found graphics layer " << layerNum << " with texttype " << textTypeNum << " in structure " << CurrentStrName << endl;
    }
}


// true or false .. does this graphics layer exist in the stream file?
int GDSFILE::gLayer(short layerNum)
{
    if ((layerNum < NUMGDSLAYERS) && Glayers[layerNum]) return TRUE;
    else return FALSE;
}


// true or false .. does this graphicsLayer/dataType exist in the stream file?
int GDSFILE::layerDataType(short layerNum, short dataType)
{
    if ((layerNum < NUMGDSLAYERS) && (dataType < NUMGDSLAYERS) &&
             LayerDataTypes[layerNum][dataType]) return TRUE;
    else return FALSE;
}


// true or false .. does this text layer exist in the stream file?
int GDSFILE::tLayer(short layerNum)
{
    if ((layerNum < NUMGDSLAYERS) && Tlayers[layerNum]) return TRUE;
    else return FALSE;
}


// true or false .. does this textLayer/textType exist in the stream file?
int GDSFILE::layerTextType(short layerNum, short textType)
{
    if ((layerNum < NUMGDSLAYERS) && (textType < NUMGDSLAYERS) &&
             LayerTextTypes[layerNum][textType]) return TRUE;
    else return FALSE;
}


//
void GDSFILE::copy(char src_rec[],    // source record
                   char dst_rec[],    // destination record
                   int  num)          // number of chars for copy
{
    for(int i = 0; i < num; i++)  dst_rec[i] = src_rec[i];
}

/// PUTS AREF IN LIBRARY {like Calma's AREF command}
void GDSFILE::putAref(
         char*  sname,
         unsigned short ref,   // 1 for reflection, 0 for no reflection
         double mag,
         double angle,
         short  col,
         short  row,
         double x1, double y1, double x2, double y2, double x3, double y3, // x1y1:Origin, x2y2:Column, x3y3:Row
         int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING],
         double dbu_uu
)
{
    int index;
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE AREF
    Length = 0;
    Rectyp = AREF;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE SNAME
    strcpy(Record, sname);
    Length = strlen(Record);
    if (Length%2)
    {
        Record[Length] = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = SNAME;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE COLROW
    Length = 4;
    Rectyp = COLROW;
    Dattyp = INTEGER_2;
    putI16(col, 0);
    putI16(row, 2);
    wrstrm();

    // WRITE XY
    Length = 24;
    Rectyp = XY;
    Dattyp = INTEGER_4;

    if (x1 >= 0.0) { putI32((long int)((x1 + epsilon) * factor), 0); }
    else           { putI32((long int)((x1 - epsilon) * factor), 0); }

    if (y1 >= 0.0) { putI32((long int)((y1 + epsilon) * factor), 4); }
    else           { putI32((long int)((y1 - epsilon) * factor), 4); }

    if (x2 >= 0.0) { putI32((long int)((x2 + epsilon) * factor), 8); }
    else           { putI32((long int)((x2 - epsilon) * factor), 8); }

    if (y2 >= 0.0) { putI32((long int)((y2 + epsilon) * factor), 12); }
    else           { putI32((long int)((y2 - epsilon) * factor), 12); }

    if (x3 >= 0.0) { putI32((long int)((x3 + epsilon) * factor), 16); }
    else           { putI32((long int)((x3 - epsilon) * factor), 16); }

    if (y3 >= 0.0) { putI32((long int)((y3 + epsilon) * factor), 20); }
    else           { putI32((long int)((y3 - epsilon) * factor), 20); }

    wrstrm();

    for (index=0; index <= propIndex; index++)
    {
        Length = 2;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPATTR;
        Dattyp = INTEGER_2;
        putI16(propNumArray[index], 0);
        wrstrm();

        Length = 4;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPVALUE;
        Dattyp = ACSII_STRING;
        strcpy(Record, propValueArray[index]);
        Length = strlen(Record);
        if (Length%2) {
            Record[Length]     = '\0';
            Record[Length + 1] = '\0';
            Length++;
        }
        wrstrm();
    }

    // WRITE ENDEL
    endEl();
}


/// PUTS AREF IN LIBRARY {like Calma's AREF command}
void GDSFILE::putAref(
         char*  sname,
         unsigned short ref,   // 1 for reflection, 0 for no reflection
         double mag,
         double angle,
         short  col,
         short  row,
         double x1, double y1, double x2, double y2, double x3, double y3, // x1y1:Origin, x2y2:Column, x3y3:Row
         double dbu_uu
)
{
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE AREF
    Length = 0;
    Rectyp = AREF;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE SNAME
    strcpy(Record, sname);
    Length = strlen(Record);
    if (Length%2)
    {
        Record[Length] = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = SNAME;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE COLROW
    Length = 4;
    Rectyp = COLROW;
    Dattyp = INTEGER_2;
    putI16(col, 0);
    putI16(row, 2);
    wrstrm();

    // WRITE XY
    Length = 24;
    Rectyp = XY;
    Dattyp = INTEGER_4;

    if (x1 >= 0.0) { putI32((long int)((x1 + epsilon) * factor), 0); }
    else           { putI32((long int)((x1 - epsilon) * factor), 0); }

    if (y1 >= 0.0) { putI32((long int)((y1 + epsilon) * factor), 4); }
    else           { putI32((long int)((y1 - epsilon) * factor), 4); }

    if (x2 >= 0.0) { putI32((long int)((x2 + epsilon) * factor), 8); }
    else           { putI32((long int)((x2 - epsilon) * factor), 8); }

    if (y2 >= 0.0) { putI32((long int)((y2 + epsilon) * factor), 12); }
    else           { putI32((long int)((y2 - epsilon) * factor), 12); }

    if (x3 >= 0.0) { putI32((long int)((x3 + epsilon) * factor), 16); }
    else           { putI32((long int)((x3 - epsilon) * factor), 16); }

    if (y3 >= 0.0) { putI32((long int)((y3 + epsilon) * factor), 20); }
    else           { putI32((long int)((y3 - epsilon) * factor), 20); }

    wrstrm();

    // WRITE ENDEL
    endEl();
}


/// PUTS AREF IN LIBRARY {like Calma's AREF command}
void GDSFILE::putAref(
         char*  sname,
         unsigned short ref,   // 1 for reflection, 0 for no reflection
         double mag,
         double angle,
         short  col,
         short  row,
         double x1, double y1, double x2, double y2, double x3, double y3 // x1y1:Origin, x2y2:Column, x3y3:Row
)
{
    // WRITE AREF
    Length = 0;
    Rectyp = AREF;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE SNAME
    strcpy(Record, sname);
    Length = strlen(Record);
    if (Length%2)
    {
        Record[Length] = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = SNAME;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE COLROW
    Length = 4;
    Rectyp = COLROW;
    Dattyp = INTEGER_2;
    putI16(col, 0);
    putI16(row, 2);
    wrstrm();

    // WRITE XY
    Length = 24;
    Rectyp = XY;
    Dattyp = INTEGER_4;

    if (x1 >= 0.0) { putI32((long int)((x1 + G_epsilon) * 1000), 0); }
    else           { putI32((long int)((x1 - G_epsilon) * 1000), 0); }

    if (y1 >= 0.0) { putI32((long int)((y1 + G_epsilon) * 1000), 4); }
    else           { putI32((long int)((y1 - G_epsilon) * 1000), 4); }

    if (x2 >= 0.0) { putI32((long int)((x2 + G_epsilon) * 1000), 8); }
    else           { putI32((long int)((x2 - G_epsilon) * 1000), 8); }

    if (y2 >= 0.0) { putI32((long int)((y2 + G_epsilon) * 1000), 12); }
    else           { putI32((long int)((y2 - G_epsilon) * 1000), 12); }

    if (x3 >= 0.0) { putI32((long int)((x3 + G_epsilon) * 1000), 16); }
    else           { putI32((long int)((x3 - G_epsilon) * 1000), 16); }

    if (y3 >= 0.0) { putI32((long int)((y3 + G_epsilon) * 1000), 20); }
    else           { putI32((long int)((y3 - G_epsilon) * 1000), 20); }

    wrstrm();

    // WRITE ENDEL
    endEl();
}


// PUTDBL - puts 64 bit real in stream output buffer
void GDSFILE::putDbl(double  dbl, int offset)
{
    int     negative;
    double  r;
    int     expon;
    int     i;
    int     byte;

    if(dbl < 0.0)
    {
        negative = TRUE;
        r = -dbl;
    }
    else
    {
        negative = FALSE;
        r = dbl;
    }

    expon = 0;
    while(r >= 1.0)
    {
        expon++;
        r = r / 16.0;
    }

    if (r != 0)
    {
        while(r < 0.0625)
        {
            expon--;
            r = r * 16.0;
        }
    }

    if(negative == 1) expon = 192 + expon;
    else              expon =  64 + expon;

    *(Record + offset) = (char)expon;

    for(i = 1; i <= 7; i++)
    {
        byte = (int)(r*256.0);
        *(Record + offset + i) = (char)byte;
        r = r * 256.0 - (double)byte;
    }
}


// PUTI32 - puts 32 bit integer in stream output buffer
void GDSFILE::putI32(int i32, int offset)
{
    int     negative,
            rem,
            i,
            byte,
            fact;

    if(i32 < 0)
    {
        negative = TRUE;
        rem = -i32 -1;
    }
    else
    {
        negative = FALSE;
        rem = i32;
    }

    fact = 256 * 256 * 256;
    for(i = 3; i >= 0; i--)
    {
        byte = rem / fact;
        rem = rem - byte * fact;
        if(negative == 1) byte = 255 - byte;
        *(Record + offset + 3 - i) = (char)byte;
        fact = fact / 256;
    }
}


// PUTI16 - puts 16 bit integer in stream output buffer
void GDSFILE::putI16(unsigned short i16, int offset)
{
    unsigned short     rem;

    rem = i16;
    *(Record + offset)     = (char)(rem / 256);
    *(Record + offset + 1) = (char)(rem % 256);
}

// INITIALIZES LIBRARY HEADER, BGNLIB, LIBNAME, UNITS
//    {for Calma compabability end the library name with ".DB"}
void GDSFILE::initLib(char *library, double dbu_uu, double dbu_m,
    int myear, int mmon, int mmday, int mhour, int mmin, int msec,
    int ayear, int amon, int amday, int ahour, int amin, int asec,
    int version)
{
    // WRITE HEADER
    Length = 2;
    Rectyp = HEADER;
    Dattyp = INTEGER_2;
    putI16(version, 0);    // writing release 3 type stuff
    wrstrm();

    // WRITE BGNLIB
    Length = 24;
    Rectyp = BGNLIB;
    Dattyp = INTEGER_2;
    if (myear > 1900) myear -= 1900;
    putI16(myear, 0);   // modification time
    putI16(mmon,  2);
    putI16(mmday, 4);
    putI16(mhour, 6);
    putI16(mmin,  8);
    putI16(msec, 10);
    if (ayear > 1900) ayear -= 1900;
    putI16(ayear, 12);  // last access time
    putI16(amon,  14);
    putI16(amday, 16);
    putI16(ahour, 18);
    putI16(amin,  20);
    putI16(asec,  22);
    wrstrm();

    // WRITE LIBNAME
    strcpy(Record, library);
    Length = strlen(Record);
    if (Length%2)
    {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = LIBNAME;
    Dattyp = ACSII_STRING;
    wrstrm();
    libName(library);

    // WRITE UNITS
    Length = 16;
    Rectyp = UNITS;
    Dattyp = REAL_8;
    putDbl(dbu_uu, 0);     // Calma default is 1.0e-3
    putDbl(dbu_m,  8);     // Calma default is 1.0e-9
    wrstrm();
}


// INITIALIZES LIBRARY HEADER, BGNLIB, LIBNAME, UNITS
//    {for Calma compabability end the library name with ".DB"}
void GDSFILE::initLib(char *library, double dbu_uu, double dbu_m,
    int myear, int mmon, int mmday, int mhour, int mmin, int msec,
    int ayear, int amon, int amday, int ahour, int amin, int asec)
{
    GDSFILE::initLib(library,dbu_uu,dbu_m,myear,mmon,mmday,mhour,mmin,msec,ayear,amon,amday,ahour,amin,asec,3);
}


// INITIALIZES LIBRARY HEADER, BGNLIB, LIBNAME, UNITS
//    {for Calma compabability end the library name with ".DB"}
void GDSFILE::initLib(char *library, double dbu_uu, double dbu_m)
{
    struct  tm   *ts;
    time_t       time_val;

    // WRITE HEADER
    Length = 2;
    Rectyp = HEADER;
    Dattyp = INTEGER_2;
    putI16(3, 0);    // writing release 3 type stuff
    wrstrm();

    // WRITE BGNLIB
    Length = 24;
    Rectyp = BGNLIB;
    Dattyp = INTEGER_2;
    time(&time_val);
    ts = localtime(&time_val);
    putI16(ts->tm_year,    0);   // modification time
    putI16(ts->tm_mon + 1, 2);
    putI16(ts->tm_mday,    4);
    putI16(ts->tm_hour,    6);
    putI16(ts->tm_min,     8);
    putI16(ts->tm_sec,     10);
    putI16(ts->tm_year,    12);  // last access time
    putI16(ts->tm_mon + 1, 14);
    putI16(ts->tm_mday,    16);
    putI16(ts->tm_hour,    18);
    putI16(ts->tm_min,     20);
    putI16(ts->tm_sec,     22);
    wrstrm();

    // WRITE LIBNAME
    strcpy(Record, library);
    Length = strlen(Record);
    if (Length%2)
    {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = LIBNAME;
    Dattyp = ACSII_STRING;
    wrstrm();
    libName(library);

    // WRITE UNITS
    Length = 16;
    Rectyp = UNITS;
    Dattyp = REAL_8;
    putDbl(dbu_uu, 0);     // Calma default is 1.0e-3
    putDbl(dbu_m,  8);     // Calma default is 1.0e-9
    wrstrm();
}


// INITIALIZES LIBRARY HEADER, BGNLIB, LIBNAME, UNITS
//    {for Calma compabability end the library name with ".DB"}
void GDSFILE::initLib(char *library)
{
    struct  tm   *ts;
    time_t       time_val;

    // WRITE HEADER
    Length = 2;
    Rectyp = HEADER;
    Dattyp = INTEGER_2;
    putI16(3, 0);    // writing release 3 type stuff
    wrstrm();

    // WRITE BGNLIB
    Length = 24;
    Rectyp = BGNLIB;
    Dattyp = INTEGER_2;
    time(&time_val);
    ts = localtime(&time_val);
    putI16(ts->tm_year,    0);   // modification time
    putI16(ts->tm_mon + 1, 2);
    putI16(ts->tm_mday,    4);
    putI16(ts->tm_hour,    6);
    putI16(ts->tm_min,     8);
    putI16(ts->tm_sec,     10);
    putI16(ts->tm_year,    12);  // last access time
    putI16(ts->tm_mon + 1, 14);
    putI16(ts->tm_mday,    16);
    putI16(ts->tm_hour,    18);
    putI16(ts->tm_min,     20);
    putI16(ts->tm_sec,     22);
    wrstrm();

    // WRITE LIBNAME
    strcpy(Record, library);
    Length = strlen(Record);
    if (Length%2)
    {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = LIBNAME;
    Dattyp = ACSII_STRING;
    wrstrm();
    libName(library);

    // WRITE UNITS
    Length = 16;
    Rectyp = UNITS;
    Dattyp = REAL_8;
    putDbl(1.0e-3, 0);     // Calma default is 1.0e-3
    putDbl(1.0e-9, 8);     // Calma default is 1.0e-9
    wrstrm();
}


// PUTS A RECTANGULAR BOUNDARY ON A SPECIFIED LAYER {like calma's RT command}
void GDSFILE::putRt(int    layer,
       int    datatyp,
       double minX, double minY,
       double maxX, double maxY,
       double dbu_uu)
{
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE BOUNDARY
    Length = 0;
    Rectyp = BOUNDARY;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE XY
    Length = 40;                    // 10 4 byte records
    Rectyp = XY;
    Dattyp = INTEGER_4;

    if (minX >= 0.0) { putI32((long int)((minX + epsilon) * factor),  0); }   // lower left corner
    else             { putI32((long int)((minX - epsilon) * factor),  0); }   // lower left corner

    if (minY >= 0.0) { putI32((long int)((minY + epsilon) * factor),  4); }
    else             { putI32((long int)((minY - epsilon) * factor),  4); }

    if (minX >= 0.0) { putI32((long int)((minX + epsilon) * factor),  8); }   // upper left corner
    else             { putI32((long int)((minX - epsilon) * factor),  8); }   // upper left corner

    if (maxY >= 0.0) { putI32((long int)((maxY + epsilon) * factor), 12); }
    else             { putI32((long int)((maxY - epsilon) * factor), 12); }

    if (maxX >= 0.0) { putI32((long int)((maxX + epsilon) * factor), 16); }   // upper right corner
    else             { putI32((long int)((maxX - epsilon) * factor), 16); }   // upper right corner

    if (maxY >= 0.0) { putI32((long int)((maxY + epsilon) * factor), 20); }
    else             { putI32((long int)((maxY - epsilon) * factor), 20); }

    if (maxX >= 0.0) { putI32((long int)((maxX + epsilon) * factor), 24); }   // lower right corner
    else             { putI32((long int)((maxX - epsilon) * factor), 24); }   // lower right corner

    if (minY >= 0.0) { putI32((long int)((minY + epsilon) * factor), 28); }
    else             { putI32((long int)((minY - epsilon) * factor), 28); }

    if (minX >= 0.0) { putI32((long int)((minX + epsilon) * factor), 32); }   // lower left corner again
    else             { putI32((long int)((minX - epsilon) * factor), 32); }   // lower left corner again

    if (minY >= 0.0) { putI32((long int)((minY + epsilon) * factor), 36); }
    else             { putI32((long int)((minY - epsilon) * factor), 36); }

    wrstrm();

    // WRITE ENDEL
    endEl();
}


// PUTS A RECTANGULAR BOUNDARY ON A SPECIFIED LAYER {like calma's RT command}
void GDSFILE::putRt(int    layer,
       int    datatyp,
       double minX, double minY,
       double maxX, double maxY)
{
    // WRITE BOUNDARY
    Length = 0;
    Rectyp = BOUNDARY;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE XY
    Length = 40;                    // 10 4 byte records
    Rectyp = XY;
    Dattyp = INTEGER_4;

    if (minX >= 0.0) { putI32((long int)((minX + G_epsilon) * 1000),  0); }   // lower left corner
    else             { putI32((long int)((minX - G_epsilon) * 1000),  0); }   // lower left corner

    if (minY >= 0.0) { putI32((long int)((minY + G_epsilon) * 1000),  4); }
    else             { putI32((long int)((minY - G_epsilon) * 1000),  4); }

    if (minX >= 0.0) { putI32((long int)((minX + G_epsilon) * 1000),  8); }   // upper left corner
    else             { putI32((long int)((minX - G_epsilon) * 1000),  8); }   // upper left corner

    if (maxY >= 0.0) { putI32((long int)((maxY + G_epsilon) * 1000), 12); }
    else             { putI32((long int)((maxY - G_epsilon) * 1000), 12); }

    if (maxX >= 0.0) { putI32((long int)((maxX + G_epsilon) * 1000), 16); }   // upper right corner
    else             { putI32((long int)((maxX - G_epsilon) * 1000), 16); }   // upper right corner

    if (maxY >= 0.0) { putI32((long int)((maxY + G_epsilon) * 1000), 20); }
    else             { putI32((long int)((maxY - G_epsilon) * 1000), 20); }

    if (maxX >= 0.0) { putI32((long int)((maxX + G_epsilon) * 1000), 24); }   // lower right corner
    else             { putI32((long int)((maxX - G_epsilon) * 1000), 24); }   // lower right corner

    if (minY >= 0.0) { putI32((long int)((minY + G_epsilon) * 1000), 28); }
    else             { putI32((long int)((minY - G_epsilon) * 1000), 28); }

    if (minX >= 0.0) { putI32((long int)((minX + G_epsilon) * 1000), 32); }   // lower left corner again
    else             { putI32((long int)((minX - G_epsilon) * 1000), 32); }   // lower left corner again

    if (minY >= 0.0) { putI32((long int)((minY + G_epsilon) * 1000), 36); }
    else             { putI32((long int)((minY - G_epsilon) * 1000), 36); }

    // WRITE ENDEL
    endEl();
}

// PUTS SREF IN LIBRARY {like Calma's SREF command}
void GDSFILE::putSref( char*  sname,
          unsigned short ref,  // 0 or 1
          double mag,
          double angle,
          double x_coord,
          double y_coord,
          int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING],
          double dbu_uu)
{
    int index;
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE SREF
    Length = 0;
    Rectyp = SREF;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE SNAME
    strcpy(Record, sname);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = SNAME;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;

    if (x_coord >= 0.0) { putI32((long int)((x_coord + epsilon) * factor), 0); }
    else                { putI32((long int)((x_coord - epsilon) * factor), 0); }
    if (y_coord >= 0.0) { putI32((long int)((y_coord + epsilon) * factor), 4); }
    else                { putI32((long int)((y_coord - epsilon) * factor), 4); }
    wrstrm();

    for (index=0; index <= propIndex; index++)
    {
        Length = 2;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPATTR;
        Dattyp = INTEGER_2;
        putI16(propNumArray[index], 0);
        wrstrm();

        Length = 4;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPVALUE;
        Dattyp = ACSII_STRING;
        strcpy(Record, propValueArray[index]);
        Length = strlen(Record);
        if (Length%2) {
            Record[Length]     = '\0';
            Record[Length + 1] = '\0';
            Length++;
        }
        wrstrm();
    }

    // WRITE ENDEL
    endEl();
}

// PUTS SREF IN LIBRARY {like Calma's SREF command}
void GDSFILE::putSref( char*  sname,
          unsigned short ref,  // 0 or 1
          double mag,
          double angle,
          double x_coord,
          double y_coord,
          double dbu_uu)
{
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE SREF
    Length = 0;
    Rectyp = SREF;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE SNAME
    strcpy(Record, sname);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = SNAME;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;

    if (x_coord >= 0.0) { putI32((long int)((x_coord + epsilon) * factor), 0); }
    else                { putI32((long int)((x_coord - epsilon) * factor), 0); }
    if (y_coord >= 0.0) { putI32((long int)((y_coord + epsilon) * factor), 4); }
    else                { putI32((long int)((y_coord - epsilon) * factor), 4); }
    wrstrm();

    // WRITE ENDEL
    endEl();
}


// PUTS SREF IN LIBRARY {like Calma's SREF command}
void GDSFILE::putSref( char*  sname,
          unsigned short ref,  // 0 or 1
          double mag,
          double angle,
          double x_coord,
          double y_coord)
{
    // WRITE SREF
    Length = 0;
    Rectyp = SREF;
    Dattyp = NO_DATA;
    wrstrm();
    // WRITE SNAME
    strcpy(Record, sname);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = SNAME;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;
    if (x_coord >= 0.0) { putI32((long int)((x_coord + G_epsilon) * 1000), 0); }
    else                { putI32((long int)((x_coord - G_epsilon) * 1000), 0); }
    if (y_coord >= 0.0) { putI32((long int)((y_coord + G_epsilon) * 1000), 4); }
    else                { putI32((long int)((y_coord - G_epsilon) * 1000), 4); }
    wrstrm();

    // WRITE ENDEL
    endEl();
}


// INITIALIZE STRUCTURE (like Calma's BSTRUCT)
void GDSFILE::beginStr(char *str_name,
    int myear, int mmon, int mmday, int mhour, int mmin, int msec,
    int ayear, int amon, int amday, int ahour, int amin, int asec)
{
    time_t       time_val;
    struct  tm   *ts;

    // WRITE BGNSTR
    Length = 24;
    Rectyp = BGNSTR;
    Dattyp = INTEGER_2;
    time(&time_val);
    ts = localtime(&time_val);
    putI16(myear,0);   // modification time
    putI16(mmon, 2);
    putI16(mmday,4);
    putI16(mhour,6);
    putI16(mmin, 8);
    putI16(msec, 10);
    putI16(ayear,12);   // last access time
    putI16(amon, 14);
    putI16(amday,16);
    putI16(ahour,18);
    putI16(amin, 20);
    putI16(asec, 22);
    wrstrm();

    // WRITE STRNAME
    strcpy(Record, str_name);
    Length =  strlen(Record);
    if (Length%2)
    {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRNAME;
    Dattyp = ACSII_STRING;
    wrstrm();
}


// INITIALIZE STRUCTURE (like Calma's BSTRUCT)
void GDSFILE::beginStr(char *str_name)
{
    time_t       time_val;
    struct  tm   *ts;

    // WRITE BGNSTR
    Length = 24;
    Rectyp = BGNSTR;
    Dattyp = INTEGER_2;
    time(&time_val);
    ts = localtime(&time_val);
    putI16(ts->tm_year,     0);   // modification time
    putI16(ts->tm_mon + 1,  2);
    putI16(ts->tm_mday,     4);
    putI16(ts->tm_hour,     6);
    putI16(ts->tm_min,      8);
    putI16(ts->tm_sec,     10);
    putI16(ts->tm_year,    12);   // last access time
    putI16(ts->tm_mon + 1, 14);
    putI16(ts->tm_mday,    16);
    putI16(ts->tm_hour,    18);
    putI16(ts->tm_min,     20);
    putI16(ts->tm_sec,     22);
    wrstrm();

    // WRITE STRNAME
    strcpy(Record, str_name);
    Length =  strlen(Record);
    if (Length%2)
    {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRNAME;
    Dattyp = ACSII_STRING;
    wrstrm();
}


//*********** putText *********************
//* PUTS A TEXT RECORD ON SPECIFIED LAYER *
//* layer must be >=0 && <=255            *
//*****************************************
void GDSFILE::putText(unsigned short layer,
                      unsigned short ref,  // 1 or 0
                      double         mag,
                      double         angle,
                      double x, double y,
                      char*          txt,
                      int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING],
                      double dbu_uu)
{
    int      index;
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;
    unsigned short int tt = 0;

    // WRITE TEXT
    Length = 0;
    Rectyp = TEXT;
    Dattyp = NO_DATA;
    wrstrm();

     // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE TEXTTYPE
    Length = 2;
    Rectyp = TEXTTYPE;
    Dattyp = INTEGER_2;
    putI16(tt, 0);
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;
    if (x >= 0.0) { putI32((long int)((x + epsilon) * factor), 0); }
    else          { putI32((long int)((x - epsilon) * factor), 0); }
    if (y >= 0.0) { putI32((long int)((y + epsilon) * factor), 4); }
    else          { putI32((long int)((y - epsilon) * factor), 4); }
    wrstrm();

    // WRITE STRING
    strcpy(Record, txt);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRING;
    Dattyp = ACSII_STRING;
    wrstrm();

    for (index=0; index <= propIndex; index++)
    {
        Length = 2;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPATTR;
        Dattyp = INTEGER_2;
        putI16(propNumArray[index], 0);
        wrstrm();

        Length = 4;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPVALUE;
        Dattyp = ACSII_STRING;
        strcpy(Record, propValueArray[index]);
        Length = strlen(Record);
        if (Length%2) {
            Record[Length]     = '\0';
            Record[Length + 1] = '\0';
            Length++;
        }
        wrstrm();
    }

    // WRITE ENDEL
    endEl();
}


//*********** putText *********************
//* PUTS A TEXT RECORD ON SPECIFIED LAYER *
//* layer must be >=0 && <=255            *
//*****************************************
void GDSFILE::putText(unsigned short layer,
                      unsigned short ref,  // 1 or 0
                      double         mag,
                      double         angle,
                      double x, double y,
                      char*          txt,
                      double dbu_uu)
{
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;
    unsigned short int tt = 0;

    // WRITE TEXT
    Length = 0;
    Rectyp = TEXT;
    Dattyp = NO_DATA;
    wrstrm();

     // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE TEXTTYPE
    Length = 2;
    Rectyp = TEXTTYPE;
    Dattyp = INTEGER_2;
    putI16(tt, 0);
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;
    if (x >= 0.0) { putI32((long int)((x + epsilon) * factor), 0); }
    else          { putI32((long int)((x - epsilon) * factor), 0); }
    if (y >= 0.0) { putI32((long int)((y + epsilon) * factor), 4); }
    else          { putI32((long int)((y - epsilon) * factor), 4); }
    wrstrm();

    // WRITE STRING
    strcpy(Record, txt);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRING;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE ENDEL
    endEl();
}


//*********** putText *********************
//* PUTS A TEXT RECORD ON SPECIFIED LAYER *
//* layer must be >=0 && <=255            *
//*****************************************
void GDSFILE::putText(unsigned short layer,
                      unsigned short ref,  // 1 or 0
                      double         mag,
                      double         angle,
                      double x, double y,
                      char*          txt)
{
    unsigned short int tt = 0;

    // WRITE TEXT
    Length = 0;
    Rectyp = TEXT;
    Dattyp = NO_DATA;
    wrstrm();

     // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE TEXTTYPE
    Length = 2;
    Rectyp = TEXTTYPE;
    Dattyp = INTEGER_2;
    putI16(tt, 0);
    wrstrm();

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    Length = 8;
    Rectyp = MAG;
    Dattyp = REAL_8;
    putDbl(mag, 0);
    wrstrm();

    // WRITE ANGLE
    Length = 8;
    Rectyp = ANGLE;
    Dattyp = REAL_8;
    putDbl(angle, 0);
    wrstrm();

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;
    if (x >= 0.0) { putI32((long int)((x + G_epsilon) * 1000), 0); }
    else          { putI32((long int)((x - G_epsilon) * 1000), 0); }
    if (y >= 0.0) { putI32((long int)((y + G_epsilon) * 1000), 4); }
    else          { putI32((long int)((y - G_epsilon) * 1000), 4); }
    wrstrm();

    // WRITE STRING
    strcpy(Record, txt);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRING;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE ENDEL
    endEl();
}


//*********** putText *********************
//* PUTS A TEXT RECORD ON SPECIFIED LAYER *
//* layer must be >=0 && <=64K            *
//*****************************************
//#  <text>::=          TEXT [ELFLAGS] [PLEX] LAYER <textbody>                   #
//#                                                                              #
//#  <textbody>::=      TEXTTYPE [PRESENTATION] [PATHTYPE] [WIDTH] [<strans>] XY #
//#                     STRING                                                   #
//#                                                                              #
//#  <strans>::=        STRANS [MAG] [ANGLE]                                     #
void GDSFILE::putText(unsigned short layer,
                      unsigned short textType,
                      unsigned short fontType,
                      char*          textJust,
                      unsigned short pathType,
                      double         width,
                      unsigned short ref,
                      double         mag,
                      double         angle,
                      double         x,
                      double         y,
                      char*          txt,
                      int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING],
                      double dbu_uu)
{
    int    index;
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;
    unsigned short int presentation = 0;
    double   fstep;
    long int istep;

    // WRITE TEXT
    Length = 0;
    Rectyp = TEXT;
    Dattyp = NO_DATA;
    wrstrm();

     // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE TEXTTYPE
    Length = 2;
    Rectyp = TEXTTYPE;
    Dattyp = INTEGER_2;
    putI16(textType, 0);
    wrstrm();

    if (fontType)
    {
        if (fontType == 3)
        {
            presentation |= 0x30;
        }
        else if (fontType == 2)
        {
            presentation |= 0x20;
        }
        else if (fontType == 1)
        {
            presentation |= 0x10;
        }
    }
    if (strlen(textJust) == 2)
    {
        if (strncmp(textJust,"tl",2)) // default
        {
            if (textJust[0] == 'b')
            {
                presentation |= 0x8;
            }
            else if (textJust[0] == 'm')
            {
                presentation |= 0x4;
            }

            if (textJust[1] == 'r')
            {
                presentation |= 0x2;
            }
            else if (textJust[1] == 'c')
            {
                presentation |= 0x1;
            }
        }
    }
    if (presentation)
    {
        // WRITE PRESENTATION
        Length = 2;
        Rectyp = PRESENTATION;
        Dattyp = BIT_ARRAY;             // bit array
        putI16(presentation, 0);
        wrstrm();
    }
    // WRITE PATHTYPE
    Length = 2;
    Rectyp = PATHTYPE;
    Dattyp = INTEGER_2;
    putI16(pathType, 0);
    wrstrm();

    fstep = (double) (width + epsilon) * factor; // done in steps because of previous compilier bug
    istep = (long int) fstep;
    if (istep)
    {
        // WRITE WIDTH
        Length = 4;
        Rectyp = WIDTH;
        Dattyp = INTEGER_4;
        putI32(istep, 0);
        wrstrm();
    }

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    if (mag != 1.0)
    {
        Length = 8;
        Rectyp = MAG;
        Dattyp = REAL_8;
        putDbl(mag, 0);
        wrstrm();
    }

    // WRITE ANGLE
    if (angle != 0.0)
    {
        Length = 8;
        Rectyp = ANGLE;
        Dattyp = REAL_8;
        putDbl(angle, 0);
        wrstrm();
    }

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;
    if (x >= 0.0) { putI32((long int)((x + epsilon) * factor), 0); }
    else          { putI32((long int)((x - epsilon) * factor), 0); }
    if (y >= 0.0) { putI32((long int)((y + epsilon) * factor), 4); }
    else          { putI32((long int)((y - epsilon) * factor), 4); }
    wrstrm();

    // WRITE STRING
    strcpy(Record, txt);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRING;
    Dattyp = ACSII_STRING;
    wrstrm();

    for (index=0; index <= propIndex; index++)
    {
        Length = 2;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPATTR;
        Dattyp = INTEGER_2;
        putI16(propNumArray[index], 0);
        wrstrm();

        Length = 4;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPVALUE;
        Dattyp = ACSII_STRING;
        strcpy(Record, propValueArray[index]);
        Length = strlen(Record);
        if (Length%2) {
            Record[Length]     = '\0';
            Record[Length + 1] = '\0';
            Length++;
        }
        wrstrm();
    }

    // WRITE ENDEL
    endEl();
}


//*********** putText *********************
//* PUTS A TEXT RECORD ON SPECIFIED LAYER *
//* layer must be >=0 && <=64K            *
//*****************************************
//#  <text>::=          TEXT [ELFLAGS] [PLEX] LAYER <textbody>                   #
//#                                                                              #
//#  <textbody>::=      TEXTTYPE [PRESENTATION] [PATHTYPE] [WIDTH] [<strans>] XY #
//#                     STRING                                                   #
//#                                                                              #
//#  <strans>::=        STRANS [MAG] [ANGLE]                                     #
void GDSFILE::putText(unsigned short layer,
                      unsigned short textType,
                      unsigned short fontType,
                      char*          textJust,
                      unsigned short pathType,
                      double         width,
                      unsigned short ref,
                      double         mag,
                      double         angle,
                      double         x,
                      double         y,
                      char*          txt,
                      double dbu_uu)
{
    double factor = (1 / dbu_uu);
    double epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;
    unsigned short int presentation = 0;
    double   fstep;
    long int istep;

    // WRITE TEXT
    Length = 0;
    Rectyp = TEXT;
    Dattyp = NO_DATA;
    wrstrm();

     // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE TEXTTYPE
    Length = 2;
    Rectyp = TEXTTYPE;
    Dattyp = INTEGER_2;
    putI16(textType, 0);
    wrstrm();

    if (fontType)
    {
        if (fontType == 3)
        {
            presentation |= 0x30;
        }
        else if (fontType == 2)
        {
            presentation |= 0x20;
        }
        else if (fontType == 1)
        {
            presentation |= 0x10;
        }
    }
    if (strlen(textJust) == 2)
    {
        if (strncmp(textJust,"tl",2)) // default
        {
            if (textJust[0] == 'b')
            {
                presentation |= 0x8;
            }
            else if (textJust[0] == 'm')
            {
                presentation |= 0x4;
            }

            if (textJust[1] == 'r')
            {
                presentation |= 0x2;
            }
            else if (textJust[1] == 'c')
            {
                presentation |= 0x1;
            }
        }
    }
    if (presentation)
    {
        // WRITE PRESENTATION
        Length = 2;
        Rectyp = PRESENTATION;
        Dattyp = BIT_ARRAY;             // bit array
        putI16(presentation, 0);
        wrstrm();
    }
    // WRITE PATHTYPE
    Length = 2;
    Rectyp = PATHTYPE;
    Dattyp = INTEGER_2;
    putI16(pathType, 0);
    wrstrm();

    fstep = (double) (width + epsilon) * factor; // done in steps because of previous compilier bug
    istep = (long int) fstep;
    if (istep)
    {
        // WRITE WIDTH
        Length = 4;
        Rectyp = WIDTH;
        Dattyp = INTEGER_4;
        putI32(istep, 0);
        wrstrm();
    }

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    if (mag != 1.0)
    {
        Length = 8;
        Rectyp = MAG;
        Dattyp = REAL_8;
        putDbl(mag, 0);
        wrstrm();
    }

    // WRITE ANGLE
    if (angle != 0.0)
    {
        Length = 8;
        Rectyp = ANGLE;
        Dattyp = REAL_8;
        putDbl(angle, 0);
        wrstrm();
    }

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;
    if (x >= 0.0) { putI32((long int)((x + epsilon) * factor), 0); }
    else          { putI32((long int)((x - epsilon) * factor), 0); }
    if (y >= 0.0) { putI32((long int)((y + epsilon) * factor), 4); }
    else          { putI32((long int)((y - epsilon) * factor), 4); }
    wrstrm();

    // WRITE STRING
    strcpy(Record, txt);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRING;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE ENDEL
    endEl();
}


//*********** putText *********************
//* PUTS A TEXT RECORD ON SPECIFIED LAYER *
//* layer must be >=0 && <=64K            *
//*****************************************
//#  <text>::=          TEXT [ELFLAGS] [PLEX] LAYER <textbody>                   #
//#                                                                              #
//#  <textbody>::=      TEXTTYPE [PRESENTATION] [PATHTYPE] [WIDTH] [<strans>] XY #
//#                     STRING                                                   #
//#                                                                              #
//#  <strans>::=        STRANS [MAG] [ANGLE]                                     #
void GDSFILE::putText(unsigned short layer,
                      unsigned short textType,
                      unsigned short fontType,
                      char*          textJust,
                      unsigned short pathType,
                      double         width,
                      unsigned short ref,
                      double         mag,
                      double         angle,
                      double         x,
                      double         y,
                      char*          txt)
{
    unsigned short int presentation = 0;
    double   fstep;
    long int istep;

    // WRITE TEXT
    Length = 0;
    Rectyp = TEXT;
    Dattyp = NO_DATA;
    wrstrm();

     // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE TEXTTYPE
    Length = 2;
    Rectyp = TEXTTYPE;
    Dattyp = INTEGER_2;
    putI16(textType, 0);
    wrstrm();

    if (fontType)
    {
        if (fontType == 3)
        {
            presentation |= 0x30;
        }
        else if (fontType == 2)
        {
            presentation |= 0x20;
        }
        else if (fontType == 1)
        {
            presentation |= 0x10;
        }
    }
    if (strlen(textJust) == 2)
    {
        if (strncmp(textJust,"tl",2)) // default
        {
            if (textJust[0] == 'b')
            {
                presentation |= 0x8;
            }
            else if (textJust[0] == 'm')
            {
                presentation |= 0x4;
            }

            if (textJust[1] == 'r')
            {
                presentation |= 0x2;
            }
            else if (textJust[1] == 'c')
            {
                presentation |= 0x1;
            }
        }
    }
    if (presentation)
    {
        // WRITE PRESENTATION
        Length = 2;
        Rectyp = PRESENTATION;
        Dattyp = BIT_ARRAY;             // bit array
        putI16(presentation, 0);
        wrstrm();
    }
    // WRITE PATHTYPE
    Length = 2;
    Rectyp = PATHTYPE;
    Dattyp = INTEGER_2;
    putI16(pathType, 0);
    wrstrm();

    fstep = (double) (width + G_epsilon) * 1000; // done in steps because of previous compilier bug
    istep = (long int) fstep;
    if (istep)
    {
        // WRITE WIDTH
        Length = 4;
        Rectyp = WIDTH;
        Dattyp = INTEGER_4;
        putI32(istep, 0);
        wrstrm();
    }

    // WRITE STRANS
    Length = 2;
    Rectyp = STRANS;
    Dattyp = BIT_ARRAY;             // bit array
    putI16(ref * 0x8000, 0);
    wrstrm();

    // WRITE MAG
    if (mag != 1.0)
    {
        Length = 8;
        Rectyp = MAG;
        Dattyp = REAL_8;
        putDbl(mag, 0);
        wrstrm();
    }

    // WRITE ANGLE
    if (angle != 0.0)
    {
        Length = 8;
        Rectyp = ANGLE;
        Dattyp = REAL_8;
        putDbl(angle, 0);
        wrstrm();
    }

    // WRITE XY
    Length = 8;
    Rectyp = XY;
    Dattyp = INTEGER_4;
    if (x >= 0.0) { putI32((long int)((x + G_epsilon) * 1000), 0); }
    else          { putI32((long int)((x - G_epsilon) * 1000), 0); }
    if (y >= 0.0) { putI32((long int)((y + G_epsilon) * 1000), 4); }
    else          { putI32((long int)((y - G_epsilon) * 1000), 4); }
    wrstrm();

    // WRITE STRING
    strcpy(Record, txt);
    Length = strlen(Record);
    if (Length%2) {
        Record[Length]     = '\0';
        Record[Length + 1] = '\0';
        Length++;
    }
    Rectyp = STRING;
    Dattyp = ACSII_STRING;
    wrstrm();

    // WRITE ENDEL
    endEl();
}


// PUTS AN N-SIDED BOUNDARY ON A SPECIFIED LAYER
int GDSFILE::putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert, double dbu_uu)
{
    int      index;
    double   fstep;
    long int istep;
    double   factor = (1 / dbu_uu);
    double   epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE BOUNDARY
    Length = 0;
    Rectyp = BOUNDARY;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE XY
    Length = nVert * 8;             // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        if (xArray[index] >= 0.0) { fstep = (double) (xArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (xArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8));

        if (yArray[index] >= 0.0) { fstep = (double) (yArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (yArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        if (fstep >= 0.0) { fstep += epsilon; }
        else              { fstep -= epsilon; }
        istep = (long int) fstep;
        putI32(istep, (index*8) + 4);
    }
    if (xArray[0] != xArray[nVert - 1] && yArray[0] != yArray[nVert - 1])
    {
        if (xArray[0] >= 0.0) { fstep = (double) (xArray[0] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (xArray[0] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index * 8));

        if (yArray[0] >= 0.0) { fstep = (double) (yArray[0] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (yArray[0] - epsilon) * factor; } // done in steps because of previous compilier bug
        if (fstep >= 0.0) { fstep += epsilon; }
        else              { fstep -= epsilon; }
        istep = (long int) fstep;
        putI32(istep, (index * 8) + 4);
        Length += 8;
    }
    wrstrm();

    // WRITE ENDEL
    endEl();

    return 0;
}


// PUTS AN N-SIDED BOUNDARY ON A SPECIFIED LAYER
int GDSFILE::putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert)
{
    int      index;
    double   fstep;
    long int istep;

    // WRITE BOUNDARY
    Length = 0;
    Rectyp = BOUNDARY;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE XY
    Length = nVert * 8;             // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        if (xArray[index] >= 0.0) { fstep = (double) (xArray[index] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (xArray[index] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8));

        if (yArray[index] >= 0.0) { fstep = (double) (yArray[index] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (yArray[index] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8) + 4);
    }
    if (xArray[0] != xArray[nVert - 1] && yArray[0] != yArray[nVert - 1])
    {
        if (xArray[0] >= 0.0) { fstep = (double) (xArray[0] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (xArray[0] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index * 8));

        if (yArray[0] >= 0.0) { fstep = (double) (yArray[0] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (yArray[0] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index * 8) + 4);
        Length += 8;
    }
    wrstrm();

    // WRITE ENDEL
    endEl();

    return 0;
}


// PUTS AN N-SIDED BOUNDARY ON A SPECIFIED LAYER
int GDSFILE::putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING], double dbu_uu)
{
    int      index;
    double   fstep;
    long int istep;
    double   factor = (1 / dbu_uu);
    double   epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE BOUNDARY
    Length = 0;
    Rectyp = BOUNDARY;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE XY
    Length = nVert * 8;             // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        if (xArray[index] >= 0.0) { fstep = (double) (xArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (xArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8));

        if (yArray[index] >= 0.0) { fstep = (double) (yArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (yArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8) + 4);
    }
    if (xArray[0] != xArray[nVert - 1] && yArray[0] != yArray[nVert - 1])
    {
        if (xArray[0] >= 0.0) { fstep = (double) (xArray[0] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (xArray[0] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index * 8));

        if (yArray[0] >= 0.0) { fstep = (double) (yArray[0] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (yArray[0] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index * 8) + 4);
        Length += 8;
    }
    wrstrm();

    for (index=0; index <= propIndex; index++)
    {
        Length = 2;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPATTR;
        Dattyp = INTEGER_2;
        putI16(propNumArray[index], 0);
        wrstrm();

        Length = 4;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPVALUE;
        Dattyp = ACSII_STRING;
        strcpy(Record, propValueArray[index]);
        Length = strlen(Record);
        if (Length%2) {
            Record[Length]     = '\0';
            Record[Length + 1] = '\0';
            Length++;
        }
        wrstrm();
    }

    // WRITE ENDEL
    endEl();

    return 0;
}

// PUTS AN N-SIDED BOUNDARY ON A SPECIFIED LAYER
int GDSFILE::putBndDbl(int layer, int datatyp, double xArray[], double yArray[], int nVert, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING])
{
    int      index;
    double   fstep;
    long int istep;

    // WRITE BOUNDARY
    Length = 0;
    Rectyp = BOUNDARY;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE XY
    Length = nVert * 8;             // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        if (xArray[index] >= 0.0) { fstep = (double) (xArray[index] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (xArray[index] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8));

        if (yArray[index] >= 0.0) { fstep = (double) (yArray[index] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                      { fstep = (double) (yArray[index] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8) + 4);
    }
    if (xArray[0] != xArray[nVert - 1] && yArray[0] != yArray[nVert - 1])
    {
        if (xArray[0] >= 0.0) { fstep = (double) (xArray[0] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (xArray[0] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index * 8));

        if (yArray[0] >= 0.0) { fstep = (double) (yArray[0] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                  { fstep = (double) (yArray[0] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index * 8) + 4);
        Length += 8;
    }
    wrstrm();

    for (index=0; index <= propIndex; index++)
    {
        Length = 2;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPATTR;
        Dattyp = INTEGER_2;
        putI16(propNumArray[index], 0);
        wrstrm();

        Length = 4;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPVALUE;
        Dattyp = ACSII_STRING;
        strcpy(Record, propValueArray[index]);
        Length = strlen(Record);
        if (Length%2) {
            Record[Length]     = '\0';
            Record[Length + 1] = '\0';
            Length++;
        }
        wrstrm();
    }

    // WRITE ENDEL
    endEl();

    return 0;
}


// PUTS AN N-SIDED BOUNDARY ON A SPECIFIED LAYER --- note: 1.0 = 1000
int GDSFILE::putBndInt(int layer, int datatyp, int xArray[], int yArray[], int nVert)
{
    int     index;

    // WRITE BOUNDARY
    Length = 0;
    Rectyp = BOUNDARY;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE XY
    Length = nVert * 8;             // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        putI32(xArray[index], (index*8));
        putI32(yArray[index], (index*8) + 4);
    }
    if (xArray[0] != xArray[nVert - 1] && yArray[0] != yArray[nVert - 1])
    {
        putI32(xArray[0], (index*8));
        putI32(yArray[0], (index*8) + 4);
        Length += 8;
    }
    wrstrm();

    // WRITE ENDEL
    endEl();
    return 0;
}

// PUTS AN N-point PATH ON A SPECIFIED LAYER --- note: 1.0 = 1000 */
// <path>::= PATH [ELFLAGS] [PLEX] LAYER DATATYPE [PATHTYPE] [WIDTH] [BGNEXTN] [ENDEXTN] [XY]
int GDSFILE::putPathDbl(int layer, int datatyp, int pathtyp, double width, double bgnextn, double endextn, double xArray[], double yArray[], int nVert, int propIndex, int propNumArray[], char propValueArray[][LENGTHLSTRING], double dbu_uu)
{
    int      index;
    double   fstep;
    long int istep;
    double   factor = (1 / dbu_uu);
    double   epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE PATH
    Length = 0;
    Rectyp = PATH;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    if (pathtyp > 0)
    {
        // WRITE PATHTYPE
        Length = 2;
        Rectyp = PATHTYPE;
        Dattyp = INTEGER_2;
        putI16(pathtyp, 0);
        wrstrm();
    }

    // WRITE WIDTH
    Length = 4;
    Rectyp = WIDTH;
    Dattyp = INTEGER_4;
    fstep = (double) (width + epsilon) * factor; // done in steps because of previous compilier bug
    istep = (long int) fstep;
    putI32(istep, 0);
    wrstrm();

    if (pathtyp == 4)
    {
        // WRITE BGNEXTN
        Length = 4;
        Rectyp = BGNEXTN;
        Dattyp = INTEGER_4;
        if (bgnextn >= 0.0) { fstep =  (double) (bgnextn + epsilon) * factor; } // done in steps because of previous compilier bug
        else                { fstep =  (double) (bgnextn - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, 0);
        wrstrm();

        // WRITE ENDEXTN
        Length = 4;
        Rectyp = ENDEXTN;
        Dattyp = INTEGER_4;
        if (endextn >= 0.0) { fstep =  (double) (endextn + epsilon) * factor; } // done in steps because of previous compilier bug
        else                { fstep =  (double) (endextn - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, 0);
        wrstrm();
    }

    // WRITE XY
    Length = nVert * 8;       // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        if (xArray[index] >= 0.0) { fstep =  (double) (xArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep =  (double) (xArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8));

        if (yArray[index] >= 0.0) { fstep =  (double) (yArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep =  (double) (yArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8) + 4);
    }
    wrstrm();

    for (index=0; index <= propIndex; index++)
    {
        Length = 2;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPATTR;
        Dattyp = INTEGER_2;
        putI16(propNumArray[index], 0);
        wrstrm();

        Length = 4;             // N, 4 byte records of xArray, and yArray
        Rectyp = PROPVALUE;
        Dattyp = ACSII_STRING;
        strcpy(Record, propValueArray[index]);
        Length = strlen(Record);
        if (Length%2) {
            Record[Length]     = '\0';
            Record[Length + 1] = '\0';
            Length++;
        }
        wrstrm();
    }

    // WRITE ENDEL
    endEl();

    return 0;
}


// PUTS AN N-point PATH ON A SPECIFIED LAYER --- note: 1.0 = 1000 */
// <path>::= PATH [ELFLAGS] [PLEX] LAYER DATATYPE [PATHTYPE] [WIDTH] [BGNEXTN] [ENDEXTN] [XY]
int GDSFILE::putPathDbl(int layer, int datatyp, int pathtyp, double width, double bgnextn, double endextn, double xArray[], double yArray[], int nVert, double dbu_uu)
{
    int      index;
    double   fstep;
    long int istep;
    double   factor = (1 / dbu_uu);
    double   epsilon = (dbu_uu / 20.0);
    if (G_epsilon < epsilon) epsilon = G_epsilon;

    // WRITE PATH
    Length = 0;
    Rectyp = PATH;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    if (pathtyp > 0)
    {
        // WRITE PATHTYPE
        Length = 2;
        Rectyp = PATHTYPE;
        Dattyp = INTEGER_2;
        putI16(pathtyp, 0);
        wrstrm();
    }

    // WRITE WIDTH
    Length = 4;
    Rectyp = WIDTH;
    Dattyp = INTEGER_4;
    fstep = (double) (width + epsilon) * factor; // done in steps because of previous compilier bug
    istep = (long int) fstep;
    putI32(istep, 0);
    wrstrm();

    if (pathtyp == 4)
    {
        // WRITE BGNEXTN
        Length = 4;
        Rectyp = BGNEXTN;
        Dattyp = INTEGER_4;
        if (bgnextn >= 0.0) { fstep =  (double) (bgnextn + epsilon) * factor; } // done in steps because of previous compilier bug
        else                { fstep =  (double) (bgnextn - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, 0);
        wrstrm();

        // WRITE ENDEXTN
        Length = 4;
        Rectyp = ENDEXTN;
        Dattyp = INTEGER_4;
        if (endextn >= 0.0) { fstep =  (double) (endextn + epsilon) * factor; } // done in steps because of previous compilier bug
        else                { fstep =  (double) (endextn - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, 0);
        wrstrm();
    }

    // WRITE XY
    Length = nVert * 8;       // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        if (xArray[index] >= 0.0) { fstep =  (double) (xArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep =  (double) (xArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8));

        if (yArray[index] >= 0.0) { fstep =  (double) (yArray[index] + epsilon) * factor; } // done in steps because of previous compilier bug
        else                      { fstep =  (double) (yArray[index] - epsilon) * factor; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8) + 4);
    }
    wrstrm();

    // WRITE ENDEL
    endEl();

    return 0;
}


// PUTS AN N-point PATH ON A SPECIFIED LAYER --- note: 1.0 = 1000 */
// <path>::= PATH [ELFLAGS] [PLEX] LAYER DATATYPE [PATHTYPE] [WIDTH] [BGNEXTN] [ENDEXTN] [XY]
int GDSFILE::putPathDbl(int layer, int datatyp, int pathtyp, double width, double bgnextn, double endextn, double xArray[], double yArray[], int nVert)
{
    int      index;
    double   fstep;
    long int istep;

    // WRITE PATH
    Length = 0;
    Rectyp = PATH;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    if (pathtyp > 0)
    {
        // WRITE PATHTYPE
        Length = 2;
        Rectyp = PATHTYPE;
        Dattyp = INTEGER_2;
        putI16(pathtyp, 0);
        wrstrm();
    }

    // WRITE WIDTH
    Length = 4;
    Rectyp = WIDTH;
    Dattyp = INTEGER_4;
    fstep = (double) (width + G_epsilon) * 1000; // done in steps because of previous compilier bug
    istep = (long int) fstep;
    putI32(istep, 0);
    wrstrm();

    if (pathtyp == 4)
    {
        // WRITE BGNEXTN
        Length = 4;
        Rectyp = BGNEXTN;
        Dattyp = INTEGER_4;
        if (bgnextn >= 0.0) { fstep =  (double) (bgnextn + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                { fstep =  (double) (bgnextn - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, 0);
        wrstrm();

        // WRITE ENDEXTN
        Length = 4;
        Rectyp = ENDEXTN;
        Dattyp = INTEGER_4;
        if (endextn >= 0.0) { fstep =  (double) (endextn + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                { fstep =  (double) (endextn - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, 0);
        wrstrm();
    }

    // WRITE XY
    Length = nVert * 8;       // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++)
    {
        if (xArray[index] >= 0.0) { fstep =  (double) (xArray[index] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                      { fstep =  (double) (xArray[index] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8));

        if (yArray[index] >= 0.0) { fstep =  (double) (yArray[index] + G_epsilon) * 1000; } // done in steps because of previous compilier bug
        else                      { fstep =  (double) (yArray[index] - G_epsilon) * 1000; } // done in steps because of previous compilier bug
        istep = (long int) fstep;
        putI32(istep, (index*8) + 4);
    }
    wrstrm();

    // WRITE ENDEL
    endEl();

    return 0;
}


// PUTS AN N-point PATH ON A SPECIFIED LAYER --- note: 1.0 = 1000 */
int GDSFILE::putPathInt(int layer, int datatyp, int width, int xArray[], int yArray[], int nVert)
{
    int  index;

    // WRITE PATH
    Length = 0;
    Rectyp = PATH;
    Dattyp = NO_DATA;
    wrstrm();

    // WRITE LAYER
    Length = 2;
    Rectyp = LAYER;
    Dattyp = INTEGER_2;
    putI16(layer, 0);
    wrstrm();

    // WRITE DATATYPE
    Length = 2;
    Rectyp = DATATYPE;
    Dattyp = INTEGER_2;
    putI16(datatyp, 0);
    wrstrm();

    // WRITE WIDTH
    Length = 4;
    Rectyp = WIDTH;
    Dattyp = INTEGER_4;
    putI32(width, 0);
    wrstrm();

    // WRITE XY
    Length = nVert * 8;       // N, 4 byte records of xArray, and yArray
    Rectyp = XY;
    Dattyp = INTEGER_4;
    for (index=0; index < nVert; index++) {
        putI32(xArray[index], (index*8));
        putI32(yArray[index], (index*8) + 4);
    }
    wrstrm();

    // WRITE ENDEL
    endEl();

    return 0;
}


void GDSFILE::copyRecord(char* copy)    //copy current record to "copy"
{
    char* ptr = record();
    for(int i=0; i<204800; i++) copy[i] = *ptr++;
}


////////////////////////////////////////////////
int GDSFILE::roundInt(int input, int round2grid)
{
    if (input > 0 ) return(((input + (round2grid/2))/round2grid ) * round2grid);
    else            return(((input - (round2grid/2))/round2grid ) * round2grid);
}

