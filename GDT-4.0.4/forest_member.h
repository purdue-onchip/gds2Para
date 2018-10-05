// Copyright 1995-2010 by Ken Schumack (Schumack@cpan.org
// @(#) $Id: forest_member.h 67 2010-09-24 20:14:19Z schumack $ 
#include <stdio.h>
int forest_member(char* string)
{
    char c;
    
    c = string[0];
    
    switch (c)
    {
        case '0':
            return 0;
            break;
        case '1':
            return 1;
            break;
        case '2':
            return 2;
            break;
        case '3':
            return 3;
            break;
        case '4':
            return 4;
            break;
        case '5':
            return 5;
            break;
        case '6':
            return 6;
            break;
        case '7':
            return 7;
            break;
        case '8':
            return 8;
            break;
        case '9':
            return 9;
            break;
        case 'A':
        case 'a':
            return 10;
            break;
        case 'B':
        case 'b':
            return 11;
            break;
        case 'C':
        case 'c':
            return 12;
            break;
        case 'D':
        case 'd':
            return 13;
            break;
        case 'E':
        case 'e':
            return 14;
            break;
        case 'F':
        case 'f':
            return 15;
            break;
        case 'G':
        case 'g':
            return 16;
            break;
        case 'H':
        case 'h':
            return 17;
            break;
        case 'I':
        case 'i':
            return 18;
            break;
        case 'J':
        case 'j':
            return 19;
            break;
        case 'K':
        case 'k':
            return 20;
            break;
        case 'L':
        case 'l':
            return 21;
            break;
        case 'M':
        case 'm':
            return 22;
            break;
        case 'N':
        case 'n':
            return 23;
            break;
        case 'O':
        case 'o':
            return 24;
            break;
        case 'P':
        case 'p':
            return 25;
            break;
        case 'Q':
        case 'q':
            return 26;
            break;
        case 'R':
        case 'r':
            return 27;
            break;
        case 'S':
        case 's':
            return 28;
            break;
        case 'T':
        case 't':
            return 29;
            break;
        case 'U':
        case 'u':
            return 30;
            break;
        case 'V':
        case 'v':
            return 31;
            break;
        case 'W':
        case 'w':
            return 32;
            break;
        case 'X':
        case 'x':
            return 33;
            break;
        case 'Y':
        case 'y':
            return 34;
            break;
        case 'Z':
        case 'z':
            return 35;
            break;
        default:
            return 36;
            break;
    }
}

