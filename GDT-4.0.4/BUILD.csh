#!/bin/tcsh -f
# @(#) $Id: BUILD.csh 95 2014-10-28 23:03:27Z schumack $

unsetenv LD_LIBRARY_PATH
set DEBUG=0

set option=""
if ($# > 0) then
  set option=$1
endif

if ("$option" == "-debug") then
  set DEBUG=1
  echo "Building in debug mode. DEBUG=$DEBUG"
else
  echo "Building in release mode. DEBUG=$DEBUG"
endif

set path=($path /bin /usr/bin /usr/sbin /sbin) ## for getconf uname rm echo gcc
set OS=`uname -s`
rm gds2gdt.$OS >& /dev/null

#-----------------------------------------------------------------
set conf=`getconf LFS_CFLAGS`
if ("$conf" == "*-D*") then
  set CCFLAGS="$conf"
else
  set CCFLAGS="-D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64"
endif

set EXTRA_LIBS="-L."
set LIBDIR=""
if (`uname -m` == "x86_64") then #for x86_64
  set LIBDIRS=(/usr/lib64 /usr/lib/x86_64-linux-gnu)
else
  set LIBDIRS=(/usr/lib /lib/i386-linux-gnu)
endif
foreach dirName ($LIBDIRS)
  echo "checking if $dirName is valid"
  if ("$LIBDIR" == "") then
    if ((-e "$dirName/libm.a") || (-e "$dirName/libm.so")) then
      set LIBDIR=$dirName
    endif
  endif
  if (-d $dirName) then
    set EXTRA_LIBS="$EXTRA_LIBS -L$dirName"
  endif
end

if ("$OS" == "SunOS") then ## some solaris machines need this
  set CC="gcc -m32"
else
  set CC="gcc"
endif

set EXTRA_INCLUDES="-I."

set CC1=`echo $CC|sed 's| .*||'`
if ($DEBUG) then
  echo CC=$CC
  if ("$CC1" == 'gcc' ) then
    echo gcc --version == `gcc --version |head -1`
  endif
endif
if (( $CC1 == "gcc" ) && ( `gcc --version |head -1` =~ '* 4.*' )) then
  #note: not all versions of sed have -i option so use .tmp file
  foreach f (`ls *.C`)
    echo uncommenting gcc4 stuff in $f
    sed 's|^//#include <gcc4.h|#include <gcc4.h|' $f > $f.tmp
    sed 's|^//using namespace|using namespace|' $f.tmp > $f
    sed 's|^#include <iostream.h|//#include <iostream.h|' $f > $f.tmp
    sed 's|^#include <ostream.h|//#include <ostream.h|' $f.tmp > $f
    #mv $f.tmp $f
    rm $f.tmp
  end
else
  foreach f (`ls *.C`)
    echo commenting gcc4 stuff in $f
    sed 's|^#include <gcc4.h|//#include <gcc4.h|' $f > $f.tmp
    sed 's|^using namespace|//using namespace|' $f.tmp > $f
    sed 's|^//#include <iostream.h|#include <iostream.h|' $f > $f.tmp
    sed 's|^//#include <ostream.h|#include <ostream.h|' $f.tmp > $f
    #mv $f.tmp $f
    rm $f.tmp
  end
endif

set OPT='-O'
if ($DEBUG) then
  set OPT='-g -DDEBUG' ##'-g -pg' for debugging
endif

rm *.o >& /dev/null
echo ""

#------------------------------- gds2gdt ----------------------------------
echo "Compiling gds2gdt --------------------------------------------------"
set mathLibStatic="-lm $LIBDIR/libm.a -static" #prefer static version
if (! -e "$LIBDIR/libm.a") then
  echo "Did not find static math lib. Will try to use dynamically linked version"
  set mathLibStatic="$LIBDIR/libm.so"
endif

echo $CC $CCFLAGS $EXTRA_INCLUDES sRemoveTrailingZeros ====
$CC $CCFLAGS $EXTRA_INCLUDES sRemoveTrailingZeros.C -c $OPT -o sRemoveTrailingZeros.o

echo $CC $CCFLAGS $EXTRA_INCLUDES get_field =======
$CC $CCFLAGS $EXTRA_INCLUDES get_field.C -c $OPT -o get_field.o

echo $CC $CCFLAGS $EXTRA_INCLUDES stoupper =======
$CC $CCFLAGS $EXTRA_INCLUDES stoupper.C -c $OPT -o stoupper.o

echo $CC $CCFLAGS $EXTRA_INCLUDES sfind =======
$CC $CCFLAGS $EXTRA_INCLUDES sfind.C -c $OPT -o sfind.o

echo $CC $CCFLAGS $EXTRA_INCLUDES match_string =======
$CC $CCFLAGS $EXTRA_INCLUDES match_string.C -c $OPT -o match_string.o

echo $CC $CCFLAGS $EXTRA_INCLUDES sRemoveSpaces =======
$CC $CCFLAGS $EXTRA_INCLUDES sRemoveSpaces.C -c $OPT -o sRemoveSpaces.o

echo $CC $CCFLAGS $EXTRA_INCLUDES sRemoveWhiteSpace ====
$CC $CCFLAGS $EXTRA_INCLUDES sRemoveWhiteSpace.C -c $OPT -o sRemoveWhiteSpace.o

echo $CC $CCFLAGS $EXTRA_INCLUDES mystrncpy =======
$CC $CCFLAGS $EXTRA_INCLUDES mystrncpy.C -c $OPT -o mystrncpy.o

echo $CC $CCFLAGS $EXTRA_INCLUDES gdsStream =======
$CC $CCFLAGS $EXTRA_INCLUDES gdsStream.C -c $OPT -o gdsStream.o -Wno-deprecated

echo $CC $CCFLAGS $EXTRA_INCLUDES $EXTRA_LIBS gds2gdt =======
echo $CC gds2gdt.C $OPT -o gds2gdt.$OS $EXTRA_INCLUDES $EXTRA_LIBS -lstdc++ *.o -Wreturn-type -Wswitch -Wcomment -Wformat -Wchar-subscripts -Wparentheses -Wpointer-arith -Wcast-qual -Woverloaded-virtual -Wno-write-strings $mathLibStatic -Wno-deprecated
$CC gds2gdt.C $OPT -o gds2gdt.$OS $EXTRA_INCLUDES $EXTRA_LIBS -lstdc++ *.o -Wreturn-type -Wswitch -Wcomment -Wformat -Wchar-subscripts -Wparentheses -Wpointer-arith -Wcast-qual -Woverloaded-virtual -Wno-write-strings $mathLibStatic -Wno-deprecated

if (! $DEBUG) then
  strip gds2gdt.$OS
endif
rm *.o >& /dev/null
echo created gds2gdt.$OS
echo ""

#------------------------------- gdt2gds ----------------------------------
echo "Compiling gdt2gds --------------------------------------------------"
echo $CC $CCFLAGS $EXTRA_INCLUDES sRemoveWhiteSpace ====
$CC $CCFLAGS $EXTRA_INCLUDES sRemoveWhiteSpace.C -c $OPT -o sRemoveWhiteSpace.o

echo $CC $CCFLAGS $EXTRA_INCLUDES sRemoveTrailingZeros ====
$CC $CCFLAGS $EXTRA_INCLUDES sRemoveTrailingZeros.C -c $OPT -o sRemoveTrailingZeros.o

echo $CC $CCFLAGS $EXTRA_INCLUDES get_field =======
$CC $CCFLAGS $EXTRA_INCLUDES get_field.C -c $OPT -o get_field.o

echo $CC $CCFLAGS $EXTRA_INCLUDES stoupper =======
$CC $CCFLAGS $EXTRA_INCLUDES stoupper.C -c $OPT -o stoupper.o

echo $CC $CCFLAGS $EXTRA_INCLUDES sfind =======
$CC $CCFLAGS $EXTRA_INCLUDES sfind.C -c $OPT -o sfind.o

echo $CC $CCFLAGS $EXTRA_INCLUDES match_string =======
$CC $CCFLAGS $EXTRA_INCLUDES match_string.C -c $OPT -o match_string.o

echo $CC $CCFLAGS $EXTRA_INCLUDES sRemoveSpaces =======
$CC $CCFLAGS $EXTRA_INCLUDES sRemoveSpaces.C -c $OPT -o sRemoveSpaces.o

echo $CC $CCFLAGS $EXTRA_INCLUDES mystrncpy =======
$CC $CCFLAGS $EXTRA_INCLUDES mystrncpy.C -c $OPT -o mystrncpy.o

echo $CC $CCFLAGS $EXTRA_INCLUDES gdsStream =======
$CC $CCFLAGS $EXTRA_INCLUDES gdsStream.C -c $OPT -o gdsStream.o -Wno-deprecated

echo $CC $CCFLAGS $EXTRA_INCLUDES gdt2gds =======
$CC $CCFLAGS $EXTRA_INCLUDES $EXTRA_LIBS gdt2gds.C $OPT -o gdt2gds.$OS *.o -lstdc++ $mathLibStatic -Wreturn-type -Wswitch -Wcomment -Wformat -Wchar-subscripts -Wparentheses -Wpointer-arith -Wcast-qual -Woverloaded-virtual -Wno-deprecated -Wno-write-strings

if (! $DEBUG) then
  strip gdt2gds.$OS
endif
rm *.o >& /dev/null
echo "created gdt2gds.$OS"

