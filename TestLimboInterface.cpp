/**
 * @file   TestLimboInterface.cpp 
 * @author Michael R. Hayashi
 * @date   18 October 2018
 * @brief  Primary Limbo Interface Test Function
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
#include "limboint.h"
using std::cout;
using std::endl;

/* ===========================================
example to read .gds.gz 
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

EnumDataBase edb; 
boost::iostreams::filtering_istream in; 
in.push(boost::iostreams::gzip_decompressor());
in.push(boost::iostreams::file_source(argv[1]));

cout << "test enum api\n" << GdsParser::read(edb, in) << endl;
=========================================== */

/// @brief main function 
/// @param argc number of arguments 
/// @param argv values of arguments 
/// @return 0 if succeed 
int main(int argc, char** argv)
{
	if (argc > 1)
    {
        AsciiDataBase adb;
		string fName = argv[1];
		adb.setFileName(fName);
		std::ifstream inFile(fName.c_str());
		GdsParser::GdsReader adbReader(adb);
		bool adbIsGood = adbReader(inFile);
		vector<size_t> indCellPrint = {adb.getNumCell() - 1};
		adb.print(indCellPrint);

        /*EnumDataBase edb;
		GdsParser::GdsReader edbReader(edb);
		bool edbIsGood = edbReader(inFile);
		cout << "Test Enum API: " << edbIsGood << endl;
		//cout << "Test Enum API: " << GdsParser::read(edb, argv[1]) << endl;*/
    }
	else cout << "At least 1 argument is required" << endl;

	return 0;
}
