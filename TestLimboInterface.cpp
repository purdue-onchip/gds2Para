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
#include <parser-spef/parser-spef.hpp>
#include "limboint.h"
#include "spefwrite.h"
using std::cerr;
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
    if (argc == 2)
    {
        if (strcmp(argv[1], "--help") == 0)
        {
            cout << "Help for Test Limbo Interface binary" << endl;
            cout << "Usage: Test_$@ [options] file1 [file2]" << endl;
            cout << "Options:" << endl;
            cout << "  --help                Display this information." << endl;
            cout << "  --version             Print the version number." << endl;
            cout << "  -r, --read            Read given GDSII file into memory." << endl;
            cout << "  -w, --write           Write database in memory to given file." << endl;
            cout << "  -i, --imap            Read given IMAP 3D file and write GDSII file with name also given." << endl;
            cout << endl << "Comments:" << endl;
            cout << "The file passed after -r or --read must be Calma GDSII stream file." << endl;
            cout << " The file passed after -w or --write must be a blank SPEF file." << endl;
            cout << " The first file passed after -i or --imap must be a 3D description .imp file, and the second must be a blank .gds file." << endl;
            cout << endl << "Bug reporting:" << endl;
            cout << "Visit <https://github.com/purdue-onchip/gdsii-interface>" << endl;
        }
        else if (strcmp(argv[1], "--version") == 0)
        {
            cout << "Version Number for Test Limbo Interface binary: " << "1.0" << endl;
        }
        else
        {
            cerr << "Only \"--help\" or \"--version\" flags allowed without files passed" << endl;
        }
    }
    else if (argc == 3)
    {
        if ((strcmp(argv[1], "-r") == 0) || (strcmp(argv[1], "--read") == 0))
        {
            AsciiDataBase adb;
            string fName = argv[2];
            adb.setFileName(fName);
            std::ifstream inFile(fName.c_str());
            GdsParser::GdsReader adbReader(adb);
            bool adbIsGood = adbReader(inFile);
            vector<size_t> indCellPrint = { adb.getNumCell() - 1 };
            adb.print(indCellPrint);

            /*EnumDataBase edb;
            GdsParser::GdsReader edbReader(edb);
            bool edbIsGood = edbReader(inFile);
            cout << "Test Enum API: " << edbIsGood << endl;
            //cout << "Test Enum API: " << GdsParser::read(edb, argv[1]) << endl;*/
        }
        else if ((strcmp(argv[1], "-w") == 0) || (strcmp(argv[1], "--write") == 0))
        {
            // Load sample information in memory
            string design = "test_out";
            Waveforms blank;
            size_t numPorts = 13;
            vector<std::string> ports = { "inp1", "u1:a", "inp2", "u1:b", "out", "u3:o", "u1:o", "u4:a", "u4:o", "f1:d", "f1:a", "u2:a", "u4:b" };
            vector<char> portDir = { 'I', 'I', 'I', 'I', 'O', 'O', 'O', 'I', 'O', 'I', 'O', 'I', 'I' };
            vector<vector<double>> matG = { { 1. / (10.5e3), -1. / (10.5e3), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },{ -1. / (10.5e3), 1. / (10.5e3), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
            { 0., 0., 1. / (4.5e3), -1. / (4.5e3), 0., 0., 0., 0., 0., 0., 0., 0., 0. },{ 0., 0., -1. / (4.5e3), 1. / (4.5e3), 0., 0., 0., 0., 0., 0., 0., 0., 0. },
            { 0., 0., 0., 0., 1. / (1.4e3), -1. / (1.4e3), 0., 0., 0., 0., 0., 0., 0. },{ 0., 0., 0., 0., -1. / (1.4e3), 1. / (1.4e3), 0., 0., 0., 0., 0., 0., 0. },
            { 0., 0., 0., 0., 0., 0., 1. / (2.1e3), -1. / (2.1e3), 0., 0., 0., 0., 0. },{ 0., 0., 0., 0., 0., 0., -1. / (2.1e3), 1. / (2.1e3), 0., 0., 0., 0., 0. },
            { 0., 0., 0., 0., 0., 0., 0., 0., 1. / (2.1e3), -1. / (2.1e3), 0., 0., 0. },{ 0., 0., 0., 0., 0., 0., 0., 0., -1. / (2.1e3), 1. / (2.1e3), 0., 0., 0. },
            { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1. / (1.2e3 + 2.3e3 + 3.4e3) + 1. / (1.2e3 + 7.8e3 + 5.6e3), -1. / (1.2e3 + 2.3e3 + 3.4e3), -1. / (1.2e3 + 7.8e3 + 5.6e3) },
            { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1. / (3.4e3 + 2.3e3 + 1.2e3), 1. / (3.4e3 + 2.3e3 + 1.2e3) + 1. / (3.4e3 + 2.3e3 + 7.8e3 + 5.6e3), -1. / (3.4e3 + 2.3e3 + 7.8e3 + 5.6e3) },
            { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1. / (5.6e3 + 7.8e3 + 1.2e3), -1. / (5.6e3 + 7.8e3 + 2.3e3 + 3.4e3), 1. / (5.6e3 + 7.8e3 + 1.2e3) + 1. / (5.6e3 + 7.8e3 + 2.3e3 + 3.4e3) } };
            vector<vector<double>> matC = { { 2.5e-15, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },{ 0., 2.9e-5, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
            { 0., 0., 0.7e-15, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },{ 0., 0., 0., 1.3e-15, 0., 0., 0., 0., 0., 0., 0., 0., 0. },{ 0., 0., 0., 0., 0.5e-15, 0., 0., 0., 0., 0., 0., 0., 0. },
            { 0., 0., 0., 0., 0., 0.2e-15, 0., 0., 0., 0., 0., 0., 0. },{ 0., 0., 0., 0., 0., 0., 0.35e-15, 0., 0., 0., 0., 0., 0. },{ 0., 0., 0., 0., 0., 0., 0., 0.65e-15, 0., 0., 0., 0., 0. },
            { 0., 0., 0., 0., 0., 0., 0., 0., 0.7e-15, 0., 0., 0., 0. },{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5e-15, 0., 0., 0. },{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 8.9e-15, 0., 0. },
            { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 6.7e-15, 0. },{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 7.8e-15 } };
            Parasitics sample(numPorts, ports, portDir, matG, matC);
            SolverDataBase sdb(design, blank, sample);

            // Prepare to write to file
            string fName = argv[2];
            sdb.setOutSPEF(fName);
            std::ifstream inFile(fName.c_str());
            //GdsParser::GdsReader adbReader(sdb);
            //bool sdbIsGood = adbReader(inFile);
            bool couldDump = sdb.printDump();
        }
        else
        {
            cerr << "Must pass a file after \"-r\" or \"-w\" flags, rerun with \"--help\" flag for details" << endl;
        }
    }
    else if (argc == 4)
    {
        if ((strcmp(argv[1], "-i") == 0) || (strcmp(argv[1], "--imap") == 0))
        {
            // Read IMAP 3D description file and write to GDSII file
            SolverDataBase sdb;
            string inIMAPFile = argv[2];
            string outGDSIIFile = argv[3];
            bool sdbIsGood = sdb.readIMAPwriteGDSII(inIMAPFile, outGDSIIFile);
            cout << "File ready at " << outGDSIIFile << endl;
        }
        else
        {
            cerr << "Must pass a IMAP file to read and blank GDSII file to write after \"-i\" flags, rerun with \"--help\" flag for details" << endl;
        }
    }
    else
    {
        cerr << "Either 1, 2, or 3 arguments are required, use \"--help\" flag for details" << endl;
    }

    return 0;
}
