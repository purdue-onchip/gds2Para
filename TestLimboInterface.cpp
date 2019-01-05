/**
 * @file   TestLimboInterface.cpp 
 * @author Michael R. Hayashi
 * @date   18 October 2018
 * @brief  Primary Limbo Interface Test Function
 */

#define _USE_MATH_DEFINES // Place before including <cmath> for e, log2(e), log10(e), ln(2), ln(10), pi, pi/2, pi/4, 1/pi, 2/pi, 2/sqrt(pi), sqrt(2), and 1/sqrt(2)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
#include <parser-spef/parser-spef.hpp>
#include <Eigen/Sparse>
#include "limboint.h"
#include "spefwrite.h"
using std::cerr;
using std::cout;
using std::endl;

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
            cout << "Usage: Test_$@ [options] file1 [file2..file3]" << endl;
            cout << "Options:" << endl;
            cout << "  --help                Display this information." << endl;
            cout << "  --version             Print the version number." << endl;
            cout << "  -r, --read            Read given GDSII file into memory." << endl;
            cout << "  -p, --parrot          Immediately output given GDSII file after reading." << endl;
            cout << "  -w, --write           Write database in memory to given SPEF file." << endl;
            cout << "  -i, --imp             Read given interconnect modeling platform file and write GDSII file with name also given." << endl;
            cout << "  -s, --simulate        Read GDSII file and sim input file into memory, simulate, and write solution to SPEF file." << endl;
            cout << endl << "Comments:" << endl;
            cout << "The file passed after -r, --read, -p, or --parrot must be a Calma GDSII stream file." << endl;
            cout << " The file passed after -w or --write must be a blank SPEF file." << endl;
            cout << " The first file passed after -i or --imp must be a 3D description .imp file, and the second must be a blank .gds file." << endl;
            cout << " The first file passed after -s or --simulate must be a Calma GDSII stream file, the second must be a sim_input file, and the third must be a blank SPEF file." << endl;
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
        else if ((strcmp(argv[1], "-p") == 0) || (strcmp(argv[1], "--parrot") == 0))
        {
            // Read and print existing file
            AsciiDataBase adb;
            string fName = argv[2];
            size_t indExtension = fName.find(".", 1);
            adb.setFileName(fName.substr(0, indExtension) + "_parrot" + fName.substr(indExtension, string::npos));
            GdsParser::GdsReader adbReader(adb);
            bool adbIsGood = adbReader(fName.c_str());
            adb.print({ adb.getNumCell() - 1 });

            // Dump to parroted file immediately
            adb.dump();
            cout << "Dumped parroted file" << endl;
        }
        else if ((strcmp(argv[1], "-w") == 0) || (strcmp(argv[1], "--write") == 0))
        {
            // Load sample information in memory
            string design = "test_out";
            Waveforms blank;

            // Setup the port information vector
            vector<Port> ports = {};
            ports.emplace_back(Port("inp1", 'I', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u1:a", 'I', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("inp2", 'I', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u1:b", 'I', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("out", 'O', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u3:o", 'O', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u1:o", 'O', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u4:a", 'I', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u4:o", 'O', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("f1:d", 'I', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("f1:a", 'O', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u2:a", 'I', 50., vector<double>(6, 0)));
            ports.emplace_back(Port("u4:b", 'I', 50., vector<double>(6, 0)));

            // Setup the Eigen sparse conductance matrix
            spMat matG(13, 13); // Initialize sparse conductance matrix (S)
            vector<dTriplet> listG; // Initialize triplet list for conductance matrix
            listG.reserve(3); // Reserve room for 3 nonzero entries in each row
            listG.push_back(dTriplet(0, 0, 1. / (10.5e3)));
            listG.push_back(dTriplet(0, 1, -1. / (10.5e3)));
            listG.push_back(dTriplet(1, 0, -1. / (10.5e3)));
            listG.push_back(dTriplet(1, 1, 1. / (10.5e3)));
            listG.push_back(dTriplet(2, 2, 1. / (4.5e3)));
            listG.push_back(dTriplet(2, 3, -1. / (4.5e3)));
            listG.push_back(dTriplet(3, 2, -1. / (4.5e3)));
            listG.push_back(dTriplet(3, 3, 1. / (4.5e3)));
            listG.push_back(dTriplet(4, 4, 1. / (1.4e3)));
            listG.push_back(dTriplet(4, 5, -1. / (1.4e3)));
            listG.push_back(dTriplet(5, 4, -1. / (1.4e3)));
            listG.push_back(dTriplet(5, 5, 1. / (1.4e3)));
            listG.push_back(dTriplet(6, 6, 1. / (2.1e3)));
            listG.push_back(dTriplet(6, 7, -1. / (2.1e3)));
            listG.push_back(dTriplet(7, 6, -1. / (2.1e3)));
            listG.push_back(dTriplet(7, 7, 1. / (2.1e3)));
            listG.push_back(dTriplet(8, 8, 1. / (2.1e3)));
            listG.push_back(dTriplet(8, 9, -1. / (2.1e3)));
            listG.push_back(dTriplet(9, 8, -1. / (2.1e3)));
            listG.push_back(dTriplet(9, 9, 1. / (2.1e3)));
            listG.push_back(dTriplet(10, 10, 1. / (1.2e3 + 2.3e3 + 3.4e3) + 1. / (1.2e3 + 7.8e3 + 5.6e3)));
            listG.push_back(dTriplet(10, 11, -1. / (1.2e3 + 2.3e3 + 3.4e3)));
            listG.push_back(dTriplet(10, 12, -1. / (1.2e3 + 7.8e3 + 5.6e3)));
            listG.push_back(dTriplet(11, 10, -1. / (3.4e3 + 2.3e3 + 1.2e3)));
            listG.push_back(dTriplet(11, 11, 1. / (3.4e3 + 2.3e3 + 1.2e3) + 1. / (3.4e3 + 2.3e3 + 7.8e3 + 5.6e3)));
            listG.push_back(dTriplet(11, 12, -1. / (3.4e3 + 2.3e3 + 7.8e3 + 5.6e3)));
            listG.push_back(dTriplet(12, 10, -1. / (5.6e3 + 7.8e3 + 1.2e3)));
            listG.push_back(dTriplet(12, 11, -1. / (5.6e3 + 7.8e3 + 2.3e3 + 3.4e3)));
            listG.push_back(dTriplet(12, 12, 1. / (5.6e3 + 7.8e3 + 1.2e3) + 1. / (5.6e3 + 7.8e3 + 2.3e3 + 3.4e3)));
            matG.setFromTriplets(listG.begin(), listG.end()); // Assign nonzero entries to sparse conductance matrix
            matG.makeCompressed(); // Conductance matrix in compressed sparse row (CSR) format

            // Setup the Eigen sparse capacitance matrix
            spMat matC(13, 13); // Initialize sparse capacitance matrix (F)
            vector<dTriplet> listC; // Initialize triplet list for capacitance matrix
            listC.reserve(1); // Reserve room for 1 nonzero entry in each row
            listC.push_back(dTriplet(0, 0, 2.5e-15));
            listC.push_back(dTriplet(1, 1, 2.9e-15));
            listC.push_back(dTriplet(2, 2, 0.7e-15));
            listC.push_back(dTriplet(3, 3, 1.3e-15));
            listC.push_back(dTriplet(4, 4, 0.5e-15));
            listC.push_back(dTriplet(5, 5, 0.2e-15));
            listC.push_back(dTriplet(6, 6, 0.35e-15));
            listC.push_back(dTriplet(7, 7, 0.65e-15));
            listC.push_back(dTriplet(8, 8, 0.7e-15));
            listC.push_back(dTriplet(9, 9, 0.5e-15));
            listC.push_back(dTriplet(10, 10, 8.9e-15));
            listC.push_back(dTriplet(11, 11, 6.7e-15));
            listC.push_back(dTriplet(12, 12, 7.8e-15));
            matC.setFromTriplets(listC.begin(), listC.end()); // Assign nonzero entries to sparse capacitance matrix
            matC.makeCompressed(); // Capactiance matrix in compressed sparse row (CSR) format

            // Create variables of custom classes
            Parasitics sample(ports, matG, matC);
            SolverDataBase sdb(design, blank, sample);

            // Prepare to write to file
            string fName = argv[2];
            sdb.setOutSPEF(fName);
            bool couldDump = sdb.printDump();
        }
        else
        {
            cerr << "Must pass a file after \"-r\", \"-p\", or \"-w\" flags, rerun with \"--help\" flag for details" << endl;
        }
    }
    else if (argc == 4)
    {
        if ((strcmp(argv[1], "-i") == 0) || (strcmp(argv[1], "--imp") == 0))
        {
            // Read interconnect modeling platform (IMP) file and write to GDSII file
            SolverDataBase sdb;
            string inIMPFile = argv[2];
            string outGDSIIFile = argv[3];
            bool sdbIsGood = sdb.readIMPwriteGDSII(inIMPFile, outGDSIIFile);
            cout << "File ready at " << outGDSIIFile << endl;
        }
        else
        {
            cerr << "Must pass a .imp file to read and blank GDSII file to write after \"-i\" flag, rerun with \"--help\" flag for details" << endl;
        }
    }
    else if (argc == 5)
    {
        if ((strcmp(argv[1], "-s") == 0) || (strcmp(argv[1], "--simulate") == 0))
        {
            // Read simulation input file
            SolverDataBase sdb;
            string inSimFile = argv[3];
            bool sdbIsGood = sdb.readSimInput(inSimFile);
            cout << "Simulation input file read" << endl;

            // Initialize mesh and other useful variables
            fdtdMesh sys;
            int status;
            clock_t t1 = clock();
            size_t indExtension = inSimFile.find(".", 1);
            string inStackFile = inSimFile.substr(0, indExtension) + "_stack.txt";
            //string inStackFile = inSimFile; // Use stack information contained within sim_input file
            //string inPolyFile = inSimFile.substr(0, indExtension) + "_polygon.txt";

            // Read GDSII file
            AsciiDataBase adb;
            string inGDSIIFile = argv[2];
            adb.setFileName(inGDSIIFile);
            GdsParser::GdsReader adbReader(adb);
            bool adbIsGood = adbReader(inGDSIIFile.c_str());
            vector<size_t> indCellPrint = { adb.getNumCell() - 1 };
            adb.print(indCellPrint, &sys);
            cout << "GDSII file read" << endl;

            // Set the number of input conductors
            sys.numCdtRow = adb.getNumCdtIn();

            // Read the input file
            unordered_map<double, int> xi, yi, zi;
            status = readInput(inStackFile.c_str(), &sys, xi, yi, zi);
            if (status == 0)
                cout << "readInput Success!" << endl;
            else {
                cout << "readInput Fail!" << endl;
                return status;
            }

            // Set D_eps and D_sig
            status = matrixConstruction(&sys);
            if (status == 0) {
                cout << "matrixConstruction Success!" << endl;
            }
            else {
                cout << "matrixConstruction Fail!" << endl;
                return status;
            }

            // Set port
            status = portSet(&sys, xi, yi, zi);
            if (status == 0)
                cout << "portSet Success!\n";
            else {
                cout << "portSet Fail!\n";
                return status;
            }

            // Parameter generation
            status = paraGenerator(&sys, xi, yi, zi);
            if (status == 0)
                cout << "paraGenerator Success!" << endl;
            else {
                cout << "paraGenerator Fail!" << endl;
                return status;
            }
            cout << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;

            // Parameter storage
            // Setup the Eigen sparse conductance matrix
            spMat matG(sys.numPorts, sys.numPorts); // Initialize sparse conductance matrix (S)
            spMat matC(sys.numPorts, sys.numPorts); // Initialize sparse capacitance matrix (F)
            vector<dTriplet> listG; // Initialize triplet list for conductance matrix
            vector<dTriplet> listC;
            listG.reserve(sys.numPorts); // Reserve room so matrix could be dense
            listC.reserve(sys.numPorts);
            for (size_t indi = 0; indi < sys.numPorts; indi++) // Loop over excitation ports
            {
                double sumG = 0.;
                double sumC = 0.;
                for (size_t indj = 0; indj < sys.numPorts; indj++) // Loop over response ports
                {
                    if (indi != indj) // Off-diagonal entries are negated node-to-node admittances
                    {
                        listG.push_back(dTriplet(indi, indj, sys.Y[indi * sys.numPorts + indj].real()));
                        listC.push_back(dTriplet(indi, indj, sys.Y[indi * sys.numPorts + indj].imag() / (2 * M_PI* sys.freqStart * sys.freqUnit)));
                    }
                    
                    // Diagonal entries need to have sum of all attached admittances, so note them
                    sumG += sys.Y[indi * sys.numPorts + indj].real();
                    sumC += sys.Y[indi * sys.numPorts + indj].imag() / (2 * M_PI* sys.freqStart * sys.freqUnit);
                }

                // Record diagonal entries
                listG.push_back(dTriplet(indi, indi, sumG));
                listC.push_back(dTriplet(indi, indi, sumC));
            }
            matG.setFromTriplets(listG.begin(), listG.end()); // Assign nonzero entries to sparse conductance matrix
            matC.setFromTriplets(listC.begin(), listC.end()); // Do not put in compressed sparse row (CSR) format due to density
            Parasitics oldPara = sdb.getParasitics(); // Get outdated parastics structure to update
            sdb.setParasitics(Parasitics(oldPara.getPorts(), matG, matC));

            // Output SPEF file
            string outSPEFFile = argv[4];
            sdb.setDesignName(adb.findNames().back());
            sdb.setOutSPEF(outSPEFFile);
            bool sdbCouldDump = sdb.printDump();
            cout << "File ready at " << outSPEFFile << endl;
        }
        else
        {
            cerr << "Must pass a GDSII file, sim_input file, and blank SPEF file to write after \"-s\" flag, rerun with \"--help\" flag for details" << endl;
        }
    }
    else
    {
        cerr << "Between 1 and 4 arguments are required, use \"--help\" flag for details" << endl;
    }

    return 0;
}
