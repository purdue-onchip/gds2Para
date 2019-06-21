/**
 * @file   TestMain.cpp 
 * @author Michael R. Hayashi
 * @date   18 October 2018
 * @brief  Primary Function for Input/Solver/Output Control Modes
 */

#define _USE_MATH_DEFINES // Place before including <cmath> for e, log2(e), log10(e), ln(2), ln(10), pi, pi/2, pi/4, 1/pi, 2/pi, 2/sqrt(pi), sqrt(2), and 1/sqrt(2)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <unordered_set>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
#include <parser-spef/parser-spef.hpp>
#include <Eigen/Sparse>
#include "limboint.h"
#include "solnoutclass.h"

// Debug testing macros (comment out if not necessary)
#define SKIP_GENERATE_STIFF

// Manipulate namespace
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
            cout << "Help for LayoutAnalyzer binary (with main testing features)" << endl;
            cout << "Usage: LayoutAnalyzer [options] file1 [file2..file3]" << endl;
            cout << "Options:" << endl;
            cout << "  --help                Display this information." << endl;
            cout << "  --version             Print the version number." << endl;
            cout << "  -r, --read            Read given GDSII file into memory." << endl;
            cout << "  -p, --parrot          Immediately output given GDSII file after reading." << endl;
            cout << "  -w, --write           Write database in memory to given SPEF file." << endl;
            cout << "  -i, --imp             Read given interconnect modeling platform file and write GDSII file with name also given." << endl;
            cout << "  -s, --simulate        Read GDSII file and sim input file into memory, simulate, and write solution to Xyce (SPICE) subcircuit file." << endl;
            cout << "  -sx, --xyce           Identical to \"-s\"." << endl;
            cout << "  -sp, --spef           Read GDSII file and sim input file into memory, simulate, and write solution to SPEF file." << endl;
            cout << endl << "Comments:" << endl;
            cout << " The file passed after -r, --read, -p, or --parrot must be a Calma GDSII stream file." << endl;
            cout << " The file passed after -w or --write must be a blank SPEF file." << endl;
            cout << " The first file passed after -i or --imp must be a 3D description .imp file, and the second must be a blank .gds file." << endl;
            cout << " The first file passed after -s or --simulate (or -sx or --xyce) must be a Calma GDSII stream file, the second must be a sim_input file, and the third must be a blank Xyce file." << endl;
            cout << " The first file passed after -sp or --spef must be a Calma GDSII stream file, the second must be a sim_input file, and the third must be a blank SPEF file." << endl;
            cout << endl << "Bug reporting:" << endl;
            cout << "Visit <https://github.com/purdue-onchip/gds2Para>" << endl;
        }
        else if (strcmp(argv[1], "--version") == 0)
        {
            cout << "Version Number for LayoutAnalyzer binary (beta with main testing features): " << "0.1" << endl;
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
        }
        else if ((strcmp(argv[1], "-p") == 0) || (strcmp(argv[1], "--parrot") == 0))
        {
            // Read and print existing GDSII file
            AsciiDataBase adb;
            string fName = argv[2];
            size_t indExtension = fName.find_last_of(".");
            adb.setFileName(fName.substr(0, indExtension) + "_parrot" + fName.substr(indExtension, string::npos));
            GdsParser::GdsReader adbReader(adb);
            bool adbIsGood = adbReader(fName.c_str());
            adb.print({ });

            // Dump to parroted file immediately
            adb.dump();
            cout << "Dumped parroted file" << endl;

            // Read GDSII file of outline for PSLG purposes
            /*AsciiDataBase adbOutline;
            adbOutline.setFileName(fName.substr(0, indExtension) + "_BACK_DRILL.gds"); // Taking specific file name (to be changed later)
            std::ifstream outlineFile(adbOutline.getFileName().c_str());
            GdsParser::GdsReader adbOutlineReader(adbOutline);
            bool adbOutlineGood = adbOutlineReader(outlineFile);
            vector<complex<double>> outlinePt = adbOutline.findPoints(adbOutline.getCell(0).getCellName(), 0., 0.);*/
            //vector<complex<double>> outlinePt = { complex<double>(+150e-6, -48.0e-6), complex<double>(+150e-6, +121e-6), complex<double>(-1.00e-6, +121e-6), complex<double>(-1.00e-6, -48.0e-6) }; // nand2 outline
            //vector<complex<double>> outlinePt = { complex<double>(+12.77e-6, -0.230e-6), complex<double>(+12.77e-6, +3.03e-6), complex<double>(-0.230e-6, +3.03e-6), complex<double>(-0.230e-6, -0.230e-6) }; // SDFFRS_X2 outline
            //vector<complex<double>> outlinePt = { complex<double>(+2.82e-3, +5.00e-5), complex<double>(+2.825e-3, +3.87e-3), complex<double>(+3.00e-5, +3.87e-3), complex<double>(+3.00e-5, +5.00e-5) }; // 4004 outline

            // Convert to planar straight-line graph (PSLG) file for external meshing
            /*adb.setFileName("examples/4004.gds");
            vector<int> layers = adb.findLayers();
            for (size_t indLayer = 0; indLayer < layers.size(); indLayer++)
            {
                adb.convertPSLG(adb.getCell(adb.getNumCell() - 1).getCellName(), layers[indLayer], outlinePt);
            }
            cout << "Created PSLG file for each layer" << endl;*/
        }
        else if ((strcmp(argv[1], "-w") == 0) || (strcmp(argv[1], "--write") == 0))
        {
            // Load sample information in memory
            string design = "test_out";
            Waveforms blank;

            // Setup the port information vector
            vector<Port> ports = {};
            ports.emplace_back(Port("inp1", 'I', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u1:a", 'I', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("inp2", 'I', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u1:b", 'I', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("out", 'O', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u3:o", 'O', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u1:o", 'O', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u4:a", 'I', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u4:o", 'O', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("f1:d", 'I', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("f1:a", 'O', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u2:a", 'I', 50., vector<double>(6, 0), -1));
            ports.emplace_back(Port("u4:b", 'I', 50., vector<double>(6, 0), -1));

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
            bool couldDump = sdb.printDumpSPEF();
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
        if ((strcmp(argv[1], "-s") == 0) || (strcmp(argv[1], "--simulate") == 0) || (strcmp(argv[1], "-sx") == 0) || (strcmp(argv[1], "--xyce") == 0) || (strcmp(argv[1], "-sp") == 0) || (strcmp(argv[1], "--spef") == 0))
        {
            // Initialize SolverDataBase, mesh, and set variables for performance tracking
            clock_t t1 = clock();
            SolverDataBase sdb;
            fdtdMesh sys;
            int status = 0; // Initialize as able to return successfully
            bool adbIsGood, sdbIsGood, sdbCouldDump;

            // Get file names
            string inGDSIIFile = argv[2];
            string inSimFile = argv[3];
            size_t indExtension = inGDSIIFile.find_last_of(".");

            // Read GDSII file
            AsciiDataBase adb;
            adb.setFileName(inGDSIIFile);
            GdsParser::GdsReader adbReader(adb);
            adbIsGood = adbReader(inGDSIIFile.c_str());
            if (adbIsGood)
            {
                string topCellName = adb.getCell(adb.getNumCell() - 1).getCellName();
                vector<size_t> indCellPrint = {}; // { adb.getNumCell() - 1 };
                adb.saveToMesh(topCellName, { 0., 0. }, strans(), &sys); // Recursively save GDSII conductor information to sys
                adb.print(indCellPrint);
                cout << "GDSII file read" << endl;
            }
            else
            {
                cerr << "Unable to read in GDSII file" << endl;
                status = 1;
                return status;
            }

            // Read simulation input file
            sdbIsGood = sdb.readSimInput(inSimFile);
            if (sdbIsGood)
            {
                cout << "Simulation input file read" << endl;
            }
            else
            {
                cerr << "Unable to read in simulation input file" << endl;
                status = 1;
                return status;
            }

            // Append information so far to fdtdMesh
            unordered_set<double> portCoorx, portCoory;
            sdb.convertToFDTDMesh(&sys, adb.getNumCdtIn(), &portCoorx, &portCoory);

            // Mesh the domain and mark conductors
            unordered_map<double, int> xi, yi, zi;
            clock_t t2 = clock();
            status = meshAndMark(&sys, xi, yi, zi, &portCoorx, &portCoory);
            if (status == 0)
            {
                cout << "meshAndMark Success!" << endl;
                cout << "meshAndMark time is " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
            }
            else
            {
                cerr << "meshAndMark Fail!" << endl;
                return status;
            }
            //sys.print();

            // Set D_eps and D_sig
            clock_t t3 = clock();
            status = matrixConstruction(&sys);
            if (status == 0)
            {
                cout << "matrixConstruction Success!" << endl;
                cout << "matrixConstruction time is " << (clock() - t3) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
            }
            else {
                cerr << "matrixConstruction Fail!" << endl;
                return status;
            }
            //sys.print();

            // Set port
            clock_t t4 = clock();
            status = portSet(&sys, xi, yi, zi);
            if (status == 0)
            {
                cout << "portSet Success!" << endl;
                cout << "portSet time is " << (clock() - t4) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
            }
            else
            {
                cerr << "portSet Fail!" << endl;
                return status;
            }
            //sys.print();

            // Generate Stiffness Matrix
#ifndef SKIP_GENERATE_STIFF
            clock_t t5 = clock();
            status = generateStiff(&sys);
            if (status == 0)
            {
                cout << "generateStiff Success!" << endl;
                cout << "generateStiff time is " << (clock() - t5) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
            }
            else
            {
                cerr << "generateStiff Fail!" << endl;
                return status;
            }
#endif
            cout << "Start to generate parameter!\n";
            // Parameter generation
            clock_t t6 = clock();
            status = paraGenerator(&sys, xi, yi, zi);
            if (status == 0)
            {
                cout << "paraGenerator Success!" << endl;
                cout << "paraGenerator time is " << (clock() - t6) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
            }
            else
            {
                cerr << "paraGenerator Fail!" << endl;
                return status;
            }
            cout << "Engine time to this point: " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
            cout << "Total time to this point: " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;

            // Parameter storage
            spMat matG(sys.numPorts, sys.numPorts); // Initialize Eigen sparse conductance matrix (S)
            spMat matC(sys.numPorts, sys.numPorts); // Initialize Eigen sparse capacitance matrix (F)
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
                    // Find Y-parameters from Z-parameters (temporary measure for diagonal elements only)]
                    double gij = sys.x[indi * sys.numPorts + indj].real() / norm(sys.x[indi * sys.numPorts + indj]); // Re(Z^-1) = Re(Z) / |Z|^2
                    double gji = sys.x[indj * sys.numPorts + indi].real() / norm(sys.x[indj * sys.numPorts + indi]);
                    double bij = -sys.x[indi * sys.numPorts + indj].imag() / norm(sys.x[indi * sys.numPorts + indj]); // Im(Z^-1) = -Im(Z) / |Z|^2
                    double bji = -sys.x[indj * sys.numPorts + indi].imag() / norm(sys.x[indj * sys.numPorts + indi]);

                    // Symmetrize diagonal components of admittance matrix before saving entry
                    double symgij = 0.5 * (gij + gji);
                    double symcij = 0.5 * (bij + bji) / (2. * M_PI * sys.freqStart * sys.freqUnit);
                    //double symgij = 0.5 * (sys.Y[indi * sys.numPorts + indj].real() + sys.Y[indj * sys.numPorts + indi].real());
                    //double symcij = 0.5 * (sys.Y[indi * sys.numPorts + indj].imag() + sys.Y[indj * sys.numPorts + indi].imag()) / (2. * M_PI * sys.freqStart * sys.freqUnit);
                    listG.push_back(dTriplet(indj, indi, symgij));
                    listC.push_back(dTriplet(indj, indi, symcij));
                    /*if (indi != indj) // Off-diagonal entries are negated node-to-node admittances
                    {
                        listG.push_back(dTriplet(indi, indj, sys.Y[indi * sys.numPorts + indj].real()));
                        listC.push_back(dTriplet(indi, indj, sys.Y[indi * sys.numPorts + indj].imag() / (2. * M_PI* sys.freqStart * sys.freqUnit)));
                    }
                    
                    // Diagonal entries of bus matrix need to have sum of all attached admittances, so note them
                    sumG += sys.Y[indi * sys.numPorts + indj].real();
                    sumC += sys.Y[indi * sys.numPorts + indj].imag() / (2 * M_PI* sys.freqStart * sys.freqUnit);*/
                }

                // Record diagonal entries
                listG.push_back(dTriplet(indi, indi, sumG));
                listC.push_back(dTriplet(indi, indi, sumC));
            }
            matG.setFromTriplets(listG.begin(), listG.end()); // Assign nonzero entries to sparse conductance matrix
            matC.setFromTriplets(listC.begin(), listC.end()); // Do not put in compressed sparse row (CSR) format due to density
            Parasitics oldPara = sdb.getParasitics(); // Get outdated parastics structure to update
            sdb.setParasitics(Parasitics(oldPara.getPorts(), matG, matC));

            // Select Output File Based on Control Mode
            if ((strcmp(argv[1], "-sp") == 0) || (strcmp(argv[1], "--spef") == 0))
            {
                // Output SPEF file
                string outSPEFFile = argv[4];
                sdb.setDesignName(adb.findNames().back());
                sdb.setOutSPEF(outSPEFFile);
                bool sdbCouldDump = sdb.printDumpSPEF();
                cout << "File ready at " << outSPEFFile << endl;
            }
            else
            {
                // Output Xyce subcircuit file
                string outXyceFile = argv[4];
                vector<size_t> indLayerPrint = {0, sdb.getNumLayer() / 2, sdb.getNumLayer() - 1}; // {}; // Can use integer division
                sdb.setDesignName(adb.findNames().back());
                sdb.setOutXyce(outXyceFile);
                sdbCouldDump = sdb.printDumpXyce(indLayerPrint);
                cout << "File ready at " << outXyceFile << endl;
            }
        }
        else
        {
            cerr << "Must pass a GDSII file, sim_input file, and blank Xyce file to write after \"-s\" or \"-sx\" flag, rerun with \"--help\" flag for details" << endl;
            cerr << "Must pass a GDSII file, sim_input file, and blank SPEF file to write after \"-sp\" flag, rerun with \"--help\" flag for details" << endl;
        }
    }
    else
    {
        cerr << "Between 1 and 4 arguments are required, use \"--help\" flag for details" << endl;
    }

    return 0;
}
