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
#include "limboint.hpp"
#include "solnoutclass.hpp"
#include "layeredFdtd.hpp"

// Debug testing macros (comment out if not necessary)
//#define SKIP_GENERATE_STIFF
//#define SKIP_WRITE_SYS_TO_FILE

// Manipulate namespace
using std::cerr;
using std::cout;
using std::endl;

int readTXT(fdtdMesh* sys, string file);
double fdtdGetValue(char *str);
int fdtdStringWord(char *s, char *word[]);

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
            cout << "Usage: mpirun LayoutAnalyzer options file1 [file2..file3]" << endl;
            cout << "Options:" << endl;
            cout << "  --help                Display this information." << endl;
            cout << "  --version             Print the version number." << endl;
            cout << "  -r, --read            Read given GDSII file into memory and display statistics." << endl;
            cout << "  -p, --parrot          Immediately output given GDSII file after reading." << endl;
            cout << "  -w, --write           Write database in memory to given SPEF file." << endl;
            cout << "  -i, --imp             Read given interconnect modeling platform file and write GDSII file with name also given." << endl;
            cout << "  -t, --pslg            Read GDSII files of design and outline and write a PSLG file for each layer." << endl;
            cout << "  -s, --simulate        Read GDSII and sim input files, simulate, and write solution to Xyce (SPICE) subcircuit." << endl;
            cout << "  -sx, --xyce           Identical to \"-s\"." << endl;
            cout << "  -sp, --spef           Read GDSII and sim input files into memory, simulate, and write solution to SPEF file." << endl;
            cout << "  -sc, --citi           Read GDSII and sim input files into memory, simulate, and write solution to CITIfile." << endl;
            cout << "  -st, --touchstone     Read GDSII and sim input files into memory, simulate, and write solution to Touchstone file." << endl;
            cout << endl << "Comments:" << endl;
            cout << " The file passed after -r, --read, -p, or --parrot must be a Calma GDSII stream file." << endl;
            cout << " The file passed after -w or --write must be a blank SPEF file." << endl;
            cout << " The first file passed after -i or --imp must be a 3D description .imp file, and the second must be a blank .gds file." << endl;
            cout << " The first file passed after -t or --pslg must be a Calma GDSII stream file of the design, and the second must be a GDSII file of just the design's outline." << endl;
            cout << " The first file passed after -s or --simulate (or -sx or --xyce) must be a Calma GDSII stream file, the second must be a sim_input file, and the third must be a blank Xyce file." << endl;
            cout << " The first file passed after -sp or --spef must be a Calma GDSII stream file, the second must be a sim_input file, and the third must be a blank SPEF file." << endl;
            cout << " The first file passed after -sc or --citi must be a Calma GDSII stream file, the second must be a sim_input file, and the third must be a blank CITIfile." << endl;
            cout << " The first file passed after -sc or --citi must be a Calma GDSII stream file, the second must be a sim_input file, and the third must be a blank Touchstone file." << endl;
            cout << endl << "Bug reporting:" << endl;
            cout << "Visit <https://github.com/purdue-onchip/gds2Para>" << endl;
        }
        else if (strcmp(argv[1], "--version") == 0)
        {
            cout << "Version Number for LayoutAnalyzer binary (beta with main testing features): " << "0.2" << endl;
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
            GdsParser::GdsReader adbReader(adb);
            bool adbIsGood = adbReader(fName.c_str());
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
            adb.print({});

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
            ports.emplace_back(Port("inp1", 'I', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u1:a", 'I', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("inp2", 'I', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u1:b", 'I', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("out", 'O', 50.0, 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u3:o", 'O', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u1:o", 'O', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u4:a", 'I', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u4:o", 'O', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("f1:d", 'I', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("f1:a", 'O', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u2:a", 'I', 50., 1, vector<double>(6, 0), -1));
            ports.emplace_back(Port("u4:b", 'I', 50., 1, vector<double>(6, 0), -1));

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
            matC.makeCompressed(); // Capacitance matrix in compressed sparse row (CSR) format

            // Create variables of custom classes
            Parasitics sample(ports, matG, matC, { 1000. });
            SolverDataBase sdb(design, blank, sample);

            // Prepare to write to file
            string fName = argv[2];
            sdb.setOutSPEF(fName);
            sdb.print({});
            bool couldDump = sdb.dumpSPEF();
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
        else if ((strcmp(argv[1], "-t") == 0) || (strcmp(argv[1], "--pslg") == 0))
        {
            // Read and print existing GDSII file
            AsciiDataBase adbDesign;
            string designFileName = argv[2];
            adbDesign.setFileName(designFileName);
            GdsParser::GdsReader adbReader(adbDesign);
            bool adbDesignGood = adbReader(designFileName.c_str());
            adbDesign.print({});

            // Read GDSII file of outline for PSLG purposes
            AsciiDataBase adbOutline;
            string outlineFileName = argv[2];
            adbOutline.setFileName(outlineFileName);
            GdsParser::GdsReader adbOutlineReader(adbOutline);
            bool adbOutlineGood = adbOutlineReader(outlineFileName.c_str());
            vector<complex<double>> outlinePt = adbOutline.findPoints(adbOutline.getCell(0).getCellName(), { 0., 0. }, strans());
            //vector<complex<double>> outlinePt = { complex<double>(+150e-6, -48.0e-6), complex<double>(+150e-6, +121e-6), complex<double>(-1.00e-6, +121e-6), complex<double>(-1.00e-6, -48.0e-6) }; // nand2 outline
            //vector<complex<double>> outlinePt = { complex<double>(+12.77e-6, -0.230e-6), complex<double>(+12.77e-6, +3.03e-6), complex<double>(-0.230e-6, +3.03e-6), complex<double>(-0.230e-6, -0.230e-6) }; // SDFFRS_X2 outline
            //vector<complex<double>> outlinePt = { complex<double>(+2.82e-3, +5.00e-5), complex<double>(+2.825e-3, +3.87e-3), complex<double>(+3.00e-5, +3.87e-3), complex<double>(+3.00e-5, +5.00e-5) }; // 4004 outline

            // Convert to planar straight-line graph (PSLG) file for external meshing
            vector<int> layers = adbDesign.findLayers();
            for (size_t indLayer = 0; indLayer < layers.size(); indLayer++)
            {
                adbDesign.convertPSLG(adbDesign.getCell(adbDesign.getNumCell() - 1).getCellName(), layers[indLayer], outlinePt);
            }
            cout << "Created PSLG file for each layer" << endl;
        }
		else if ((strcmp(argv[1], "-rt") == 0) || (strcmp(argv[1], "--readtxt") == 0))
		{   // read the txt file for all the inputs and output citi file
			// Initialize SolverDataBase, mesh, and set variables for performance tracking
			clock_t t1 = clock();
			fdtdMesh sys;
			int status = 0; // Initialize as able to return successfully
			
			readTXT(&sys, argv[2]);

			// Mesh the domain and mark conductors
			unordered_map<double, int> xi, yi, zi;
			unordered_set<double> portCoorx, portCoory;

			clock_t t2 = clock();
			status = meshAndMark(&sys, xi, yi, zi, &portCoorx, &portCoory);
			if (status == 0)
			{
				cout << "meshAndMark Success!" << endl;
				cout << "meshAndMark time is " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
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
				cout << "matrixConstruction time is " << (clock() - t3) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
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
				cout << "portSet time is " << (clock() - t4) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
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
				cout << "generateStiff time is " << (clock() - t5) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
			}
			else
			{
				cerr << "generateStiff Fail!" << endl;
				return status;
			}
#endif

			// Write object sys to files
#ifndef SKIP_WRITE_SYS_TO_FILE
			WriteSysToFile(sys);
#endif

			// Parameter generation
			clock_t t6 = clock();
			status = paraGenerator(&sys, xi, yi, zi);
			if (status == 0)
			{
				cout << "paraGenerator Success!" << endl;
				cout << "paraGenerator time is " << (clock() - t6) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
			}
			else
			{
				cerr << "paraGenerator Fail!" << endl;
				return status;
			}
			cout << "Engine time to this point: " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
			cout << "Total time to this point: " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;

		}
        else
        {
            cerr << "Must pass a .imp file to read and blank GDSII file to write after \"-i\" flag, rerun with \"--help\" flag for details" << endl;
            cerr << "Must pass two related GDSII files to read after \"-t\" flag, rerun with \"--help\" flag for details" << endl;
        }
    }
    else if (argc == 5)
    {
        if ((strcmp(argv[1], "-s") == 0) || (strcmp(argv[1], "--simulate") == 0) || (strcmp(argv[1], "-sx") == 0) || (strcmp(argv[1], "--xyce") == 0) || (strcmp(argv[1], "-sp") == 0) || (strcmp(argv[1], "--spef") == 0) || (strcmp(argv[1], "-sc") == 0) || (strcmp(argv[1], "--citi") == 0) || (strcmp(argv[1], "-st") == 0) || (strcmp(argv[1], "--touchstone") == 0))
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
                vector<size_t> indCellPrint = {}; // { adb.getNumCell() - 1 };
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
            string topCellName = adb.getCell(adb.getNumCell() - 1).getCellName();
            adb.saveToMesh(topCellName, { 0., 0. }, strans(), &sys, sdb.findLayerIgnore()); // Recursively save GDSII conductor information to sys
            sdb.convertToFDTDMesh(&sys, adb.getNumCdtIn(), &portCoorx, &portCoory); // Save simulation input information to sys

            // Mesh the domain and mark conductors
            unordered_map<double, int> xi, yi, zi;
            clock_t t2 = clock();
            status = meshAndMark(&sys, xi, yi, zi, &portCoorx, &portCoory);
            if (status == 0)
            {
                cout << "meshAndMark Success!" << endl;
                cout << "meshAndMark time is " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
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
                cout << "matrixConstruction time is " << (clock() - t3) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
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
                cout << "portSet time is " << (clock() - t4) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
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
                cout << "generateStiff time is " << (clock() - t5) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
            }
            else
            {
                cerr << "generateStiff Fail!" << endl;
                return status;
            }
#endif

            // Write object sys to files
#ifndef SKIP_WRITE_SYS_TO_FILE
            WriteSysToFile(sys);
#endif

            // Parameter generation
            clock_t t6 = clock();
            status = paraGenerator(&sys, xi, yi, zi);
            if (status == 0)
            {
                cout << "paraGenerator Success!" << endl;
                cout << "paraGenerator time is " << (clock() - t6) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
            }
            else
            {
                cerr << "paraGenerator Fail!" << endl;
                return status;
            }
            cout << "Engine time to this point: " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
            cout << "Total time to this point: " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;

            // Network parameter storage
            Parasitics newPara = sdb.getParasitics(); // Start with outdated parastics to update
            newPara.saveNetworkParam('Z', sdb.getSimSettings().getFreqsHertz(), sys.x); // Save the Z-parameters in fdtdMesh to Parasitics class
            sdb.setParasitics(newPara);

            // Select Output File Based on Control Mode
            cout << endl;
            if ((strcmp(argv[1], "-sp") == 0) || (strcmp(argv[1], "--spef") == 0))
            {
                // Output SPEF file
                string outSPEFFile = argv[4];
                vector<size_t> indLayerPrint = { 0, 1 * sdb.getNumLayer() / 3, 2 * sdb.getNumLayer() / 3, sdb.getNumLayer() - 1 }; // {}; // Can use integer division
                sdb.setDesignName(adb.findNames().back());
                sdb.setOutSPEF(outSPEFFile);
                sdb.print(indLayerPrint);
                bool sdbCouldDump = sdb.dumpSPEF();
                cout << "File ready at " << outSPEFFile << endl;
            }
            else if ((strcmp(argv[1], "-sc") == 0) || (strcmp(argv[1], "--citi") == 0))
            {
                // Convert to S-parameters
                Parasitics newPara = sdb.getParasitics(); // Start with copy of parastics to reinterpret
                newPara.convertParam('S');
                sdb.setParasitics(newPara);

                // Output Common Instrumentation Transfer and Interchange file (CITIfile)
                string outCITIFile = argv[4];
                vector<size_t> indLayerPrint = { 0, 1 * sdb.getNumLayer() / 3, 2 * sdb.getNumLayer() / 3, sdb.getNumLayer() - 1 }; // {}; // Can use integer division
                sdb.setDesignName(adb.findNames().back());
                sdb.setOutCITI(outCITIFile);
                sdb.print(indLayerPrint);
                bool sdbCouldDump = sdb.dumpCITI();
                cout << "File ready at " << outCITIFile << endl;
            }
            else if ((strcmp(argv[1], "-st") == 0) || (strcmp(argv[1], "--touchstone") == 0))
            {
                // Output Touchstone file
                string outTstoneFile = argv[4];
                vector<size_t> indLayerPrint = { 0, 1 * sdb.getNumLayer() / 3, 2 * sdb.getNumLayer() / 3, sdb.getNumLayer() - 1 }; // {}; // Can use integer division
                sdb.setDesignName(adb.findNames().back());
                sdb.setOutTouchstone(outTstoneFile);
                sdb.print(indLayerPrint);
                bool sdbCouldDump = sdb.dumpTouchstone();
                cout << "File ready at " << outTstoneFile << endl;
            }
            else
            {
                // Output Xyce subcircuit file
                string outXyceFile = argv[4];
                vector<size_t> indLayerPrint = { 0, sdb.getNumLayer() / 2, sdb.getNumLayer() - 1 }; // {}; // Can use integer division
                sdb.setDesignName(adb.findNames().back());
                sdb.setOutXyce(outXyceFile);
                sdb.print(indLayerPrint);
                sdbCouldDump = sdb.dumpXyce();
                cout << "File ready at " << outXyceFile << endl;
            }
        }
        else
        {
            cerr << "Must pass a GDSII file, sim_input file, and blank Xyce file to write after \"-s\" or \"-sx\" flag" << endl;
            cerr << "Must pass a GDSII file, sim_input file, and blank SPEF file to write after \"-sp\" flag" << endl;
            cerr << "Must pass a GDSII file, sim_input file, and blank CITI file to write after \"-sc\" flag" << endl;
            cerr << "Must pass a GDSII file, sim_input file, and blank Touchstone file to write after \"-st\" flag" << endl;
            cerr << "Rerun with \"--help\" flag for details" << endl;
        }
    }
    else
    {
        cerr << "Between 1 and 4 arguments are required, use \"--help\" flag for details" << endl;
    }

    return 0;
}

int readTXT(fdtdMesh* sys, string file) {
	/* read the txt file */
	FILE *fp, *cfp;
	char fbase[MAXC];
	char s[MAXC];
	char *word[MAXC];
	
	fp = fopen(&file[0], "r");
	if (fp == NULL) {
		perror("Failed: ");
		return -2;
	}
	
	cout << "Begin reading the txt file!" << endl;
	
	// read the conductorIn information
	fgets(s, MAXC, fp);
	fdtdStringWord(s, word);
	sys->lengthUnit = 1;
	sys->xlim1 = fdtdGetValue(word[0]);
	sys->xlim2 = fdtdGetValue(word[1]);
	sys->ylim1 = fdtdGetValue(word[2]);
	sys->ylim2 = fdtdGetValue(word[3]);
	sys->zlim1 = fdtdGetValue(word[4]);
	sys->zlim2 = fdtdGetValue(word[5]);
	while (strncmp(word[0], "FREQUENCY:", 10)) {
		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
	}
	if (!strncmp(word[0], "FREQUENCY:", 10)) {
		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
		sys->freqUnit = fdtdGetValue(word[2]);

		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
		sys->freqStart = fdtdGetValue(word[2]);

		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
		sys->freqEnd = fdtdGetValue(word[2]);

		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
		sys->nfreq = fdtdGetValue(word[2]);

		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
		sys->freqScale = fdtdGetValue(word[2]);
	}
	while (strncmp(word[0], "CONDUCTORIN:", 12)) {
		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
	}
	if (!strncmp(word[0], "CONDUCTORIN:", 12)) {
		while (fgets(s, MAXC, fp) != NULL) {
			fdtdStringWord(s, word);
			if (!strncmp(word[0], "STACK", 5)) {
				break;
			}
			fdtdOneCondct base;
			sys->conductorIn.push_back(base);
			int i = sys->conductorIn.size() - 1;
			sys->conductorIn[i].numVert = (int)fdtdGetValue(word[0]);
			sys->conductorIn[i].zmin = fdtdGetValue(word[1]);
			sys->conductorIn[i].zmax = fdtdGetValue(word[2]);
			sys->conductorIn[i].x = (double*)calloc(sys->conductorIn[i].numVert, sizeof(double));
			sys->conductorIn[i].y = (double*)calloc(sys->conductorIn[i].numVert, sizeof(double));
			for (int j = 0; j < sys->conductorIn[i].numVert; ++j) {
				sys->conductorIn[i].x[j] = fdtdGetValue(word[j * 2 + 3]);
				sys->conductorIn[i].y[j] = fdtdGetValue(word[j * 2 + 4]);
			}
		}
	}
	sys->numCdtRow = sys->conductorIn.size();

	// read the stack information
	while (strncmp(word[0], "STACK", 5)) {
		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
	}
	if (fgets(s, MAXC, fp) != NULL) {
		fdtdStringWord(s, word);
		sys->numStack = (int)fdtdGetValue(word[0]);
		//cout << "NumStack is " << word[0] << endl;
		for (int i = 0; i < sys->numStack; ++i) {
			fgets(s, MAXC, fp);
			fdtdStringWord(s, word);
			sys->stackBegCoor.push_back(fdtdGetValue(word[0]));
			sys->stackEndCoor.push_back(fdtdGetValue(word[1]));
			sys->stackEps.push_back(fdtdGetValue(word[2]));
			sys->stackSig.push_back(fdtdGetValue(word[3]));
		}
	}
	

	// read the port information
	while (strncmp(word[0], "PORT", 4)) {
		fgets(s, MAXC, fp);
		fdtdStringWord(s, word);
	}
	if (fgets(s, MAXC, fp) != NULL) {
		fdtdStringWord(s, word);
		sys->numPorts = (int)fdtdGetValue(word[0]);
		cout << "Numports is " << word[0] << endl;
		for (int i = 0; i < sys->numPorts; ++i) {
			fgets(s, MAXC, fp);
			fdtdStringWord(s, word);
			fdtdPort base;
			sys->portCoor.push_back(base);
			sys->portCoor[i].multiplicity = (int)fdtdGetValue(word[0]);
			cout << "port " << i << " multiplicity is " << word[0] << endl;
			for (int j = 0; j < sys->portCoor[i].multiplicity; ++j) {
				sys->portCoor[i].x1.push_back(fdtdGetValue(word[j * 6 + 1]));
				sys->portCoor[i].y1.push_back(fdtdGetValue(word[j * 6 + 2]));
				sys->portCoor[i].z1.push_back(fdtdGetValue(word[j * 6 + 3]));
				sys->portCoor[i].x2.push_back(fdtdGetValue(word[j * 6 + 4]));
				sys->portCoor[i].y2.push_back(fdtdGetValue(word[j * 6 + 5]));
				sys->portCoor[i].z2.push_back(fdtdGetValue(word[j * 6 + 6]));
				sys->portCoor[i].portDirection.push_back(fdtdGetValue(word[j * 6 + 7]));
			}
			
		}
	}
	return 0;
}

double fdtdGetValue(char *str)
{
	double value;                            /* converted true value      */
	sscanf(str, "%lf", &value);
	return value;
}

int fdtdStringWord(char *s, char *word[]) {
	int wno;                      /* word char counter */
	int ctr;                       /* char flag         */

	ctr = 1;
	wno = 0;
	for (; *s != '\0'; s++)
		if (*s == ' ' || *s == '\t' || *s == ',' || *s == '(' ||
			*s == ')' || *s == '\n') {
			*s = '\0';
			ctr = 1;
		}
		else {
			*s = (char)toupper(*s);
			if (ctr == 1) {
				word[wno++] = s;
				ctr = 0;
			}
		}
		return(wno);
}