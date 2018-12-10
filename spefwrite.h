/**
* File: spefwrite.h
* Author: Michael R. Hayashi
* Date: 6 November 2018
* Desc: Header for Writing SPEF Files from Memory
*/

// Make sure this file is included only once
#ifndef spefwrite_h
#define spefwrite_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <ctime>
//#include <algorithm>
#include <parser-spef/parser-spef.hpp>
using namespace std;

class Waveforms
{
private:
    std::string name;             // Placeholder name
public:
    // Default constructor
    Waveforms()
    {
        this->name = "";
    }

    // Get name
    std::string getName() const
    {
        return this->name;
    }

    // Print the waveforms information
    void print() const
    {
        cout << " ------" << endl;
        cout << " No Waveform Information" << endl;
        cout << " ------" << endl;
    }

    // Destructor
    ~Waveforms()
    {
        this->name = "";
    }
};

class Parasitics
{
private:
    size_t nPorts;                // Number of ports
    vector<std::string> ports;    // Name of each port
    vector<char> portDir;         // Direction of each port
    vector<vector<double>> matG;  // Conductance matrix (S)
    vector<vector<double>> matC;  // Capacitance matrix (F)
public:
    // Default constructor
    Parasitics()
    {
        vector<std::string> ports = {};
        vector<char> portDir = {};
        vector<vector<double>> emptMat = { {} };
        this->nPorts = 0;
        this->ports = ports;
        this->portDir = portDir;
        this->matG = emptMat;
        this->matC = emptMat;
    }

    // Parametrized constructor
    Parasitics(size_t nPorts, vector<std::string> ports, vector<char> portDir, vector<vector<double>> matG, vector<vector<double>> matC)
    {
        this->nPorts = nPorts;
        this->ports = ports;
        this->portDir = portDir;
        this->matG = matG;
        this->matC = matC;
    }

    // Get number of ports
    size_t getNPort() const
    {
        return this->nPorts;
    }

    // Get port names
    vector<std::string> getPorts() const
    {
        return this->ports;
    }

    // Get port directions
    // ('I' = input, 'O' = output, 'B' = bidirectional)
    vector<char> getPortDir() const
    {
        return this->portDir;
    }

    // Get conductance matrix
    vector<vector<double>> getGMatrix() const
    {
        return this->matG;
    }

    // Get capacitance matrix
    vector<vector<double>> getCMatrix() const
    {
        return this->matC;
    }

    // Return node-to-ground conductance
    double getGNodeGround(size_t indNode) const
    {
        double condNG = 0.0;
        for (size_t indi = 0; indi < this->getNPort(); indi++)
        {
            // Calculate by summing along given row
            condNG += ((this->matG)[indi])[0];
        }
        return condNG;
    }

    // Return total conductance represented in matrix
    double getGTotal() const
    {
        size_t nPort = this->getNPort();
        double condTot = 0.0;
        for (size_t indi = 0; indi < nPort; indi++)
        {
            condTot += getGNodeGround(indi); // Include node-to-ground conductances
            for (size_t indj = indi + 1; indj < nPort; indj++)
            {
                condTot -= (this->matG)[indi][indj]; // Include negated node-node conductance entries
            }
        }
        return condTot;
    }

    // Return node-to-ground capacitance
    double getCNodeGround(size_t indNode) const
    {
        double capNG = 0.0;
        for (size_t indi = 0; indi < this->getNPort(); indi++)
        {
            // Calculate by summing along given row
            capNG += ((this->matC)[indi])[0];
        }
        return capNG;
    }

    // Return total capacitance represented in matrix
    double getCTotal() const
    {
        size_t nPort = this->getNPort();
        double capTot = 0.0;
        for (size_t indi = 0; indi < nPort; indi++)
        {
            capTot += getCNodeGround(indi); // Include node-to-ground capacitances
            for (size_t indj = indi + 1; indj < nPort; indj++)
            {
                capTot -= (this->matC)[indi][indj]; // Include negated node-node capacitance entries
            }
        }
        return capTot;
    }

    // Print the parasitics information
    void print() const
    {
        int numPort = getNPort();
        cout << " ------" << endl;
        cout << " Parasitics Details:" << endl;
        cout << "  Port Names:" << endl;
        for (size_t indi = 0; indi < numPort; indi++) // Handle each port name
        {
            cout << "   #" << indi + 1 << ": " << (this->ports)[indi] << endl;
        }
        cout << "  Conductance Matrix (S):" << endl;
        for (size_t indi = 0; indi < numPort; indi++)
        {
            cout << "   |"; // Indent the matrix rows
            for (size_t indj = 0; indj < numPort; indj++)
            {
                char value[12];
                sprintf(value, "%-+10.2g", (this->matG)[indi][indj]);
                cout << value;
                //delete[] value; // Free array pointer
            }
            cout << "|" << endl; // End the matrix row
        }
        cout << "  Capacitance Matrix (F):" << endl;
        for (size_t indi = 0; indi < numPort; indi++)
        {
            cout << "   |"; // Indent the matrix rows
            for (size_t indj = 0; indj < numPort; indj++)
            {
                char value[12];
                sprintf(value, "%-+10.2g", (this->matC)[indi][indj]);
                cout << value;
                //delete[] value; // Free array pointer
            }
            cout << "|" << endl; // End the matrix row
        }
        cout << " ------" << endl;
    }

    // Translate parasitics to Spef struct
    spef::Spef toSPEF(std::string designName, double saveThresh)
    {
        // Initialize variables
        spef::Spef para; // Spef struct for parasitics
        char timeStr[80]; // Character array to hold formatted time
        time_t rawtime; // Timer variable
        strftime(timeStr, sizeof(timeStr), "%d-%m-%Y %H:%M:%S", localtime(&rawtime)); // Use local timezone to format string
        std::string time(timeStr); // Formatted time to string parametrized constructor
        size_t numPort = this->getNPort(); // Number of ports

        // Populate Spef struct header fields
        para.standard = "\"IEEE 1481-1998\"";
        para.design_name = "\"" + designName + "\"";
        para.date = "\"" + time + "\"";
        para.vendor = "\"DARPA ERI Contributors\"";
        para.program = "\"SPEF Writer from DARPA ERI\"";
        para.version = "\"1.0\"";
        para.design_flow = "\"NETLIST_TYPE_VERILOG\"";
        para.divider = "/";
        para.delimiter = ":";
        para.bus_delimiter = "[ ]";
        para.time_unit = "1 S";
        para.capacitance_unit = "1 F";
        para.resistance_unit = "1 OHM";
        para.inductance_unit = "1 H";

        // Populate Spef struct name map fields and ports vector
        for (size_t indi = 0; indi < numPort; indi++)
        {
            para.name_map.emplace(indi + 1, (this->ports)[indi]); // Create name map for each port
            //para.ports.emplace_back((this->ports)[indi]); // Instantiate and push new port entry by name
            para.ports.emplace_back("*" + to_string(indi + 1)); // Instantiate and push new port entry by name map
            switch ((this->portDir)[indi]) // Assign port direction
            {
            case 'O':
                para.ports.back().direction = spef::ConnectionDirection::OUTPUT;
            case 'I':
                para.ports.back().direction = spef::ConnectionDirection::INPUT;
            case 'B':
                para.ports.back().direction = spef::ConnectionDirection::INOUT;
            default:
                para.ports.back().direction = spef::ConnectionDirection::INOUT; // Treat as bidirectional if direction unclear
            }
        }

        // Populate Spef struct nets vector (single net)
        double capTot = this->getCTotal();
        double condTot = this->getGTotal();
        para.nets.emplace_back(spef::Net());
        para.nets.back().name = "all";
        para.nets.back().lcap = capTot;
        for (size_t indi = 0; indi < numPort; indi++)
        {
            para.nets.back().connections.emplace_back(spef::Connection());
            para.nets.back().connections.back().name = (this->ports)[indi];
            para.nets.back().connections.back().type = spef::ConnectionType::EXTERNAL; // All ports are external connections, not internal to cells
            para.nets.back().connections.back().direction = (para.ports[indi]).direction; // Same as port direction
            for (size_t indj = 0; indj < numPort; indj++)
            {
                if (indi == indj) // Diagonal element, so admittance to ground
                {
                    if (abs(this->getCNodeGround(indi) / capTot) >= saveThresh) // Only place sufficiently nonzero entries
                    {
                        para.nets.back().caps.emplace_back(forward_as_tuple((para.ports)[indi].name, "", this->getCNodeGround(indi)));
                    }
                    if (abs(this->getGNodeGround(indi) / condTot) >= saveThresh)
                    {
                        para.nets.back().ress.emplace_back(forward_as_tuple((para.ports)[indi].name, "", 1.0 / this->getGNodeGround(indi)));
                    }
                }
                else if (indj > indi) // Upper triangular element, so negate node-to-node admittance
                {
                    if (abs((this->matC)[indi][indj] / capTot) >= saveThresh) // Only place sufficiently nonzero entries
                    {
                        para.nets.back().caps.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[indj].name, -(this->matC)[indi][indj]));
                    }
                    if (abs((this->matG)[indi][indj] / condTot) >= saveThresh)
                    para.nets.back().ress.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[indj].name, -1.0 / ((this->matG)[indi][indj])));
                } // Skip lower triangular element, which is repeated
            }
        }

        // Return the Spef struct
        return para;
    }

    // Destructor
    ~Parasitics()
    {
        vector<std::string> ports = {};
        vector<char> portDir = {};
        vector<vector<double>> emptMat = { {} };
        this->nPorts = 0;
        this->ports = ports;
        this->portDir = portDir;
        this->matG = emptMat;
        this->matC = emptMat;
    }
};

struct SolverDataBase
{
private:
    std::string designName; // Name of design
    Waveforms wf;           // Waveforms
    Parasitics para;        // Parasitics
    std::string outSPEF;    // SPEF output file name
public:
    // Default constructor
    SolverDataBase()
    {
        Waveforms wf;
        Parasitics para;
        this->designName;
        this->wf = wf;
        this->para = para;
        this->outSPEF = "";
    }

    // Parametrized constructor
    SolverDataBase(std::string designName, Waveforms wf, Parasitics para)
    {
        this->designName = designName;
        this->wf = wf;
        this->para = para;
        this->outSPEF = "";
    }

    // Get name of design
    std::string getDesignName() const
    {
        return this->designName;
    }

    // Get waveforms
    Waveforms getWaveforms() const
    {
        return this->wf;
    }

    // Get parasitics
    Parasitics getParasitics() const
    {
        return this->para;
    }

    // Get SPEF output file name
    std::string getOutSPEF() const
    {
        return this->outSPEF;
    }

    // Set SPEF output file name
    void setOutSPEF(std::string fileName)
    {
        this->outSPEF = fileName;
    }

    // Read IMAP 3D file data into the solver database
    bool readIMAPwriteGDSII()
    {
        return true;
    }

    // Print the solver database and dump to SPEF file
    void printDump()
    {
        // Print
        cout << "Solver Database of IC Design, " << this->designName << ":" << endl;
        cout << " Waveforms:" << endl;
        (this->wf).print(); // Print the waveforms
        cout << " Parasitics in file " << this->outSPEF << ":" << endl;
        (this->para).print(); // Print the parasitics
        cout << "------" << endl;

        // Attempt to open output file
        ofstream spefFile((this->outSPEF).c_str());

        // Dump
        if (spefFile.is_open())
        {
            spef::Spef designPara = (this->para).toSPEF(this->designName, 1.0e-5);
            designPara.dump(spefFile);
        }
    }

    // Destructor
    ~SolverDataBase()
    {
        // Nothing
    }
};

#endif