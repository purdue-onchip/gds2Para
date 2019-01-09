/**
* File: solnoutclass.h
* Author: Michael R. Hayashi
* Date: 6 November 2018
* Desc: Header for solution and output custom classes
*/

// Make sure this file is included only once
#ifndef solnoutclass_h
#define solnoutclass_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <ctime>
#include <parser-spef/parser-spef.hpp>
#include <Eigen/Sparse>
#include "limboint.h"
#include "fdtd.h"
using namespace std;

// Define types for Eigen
typedef Eigen::Triplet<double, int> dTriplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> spMat;
struct myPruneFunctor
{
private:
    double reference;
public:
    // Default constructor
    myPruneFunctor()
    {
        this->reference = 0.;
    }

    // Parametrizer constructor
    myPruneFunctor(double ref)
    {
        this->reference = ref;
    }

    // Functor magic overloading parentheses
    inline bool operator() (const int& row, const int& col, const double& value) const
    {
        return (abs(value) > this->reference);
    }
};

// Custom classes for containing solution and output SPEF writer
class SimSettings
{
private:
    double lengthUnit;     // Units for lengths (m)
    vector<double> limits; // xmin, xmax, ymin, ymax, zmin, zmax (m)
    double freqUnit;       // Units for frequency (Hz)
    double freqScale;      // Frequency scaling?
    size_t nFreq;          // Number of frequencies in simulation
    vector<double> freqs;  // List of frequencies (linear or logarithmic [preferred] spacing)
public:
    // Default constructor
    SimSettings()
    {
        this->lengthUnit = 1.;
        this->limits = {0., 0., 0., 0., 0., 0.};
        this->freqUnit = 1.;
        this->freqScale = 0.;
        this->nFreq = 0;
        this->freqs = {};
    }

    // Parametrized constructor
    SimSettings(double lengthUnit, vector<double> limits, double freqUnit, double freqScale, vector<double> freqs)
    {
        this->lengthUnit = lengthUnit;
        if (limits.size() != 6)
        {
            cerr << "Must give minimum and maximum extents of design in vector of length 6. Defaulting to 0. to 0. for x, y, and z." << endl;
            this->limits = { 0., 0., 0., 0., 0., 0. };
        }
        else
        {
            vector<double> checkLims = limits;
            if (limits[0] > limits[1]) // xmin > xmax
            {
                checkLims[0] = limits[1];
                checkLims[1] = limits[0];
            }
            if (limits[2] > limits[3]) // ymin > ymax
            {
                checkLims[2] = limits[3];
                checkLims[3] = limits[2];
            }
            if (limits[4] > limits[5]) // zmin > zmax
            {
                checkLims[4] = limits[5];
                checkLims[5] = limits[4];
            }
            this->limits = checkLims;
        }
        this->freqUnit = freqUnit;
        this->freqScale = freqScale;
        this->nFreq = freqs.size();
        this->freqs = freqs;
    }

    // Get length unit (m)
    double getLengthUnit() const
    {
        return this->lengthUnit;
    }

    // Get limits (minimum extent, maximum extent) for x-dir, y-dir, and z-dir
    vector<double> getLimits() const
    {
        return this->limits;
    }

    // Get frequency unit (Hz)
    double getFreqUnit() const
    {
        return this->freqUnit;
    }

    // Get frequency scaling
    double getFreqScale() const
    {
        return this->freqScale;
    }

    // Get number of frequencies in simulation
    size_t getNFreq() const
    {
        return this->nFreq;
    }

    // Get list of frequencies (linear or logarithmic [preferred] spacing)
    vector<double> getFreqs() const
    {
        return this->freqs;
    }

    // Set length unit (m)
    void setLengthUnit(double lengthUnit)
    {
        this->lengthUnit = lengthUnit;
    }

    // Set limits: xmin, xmax, ymin, ymax, zmin, zmax (m)
    void setLimits(vector<double> limits)
    {
        if (limits.size() != 6)
        {
            cerr << "Must give minimum and maximum extents of design in vector of length 6. Defaulting to 0. to 0. for x, y, and z." << endl;
            this->limits = { 0., 0., 0., 0., 0., 0. };
        }
        else
        {
            vector<double> checkLims = limits;
            if (limits[0] > limits[1]) // xmin > xmax
            {
                checkLims[0] = limits[1];
                checkLims[1] = limits[0];
            }
            if (limits[2] > limits[3]) // ymin > ymax
            {
                checkLims[2] = limits[3];
                checkLims[3] = limits[2];
            }
            if (limits[4] > limits[5]) // zmin > zmax
            {
                checkLims[4] = limits[5];
                checkLims[5] = limits[4];
            }
            this->limits = checkLims;
        }
    }

    // Set frequency unit (Hz)
    void setFreqUnit(double freqUnit)
    {
        this->freqUnit = freqUnit;
    }

    // Set frequency scaling
    void setFreqScale(double freqScale)
    {
        this->freqScale = freqScale;
    }
    // Set list of frequencies (linear or logarithmic [preferred] spacing)
    void setFreqs(vector<double> freqs)
    {
        this->freqs = freqs;
        this->nFreq = freqs.size();
    }


    // Print the simulation settings
    void print() const
    {
        size_t numFreq = this->getNFreq();
        cout << " ------" << endl;
        cout << " Simulation Settings:" << endl;
        cout << "  Limits in x-direction: " << (this->limits)[0] << " m to " << (this->limits)[1] << " m" << endl;
        cout << "  Limits in y-direction: " << (this->limits)[2] << " m to " << (this->limits)[3] << " m" << endl;
        cout << "  Limits in z-direction: " << (this->limits)[4] << " m to " << (this->limits)[5] << " m" << endl;
        cout << "  List of " << numFreq << " frequencies to simulate with " << this->freqScale << " scaling:" << endl;
        for (size_t indi = 0; indi < numFreq; indi++)
        {
            if (numFreq - indi == 1) // Single frequency left to print
            {
                cout << "   #" << indi + 1 << " is " << (this->freqs)[indi] << " Hz" << endl;
            }
            else // Two frequencies per line
            {
                cout << "   #" << indi + 1 << " is " << (this->freqs)[indi] << " Hz, and #" << indi + 2 << " is " << (this->freqs)[indi + 1] << " Hz" << endl;
                indi++;
            }
        }
    }

    // Destructor
    ~SimSettings()
    {
        // Nothing
    }
};

class Layer
{
private:
    std::string layerName;        // Name of layer in physical stack-up
    int gdsiiNum;                 // Layer number in GDSII file
    double zStart;                // Z-coordinate of bottom of layer (m)
    double zHeight;               // Height of layer in z-direction (m)
    double epsilon_r;             // Relative permittivity of material
    double lossTan;               // Loss tangent of material
    double sigma;                 // (Real) Conductivity of material (S/m)
public:
    // Default constructor
    Layer()
    {
        this->layerName = "";
        this->gdsiiNum = -1;
        this->zStart = 0.;
        this->zHeight = 0.;
        this->epsilon_r = 1.;
        this->lossTan = 0.;
        this->sigma = 0.;
    }

    // Parametrized constructor
    Layer(std::string layerName, int gdsiiNum, double zStart, double zHeight, double epsilon_r, double lossTan, double sigma)
    {
        this->layerName = layerName;
        this->gdsiiNum = gdsiiNum;
        this->zStart = zStart;
        this->zHeight = zHeight;
        this->epsilon_r = epsilon_r;
        this->lossTan = lossTan;
        this->sigma = sigma;
    }

    // Get layer name
    std::string getLayerName() const
    {
        return this->layerName;
    }

    // Get GDSII file layer number
    // Metallic layers are nonnegative, bottom plane is 0, top plane is MAX, and -1 is used for dielectric, substrates, and undescribed planes
    int getGDSIINum() const
    {
        return this->gdsiiNum;
    }

    // Get layer bottom z-coordinate
    // Currently assigns -1.0 for layers with unknown absolute position
    double getZStart() const
    {
        return this->zStart;
    }

    // Get layer height
    double getZHeight() const
    {
        return this->zHeight;
    }

    // Get layer relative permittivity
    double getEpsilonR() const
    {
        return this->epsilon_r;
    }

    // Get layer loss tangent
    double getLossTan() const
    {
        return this->lossTan;
    }

    // Get layer conductivity
    double getSigma() const
    {
        return this->sigma;
    }

    // Set GDSII file layer number
    // Metallic layers are nonnegative, bottom plane is 0, top plane is MAX, and -1 is used for dielectric, substrates, and undescribed planes
    void setGDSIINum(int gdsiiNum)
    {
        this->gdsiiNum = gdsiiNum;
    }

    // Print the layer information
    void print() const
    {
        cout << " ------" << endl;
        cout << " Details for layer " << this->layerName << ":" << endl;
        if (this->gdsiiNum != -1)
        {
            cout << "  GDSII layer number: " << this->gdsiiNum << endl;
        }
        cout << "  Bottom z-coordinate: " << this->zStart << " m" << endl;
        cout << "  Layer height: " << this->zHeight << " m" << endl;
        cout << "  Relative permittivity: " << this->epsilon_r << endl;
        cout << "  Loss tangent: " << this->lossTan << endl;
        cout << "  Conductivity: " << this->sigma << " S/m" << endl;
        cout << " ------" << endl;
    }

    // Destructor
    ~Layer()
    {
        // Nothing
    }
};

class Port
{
private:
    std::string portName;        // Name of port
    char portDir;                // Direction of port
    double Z_source;             // Impedance of sourced attached to port (ohm)
    vector<double> coord;        // Supply and return coordinates: xsup, ysup, zsup, xret, yret, zret (m)
    int gdsiiNum;                // Layer number in GDSII file on which port exists
public:
    // Default constructor
    Port()
    {
        this->portName = "";
        this->portDir = 'B';
        this->Z_source = 0.;
        this->coord = { 0., 0., 0., 0., 0., 0. };
        this->gdsiiNum = -1;
    }

    // Parametrized constructor
    Port(std::string portName, char portDir, double Z_source, vector<double> coord, int gdsiiNum)
    {
        this->portName = portName;
        if ((portDir != 'I') && (portDir != 'O') && (portDir != 'B'))
        {
            cerr << "Port direction must be assigned as 'I' (input), 'O' (output), or 'B' (bidirectional). Defaulting to 'B'." << endl;
            this->portDir = 'B';
        }
        else
        {
            this->portDir = portDir;
        }
        this->Z_source = Z_source;
        if (coord.size() != 6)
        {
            cerr << "Must give supply then return coordinates in vector of length 6. Defaulting to origin for both." << endl;
            this->coord = { 0., 0., 0., 0., 0., 0. };
        }
        else
        {
            this->coord = coord;
        }
        this->gdsiiNum = gdsiiNum;
    }

    // Get port name
    std::string getPortName() const
    {
        return this->portName;
    }

    // Get port direction
    // ('I' = input, 'O' = output, 'B' = bidirectional)
    char getPortDir() const
    {
        return this->portDir;
    }

    // Get impedance of attached source
    double getZSource() const
    {
        return this->Z_source;
    }

    // Get supply and return coordinates
    // xsup, ysup, zsup, xret, yret, zret (m)
    vector<double> getCoord() const
    {
        return this->coord;
    }

    // Get GDSII layer number of port
    // Metallic layers are nonnegative, bottom plane is 0, top plane is MAX, and -1 is used for dielectric, substrates, and undescribed planes
    int getGDSIINum() const
    {
        return this->gdsiiNum;
    }

    // Set port name
    void setPortName(std::string name)
    {
        this->portName = name;
    }

    // Set port direction
    // ('I' = input, 'O' = output, 'B' = bidirectional)
    void setPortDir(char dir)
    {
        if ((dir != 'I') || (dir != 'O') || (dir != 'B'))
        {
            cerr << "Port direction must be assigned as 'I' (input), 'O' (output), or 'B' (bidirectional). Defaulting to 'B'." << endl;
            this->portDir = 'B';
        }
        else
        {
            this->portDir = dir;
        }
    }

    // Set impedance of attached source
    void setZSource(double Z_source)
    {
        this->Z_source = Z_source;
    }

    // Set supply and return coordinates
    // xsup, ysup, zsup, xret, yret, zret (m)
    void setCoord(vector<double> coord)
    {
        if (coord.size() != 6)
        {
            cerr << "Must give supply then return coordinates in vector of length 6. Defaulting to origin for both." << endl;
            this->coord = { 0., 0., 0., 0., 0., 0. };
        }
        else
        {
            this->coord = coord;
        }
    }

    // Set GDSII layer number on which port exists
    // Metallic layers are nonnegative, bottom plane is 0, top plane is MAX, and -1 is used for dielectric, substrates, and undescribed planes
    void setGDSIINum(int gdsiiNum)
    {
        this->gdsiiNum = gdsiiNum;
    }

    // Print the layer information
    void print() const
    {
        cout << "  ------" << endl;
        cout << "  Details for port " << this->portName << ":" << endl;
        cout << "   Port direction: ";
        switch (this->portDir)
        {
        case 'O':
            cout << "Output" << endl;
            break;
        case 'I':
            cout << "Input" << endl;
            break;
        case 'B':
            cout << "Bidirectional" << endl;
            break;
        default:
            cout << "Bidirectional" << endl; // Treat as bidirectional if direction unclear
        }
        cout << "   Attached source impedance: " << this->Z_source << " ohm" << endl;
        cout << "   Supply coordinates: (" << (this->coord)[0] << ", " << (this->coord)[1] << ", " << (this->coord)[2] << ") m" << endl;
        cout << "   Return coordinates: (" << (this->coord)[3] << ", " << (this->coord)[4] << ", " << (this->coord)[5] << ") m" << endl;
        cout << "   GDSII layer number: " << this->gdsiiNum << endl;
        cout << "  ------" << endl;
    }

    // Destructor
    ~Port()
    {
        /*this->portName = "";
        this->portDir = 'B';
        this->Z_source = 0.;
        this->coord = { 0., 0., 0., 0., 0., 0. };*/
    }
};

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
    size_t nPorts;             // Number of ports
    vector<Port> ports;        // Vector of port information
    spMat matG;                // Conductance matrix (S)
    spMat matC;                // Capacitance matrix (F)
public:
    // Default constructor
    Parasitics()
    {
        vector<Port> ports;
        spMat emptMat;
        this->nPorts = 0;
        this->ports = ports;
        this->matG = emptMat;
        this->matC = emptMat;
    }

    // Parametrized constructor
    Parasitics(vector<Port> ports, spMat matG, spMat matC)
    {
        this->nPorts = ports.size();
        this->ports = ports;
        this->matG = matG;
        this->matC = matC;
    }

    // Get number of ports
    size_t getNPort() const
    {
        return this->nPorts;
    }

    // Get vector of port information
    vector<Port> getPorts() const
    {
        return this->ports;
    }

    // Get conductance matrix
    spMat getGMatrix() const
    {
        return this->matG;
    }

    // Get capacitance matrix
    spMat getCMatrix() const
    {
        return this->matC;
    }

    // Set vector of port information
    void setPorts(vector<Port> ports)
    {
        this->ports = ports;
        this->nPorts = ports.size();
    }

    // Set conductance matrix
    void setGMatrix(spMat matG)
    {
        this->matG = matG;
    }

    // Set capacitance matrix
    void setCMatrix(spMat matC)
    {
        this->matC = matC;
    }

    // Return node-to-ground conductance
    double getGNodeGround(size_t indNode) const
    {
        // Calculate by summing along given row
        return (this->matG).innerVector(indNode).sum();
    }

    // Return total conductance represented in matrix (sum of diagonal entries)
    double getGTotal() const
    {
        return (this->matG).diagonal().sum();
    }

    // Return node-to-ground capacitance
    double getCNodeGround(size_t indNode) const
    {
        // Calculate by summing along given row
        return (this->matC).innerVector(indNode).sum();
    }

    // Return total capacitance represented in matrix (sum of diagonal entries)
    double getCTotal() const
    {
        return (this->matC).diagonal().sum();
    }

    // Find index of port by name
    // Returns index past number of ports if not found
    size_t locatePortName(std::string name) const
    {
        size_t indPort;
        for (indPort = 0; indPort < this->getNPort(); indPort++)
        {
            if (name.compare((this->ports[indPort]).getPortName()) == 0)
            {
                return indPort;
            }
        }
        return indPort;
    }

    // Return a port
    Port getPort(size_t indPort) const
    {
        return (this->ports)[indPort];
    }

    // Return all port names
    vector<string> findPortNames() const
    {
        vector<string> names;
        for (size_t indi = 0; indi < this->getNPort(); indi++)
        {
            names.push_back(((this->ports)[indi]).getPortName());
        }
        return names;
    }

    // Find ports in GDSII textboxes
    vector<Port> setPortsGDSII(size_t indCell, vector<double> center, AsciiDataBase adb)
    {
        // Recursively search cells for textboxes
        GeoCell thisCell = adb.getCell(indCell);
        vector<Port> portList;
        for (size_t indi = 0; indi < thisCell.getNumSRef(); indi++) // Follow structure references
        {
            string srefName = (thisCell.sreferences)[indi].getSRefName();
            size_t indNextCell = adb.locateCell(srefName);
            vector<double> thisSRef = (thisCell.sreferences)[indi].getSRefs();
            vector<Port> newPorts = this->setPortsGDSII(indNextCell, { (center[0] + thisSRef[0]), (center[1] + thisSRef[1]) }, adb); // Recursion step
            portList.insert(portList.end(), newPorts.begin(), newPorts.end());
        }
        for (size_t indi = 0; indi < thisCell.getNumText(); indi++) // Search cell at this level for textboxes
        {
            textbox thisTextBox = thisCell.textboxes[indi];
            vector<double> coords = {(thisTextBox.getTexts()[0] + center[0]), (thisTextBox.getTexts()[1] + center[1]), 0., (thisTextBox.getTexts()[0] + center[0]), 0., 0.}; // Assume: xret = xsup, yret = 0, zsup = zret = 0 until layer heights set
            portList.emplace_back(Port(thisTextBox.getTextStr(), 'B', 50.0, thisTextBox.getTexts(), thisTextBox.getLayer()));
        }

        // Set vector of port information
        this->ports = portList;
    }

    // Print the parasitics information
    void print() const
    {
        int numPort = getNPort();
        cout << " ------" << endl;
        cout << " Parasitics Details:" << endl;
        cout << "  List of " << numPort << " ports:" << endl;
        for (size_t indi = 0; indi < numPort; indi++) // Handle each port name
        {
            (this->ports)[indi].print();
        }
        cout << "  Conductance Matrix (S):" << endl;
        for (size_t indi = 0; indi < (this->matG).outerSize(); indi++)
        {
            for (spMat::InnerIterator it(this->matG, indi); it; ++it)
            {
                cout << "   row " << setw(4) << it.row() << ", column " << setw(4) << it.col() << ", value " << it.value() << endl;
            }
        }
        cout << "  Capacitance Matrix (F):" << endl;
        for (size_t indi = 0; indi < (this->matC).outerSize(); indi++)
        {
            for (spMat::InnerIterator it(this->matC, indi); it; ++it)
            {
                cout << "   row " << setw(4) << it.row() << ", column " << setw(4) << it.col() << ", value " << it.value() << endl;
            }
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
            para.name_map.emplace(indi + 1, (this->ports)[indi].getPortName()); // Create name map for each port
            para.ports.emplace_back("*" + to_string(indi + 1)); // Instantiate and push new port entry by name map
            switch ((this->ports)[indi].getPortDir()) // Assign port direction
            {
            case 'O':
                para.ports.back().direction = spef::ConnectionDirection::OUTPUT;
                break;
            case 'I':
                para.ports.back().direction = spef::ConnectionDirection::INPUT;
                break;
            case 'B':
                para.ports.back().direction = spef::ConnectionDirection::INOUT;
                break;
            default:
                para.ports.back().direction = spef::ConnectionDirection::INOUT; // Treat as bidirectional if direction unclear
            }
        }

        // Populate Spef struct nets vector (single net)
        double capTot = this->getCTotal();
        double condTot = this->getGTotal();
        (this->matC).prune(myPruneFunctor(saveThresh * capTot)); // Prune away nonzeros sufficiently smaller than threshold
        (this->matG).prune(myPruneFunctor(saveThresh * condTot));
        para.nets.emplace_back(spef::Net());
        para.nets.back().name = "all";
        para.nets.back().lcap = capTot;

        //spMat symC = ((this->matC) + (this->matC).transpose()) / 2.0;
        //spMat symG = ((this->matG) + (this->matG).transpose()) / 2.0;
        for (size_t indi = 0; indi < numPort; indi++)
        {
            para.nets.back().connections.emplace_back(spef::Connection());
            para.nets.back().connections.back().name = (this->ports)[indi].getPortName();
            para.nets.back().connections.back().type = spef::ConnectionType::EXTERNAL; // All ports are external connections, not internal to cells
            para.nets.back().connections.back().direction = (para.ports[indi]).direction; // Same as port direction
            for (spMat::InnerIterator it(this->matC, indi); it; ++it)
            {
                para.nets.back().caps.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[it.col()].name, it.value())); 
                /*if ((it.row() == it.col()) && (abs(this->getCNodeGround(indi)) >= saveThresh * capTot)) // Diagonal element, so sufficient admittance to ground
                {
                    para.nets.back().caps.emplace_back(forward_as_tuple((para.ports)[indi].name, "", this->getCNodeGround(indi)));
                }
                else if (it.col() > it.row()) // Upper triangular element, so negate node-to-node admittance
                {
                    para.nets.back().caps.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[it.col()].name, -it.value()));
                } // Skip lower triangular element, which is repeated */
            }
            for (spMat::InnerIterator it(this->matG, indi); it; ++it)
            {
                para.nets.back().ress.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[it.col()].name, 1.0 / (it.value())));
                /*if ((it.row() == it.col()) && (abs(this->getGNodeGround(indi)) >= saveThresh * condTot)) // Diagonal element, so sufficient admittance to ground
                {
                    para.nets.back().ress.emplace_back(forward_as_tuple((para.ports)[indi].name, "", 1.0 / this->getGNodeGround(indi)));
                }
                else if (it.col() > it.row()) // Upper triangular element, so negate node-to-node admittance
                {
                    para.nets.back().ress.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[it.col()].name, -1.0 / it.value()));
                } // Skip lower triangular element, which is repeated */
            }
        }

        // Return the Spef struct
        return para;
    }

    // Destructor
    ~Parasitics()
    {
        vector<Port> ports = {};
        spMat emptMat;
        this->nPorts = 0;
        this->ports = ports;
        this->matG = emptMat;
        this->matC = emptMat;
    }
};

struct SolverDataBase
{
private:
    std::string designName; // Name of design
    SimSettings settings;   // Simulation settings
    vector<Layer> layers;   // Layer stack-up information
    Waveforms wf;           // Waveforms
    Parasitics para;        // Parasitics
    std::string outSPEF;    // SPEF output file name
public:
    // Default constructor
    SolverDataBase()
    {
        SimSettings settings;
        vector<Layer> layers;
        Waveforms wf;
        Parasitics para;
        this->designName = "";
        this->settings = settings;
        this->layers = layers;
        this->wf = wf;
        this->para = para;
        this->outSPEF = "";
    }

    // Parametrized constructor
    SolverDataBase(std::string designName, Waveforms wf, Parasitics para)
    {
        SimSettings settings;
        vector<Layer> layers;
        this->designName = designName;
        this->settings = settings;
        this->layers = layers;
        this->wf = wf;
        this->para = para;
        this->outSPEF = "";
    }

    // Get name of design
    std::string getDesignName() const
    {
        return this->designName;
    }

    // Get simulation settings
    SimSettings getSimSettings() const
    {
        return this->settings;
    }

    // Get number of layers
    size_t getNumLayer() const
    {
        return (this->layers).size();
    }

    // Get vector of layer information
    vector<Layer> getLayers() const
    {
        return this->layers;
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

    // Set name of design
    void setDesignName(std::string designName)
    {
        this->designName = designName;
    }

    // Set simulation settings
    void setSimSettings(SimSettings newSettings)
    {
        this->settings = newSettings;
    }

    // Set layer stackup
    void setLayers(vector<Layer> layers)
    {
        this->layers = layers;
    }

    // Set waveforms
    void setWaveforms(Waveforms wf)
    {
        this->wf = wf;
    }

    // Set parasitics
    void setParasitics(Parasitics para)
    {
        this->para = para;
    }

    // Set SPEF output file name
    void setOutSPEF(std::string fileName)
    {
        this->outSPEF = fileName;
    }

    // Find index of layer by name
    // Returns index past number of layers if not found
    size_t locateLayerName(std::string name) const
    {
        size_t indLayer;
        for (indLayer = 0; indLayer < (this->layers).size(); indLayer++)
        {
            if (name.compare((this->layers[indLayer]).getLayerName()) == 0)
            {
                return indLayer;
            }
        }
        return indLayer;
    }

    // Find index of layer by GDSII number
    // Returns index past number of layers if not found
    size_t locateLayerGDSII(int gdsiiNum) const
    {
        size_t indLayer;
        for (indLayer = 0; indLayer < (this->layers).size(); indLayer++)
        {
            if ((this->layers[indLayer]).getGDSIINum() == gdsiiNum)
            {
                return indLayer;
            }
        }
        return indLayer;
    }

    // Find index of layer by starting z-coordinate
    // Returns index past number of layers if not found
    size_t locateLayerZStart(double zStart) const
    {
        size_t indLayer;
        for (indLayer = 0; indLayer < (this->layers).size(); indLayer++)
        {
            // See if layer start coordinate is within 1% of the height of that layer
            if (abs((this->layers[indLayer]).getZStart() - zStart) < 0.01 * (this->layers[indLayer]).getZHeight())
            {
                return indLayer;
            }
        }
        return indLayer;
    }

    // Return a layer
    Layer getLayer(size_t indLayer) const
    {
        return (this->layers)[indLayer];
    }

    // Return all layer names
    vector<string> findLayerNames() const
    {
        vector<string> names;
        for (size_t indi = 0; indi < (this->layers).size(); indi++)
        {
            names.push_back(((this->layers)[indi]).getLayerName());
        }
        return names;
    }

    // Read interconnect modeling platform (IMP) file data into the solver database
    bool readIMPwriteGDSII(std::string impFileName, std::string gdsiiFileName)
    {
        // Attempt to open .imp file
        ifstream impFile(impFileName.c_str());
        if (impFile.is_open())
        {
            // Time manipulations
            char timeStr[80]; // Character array to hold formatted time
            time_t rawtime; // Timer variable
            strftime(timeStr, sizeof(timeStr), "%d-%m-%Y %H:%M:%S", localtime(&rawtime)); // Use local timezone to format string
            std::string time(timeStr); // Formatted time to string parametrized constructor

            // Build single geometric cell from limboint.h to store information
            GeoCell cellIMP;
            cellIMP.cellName = impFileName.substr(0, impFileName.find(".", 1));
            cellIMP.dateCreate = time;
            cellIMP.dateMod = time;

            // Build ASCII database from limboint.h
            AsciiDataBase adbIMP;
            adbIMP.setFileName(gdsiiFileName);
            adbIMP.setDateMod(time);
            adbIMP.setDateAccess(time);

            // File is readable line-by-line
            std::string fileLine;
            getline(impFile, fileLine);

            // Save opening lines as metadata (mostly ignored right now)
            if (fileLine.compare(0, 7, "IMAP3D ") == 0)
            {
                // Save interconnect modeling platform version number
                std::string version = fileLine.substr(7, fileLine.length() - 7); // 7 characters in header syntax
            }

            // Read rest of file line-by-line
            while (!impFile.eof())
            {
                // Initialize file-scope varialbes
                static double stripLength = 0; // Initialize the strip length
                static int maxGDSIILayer = 0; // Maximum GDSII layer number encountered
                static int numFreqPts = 0; // Number of points in frequency sweep
                static vector<double> freqList = {}; // List of frequencies in sweep

                // Handle units
                if (fileLine.compare(0, 12, "LINEARUNITS ") == 0)
                {
                    // Extract unit text
                    std::string units = fileLine.substr(12, fileLine.length() - 13);

                    // Find SI multiplier for given units
                    double multSI = 1.0;
                    if (units.compare("ym") == 0) { multSI = 1.e-24; }
                    else if (units.compare("zm") == 0) { multSI = 1.e-21; }
                    else if (units.compare("am") == 0) { multSI = 1.e-18; }
                    else if (units.compare("fm") == 0) { multSI = 1.e-15; }
                    else if (units.compare("pm") == 0) { multSI = 1.e-12; }
                    else if (units.compare("nm") == 0) { multSI = 1.e-09; }
                    else if (units.compare("um") == 0) { multSI = 1.e-06; }
                    else if (units.compare("mm") == 0) { multSI = 1.e-03; }
                    else if (units.compare("cm") == 0) { multSI = 1.e-02; }
                    else if (units.compare("dm") == 0) { multSI = 1.e-01; }

                    // Propagate units information to ASCII database now
                    adbIMP.setdbUserUnits(1.);
                    adbIMP.setdbUnits(multSI);
                }
                // Handle design name and skip to analysis parameters
                else if (fileLine.compare(0, 5, "NAME ") == 0)
                {
                    this->designName = fileLine.substr(5, fileLine.length() - 6);

                    // Save file position and jump ahead to analysis section
                    int savedFilePos = impFile.tellg();
                    while (!impFile.eof())
                    {
                        // Keep reading new lines in the file
                        getline(impFile, fileLine);

                        // Handle frequency information
                        if (fileLine.compare(0, 9, "Frequency") == 0)
                        {
                            // Find frequency sweep delimiters
                            size_t indBegin = fileLine.find("begin=");
                            size_t indEnd = fileLine.find("end=");
                            size_t indNoP = fileLine.find("numberofpoints=");

                            // Save frequency sweep information to variables
                            double freqBegin = stod(fileLine.substr(indBegin + 6, indEnd - indBegin - 7)); // Length is index difference minus space
                            double freqEnd = stod(fileLine.substr(indEnd + 4, indNoP - indEnd - 5));
                            numFreqPts = stoi(fileLine.substr(indNoP + 15));
                            if (numFreqPts == 1)
                            {
                                freqList.push_back(freqBegin);
                            }
                            else if (numFreqPts == 2)
                            {
                                freqList.push_back(freqBegin);
                                freqList.push_back(freqEnd);
                            }
                            else
                            {
                                double exp10Step = log10(freqEnd / freqBegin) / (numFreqPts - 1);
                                freqList.push_back(freqBegin);
                                for (size_t indi = 1; indi < numFreqPts - 1; indi++)
                                {
                                    freqList.push_back(freqList.back() * pow(10, exp10Step));
                                }
                                freqList.push_back(freqEnd); // Ensure last frequency is exact
                            }

                            // Wait to push back simulation settings after conductors are saved (need extrema of coordinates)
                        }

                        // Handle length in this weird spot
                        if (fileLine.compare(0, 6, "Length") == 0)
                        {
                            // Finally can save stripline length
                            stripLength = stod(fileLine.substr(7)) * adbIMP.getdbUnits();
                        }
                    }

                    // Jump back to line after name section
                    impFile.clear(); // Clear the eofbit and goodbit flags
                    impFile.seekg(savedFilePos);
                    getline(impFile, fileLine);
                }
                // Handle layer stack
                if (fileLine.compare(0, 5, "STACK") == 0)
                {
                    // Move down one line
                    getline(impFile, fileLine);

                    // Keep reading until end of layer stack
                    while ((fileLine.compare(0, 10, "CONDUCTORS") != 0) && (fileLine.compare(0, 8, "BOUNDARY") != 0))
                    {
                        // Record each layer
                        if (fileLine.length() >= 3)
                        {
                            // Find layer property delimiters
                            size_t indLayNam = fileLine.find(" ");
                            size_t indZStart = fileLine.find("z=");
                            size_t indHeight = fileLine.find("h=");
                            size_t indPermit = fileLine.find("e=");
                            size_t indLosTan = fileLine.find("TanD=");
                            size_t indConduc = fileLine.find("sigma=");

                            // Save layer information to variables
                            std::string layerName = fileLine.substr(0, indLayNam);
                            int gdsiiNum = -1;
                            if (layerName.find('M') != string::npos)
                            {
                                size_t indNumber = layerName.find("M");
                                gdsiiNum = stoi(layerName.substr(indNumber + 1, layerName.length() - indNumber - 1));
                                if (gdsiiNum > maxGDSIILayer) { maxGDSIILayer = gdsiiNum; }
                            }
                            double zStart = 0.;
                            if (indZStart != string::npos)
                            {
                                zStart = stod(fileLine.substr(indZStart + 2, indHeight - indZStart - 3)) * adbIMP.getdbUnits();
                            }
                            double zHeight = stod(fileLine.substr(indHeight + 2, fileLine.find(" ", indHeight) - indHeight - 2)) * adbIMP.getdbUnits();
                            double epsilon_r = 1.;
                            if (indPermit != string::npos)
                            {
                                epsilon_r = stod(fileLine.substr(indPermit + 2, fileLine.find(" ", indPermit) - indPermit - 2));
                            }
                            double lossTan = 0.;
                            if (indLosTan != string::npos) // See if it handles being at the end of line
                            {
                                lossTan = stod(fileLine.substr(indLosTan + 5, fileLine.find(" ", indLosTan) - indLosTan - 5));
                            }
                            double sigma = 0.;
                            if (indConduc != string::npos) // See if it handles being at the end of line
                            {
                                sigma = stod(fileLine.substr(indConduc + 6, fileLine.find(" ", indConduc) - indConduc - 6));
                            }

                            // Push new layer to class vector
                            (this->layers).emplace_back(Layer(layerName, gdsiiNum, zStart, zHeight, epsilon_r, lossTan, sigma));
                        }
                        // Keep moving down the layer stack
                        getline(impFile, fileLine);
                    }
                }
                // Handle conductors
                if (fileLine.compare(0, 10, "CONDUCTORS") == 0)
                {
                    // Move down one line
                    getline(impFile, fileLine);

                    // Keep reading until end of conductor list
                    int gdsiiNum = 1; // Layer number of target GDSII file
                    double xmin = 0.;
                    double xmax = 0.;
                    double ymin = 0.;
                    double ymax = 0.;
                    while ((fileLine.compare(0, 8, "BOUNDARY") != 0) && (fileLine.compare(0, 9, "PORTTABLE") != 0))
                    {
                        // Record each conductor as a box
                        if (fileLine.length() >= 3)
                        {
                            // Find conductor property delimiters
                            size_t indCategory = fileLine.find(" ", 0);
                            size_t indCondName = fileLine.find(" ", indCategory);
                            size_t indX1 = fileLine.find("x1=", indCondName);
                            size_t indY1 = fileLine.find("y1=", indCondName);
                            size_t indZ1 = fileLine.find("z1=", indCondName);
                            size_t indX2 = fileLine.find("x2=", indCondName);
                            size_t indY2 = fileLine.find("y2=", indCondName);
                            size_t indZ2 = fileLine.find("z2=", indCondName);
                            size_t indSigma = fileLine.find("sigma=", indCondName);
                            size_t indLayer = fileLine.find("layer=", indCondName);
                            size_t indGroup = fileLine.find("group=", indCondName);

                            // Save conductor information to variables
                            std::string category = fileLine.substr(0, indCategory);
                            std::string condName = fileLine.substr(indCategory + 1, indCondName - indCategory - 1);
                            double x1 = stod(fileLine.substr(indX1 + 3, indY1 - indX1 - 4)) * adbIMP.getdbUnits(); // Length is index difference minus space
                            double y1 = stod(fileLine.substr(indY1 + 3, indZ1 - indY1 - 4)) * adbIMP.getdbUnits();
                            double x2 = stod(fileLine.substr(indX2 + 3, indY2 - indX2 - 4)) * adbIMP.getdbUnits();
                            double y2 = stod(fileLine.substr(indY2 + 3, indZ2 - indY2 - 4)) * adbIMP.getdbUnits();
                            std::string sigma = fileLine.substr(indSigma, indLayer - indSigma - 1);
                            std::string group;
                            if (indGroup != string::npos)
                            {
                                group = fileLine.substr(indGroup + 6);
                            }

                            // Check for extrema of coordinates
                            if (x1 < xmin) { xmin = x1; }
                            if (x2 > xmax) { xmax = x2; }
                            if (y1 < ymin) { ymin = y1; }
                            if (y2 > ymax) { ymax = y2; }

                            // Assign layer number based on layer stack-up
                            size_t indThisLayer = this->locateLayerName(fileLine.substr(indLayer + 6, fileLine.find(" ", indLayer)));
                            if (indThisLayer < (this->layers).size())
                            {
                                // Layer name for conductor matched name in layer stack-up
                                gdsiiNum = ((this->layers)[indThisLayer]).getGDSIINum();
                            }

                            // Push new box to the geometric cell
                            cellIMP.boxes.emplace_back(box({ x2, y1, x2, y2 + stripLength, x1, y2 + stripLength, x1, y1, x2, y1 }, gdsiiNum, { sigma, condName, category, group }, 0));

                        }
                        // Keep moving down the conductor list
                        getline(impFile, fileLine);
                    }

                    // Add on planes
                    size_t indGroundPlane = this->locateLayerName("GroundPlane"); // Start with ground (bottom) plane
                    ((this->layers)[indGroundPlane]).setGDSIINum(0); // Assign ground plane to layer 0
                    vector<std::string> propGround; // Ground plane properties
                    propGround.push_back("sigma=" + to_string(((this->layers)[indGroundPlane]).getSigma()));
                    propGround.push_back("GroundPlane");
                    propGround.push_back("plane");
                    propGround.push_back("");
                    cellIMP.boxes.emplace_back(box({ xmax, ymin, xmax, 2 * ymax + stripLength, xmin, 2 * ymax + stripLength, xmin, ymin, xmax, ymin }, ((this->layers)[indGroundPlane]).getGDSIINum(), propGround, 0));
                    size_t indTopPlane = this->locateLayerName("TopPlane"); // Next is top plane
                    ((this->layers)[indTopPlane]).setGDSIINum(++maxGDSIILayer); // Assign top plane to layer one more than maximum metallic layer
                    vector<std::string> propTop; // Top plane properties
                    propTop.push_back("sigma=" + to_string(((this->layers)[indTopPlane]).getSigma()));
                    propTop.push_back("TopPlane");
                    propTop.push_back("plane");
                    propTop.push_back("");
                    cellIMP.boxes.emplace_back(box({ xmax, ymin, xmax, 2 * ymax + stripLength, xmin, 2 * ymax + stripLength, xmin, ymin, xmax, ymin }, ((this->layers)[indTopPlane]).getGDSIINum(), propTop, 0));

                    // Save simulation settings at last
                    this->settings = SimSettings(adbIMP.getdbUnits(), { xmin, xmax, ymin, 2 * ymax + stripLength, (this->layers).back().getZStart(), (this->layers).front().getZStart() + (this->layers).front().getZHeight() }, 1.0, 0.0, freqList);
                }
                // Handle port table
                if (fileLine.compare(0, 9, "PORTTABLE") == 0)
                {
                    // Move down one line
                    getline(impFile, fileLine);

                    // Keep reading until end of port table
                    vector<Port> ports;
                    while ((fileLine.compare(0, 8, "ANALYSIS") != 0) && !impFile.eof())
                    {
                        if (fileLine.length() >= 3)
                        {
                            // Find port table delimiters
                            size_t indNameStart = fileLine.find("group=", 0) + 6;
                            size_t indNameEnd = fileLine.find(" ", indNameStart);
                            size_t indZNearStart = fileLine.find("znear=", indNameEnd) + 6;
                            size_t indZNearEnd = fileLine.find(" ", indZNearStart);
                            size_t indZFar = fileLine.find("zfar=", indNameEnd) + 5;

                            // Save port table information to variables
                            std::string groupName = fileLine.substr(indNameStart, indNameEnd - indNameStart);
                            double Z_near = stod(fileLine.substr(indZNearStart, indZNearStart - indZNearEnd));
                            double Z_far = stod(fileLine.substr(indZFar));

                            // Register a new port in vectors
                            vector<double> portSupRet;
                            size_t indCond, indLayer;
                            for (indCond = 0; indCond < cellIMP.boxes.size(); indCond++)
                            {
                                if (groupName.compare(((cellIMP.boxes[indCond]).getProps())[1]) == 0)
                                {
                                    break;
                                }
                            }
                            if (indCond >= cellIMP.boxes.size())
                            {
                                portSupRet = { 0., 0., 0., 0., 0., 0. }; // Not found, default has supply and return at origin
                            }
                            else
                            {
                                vector<double> boxCoord = (cellIMP.boxes[indCond]).getBoxes(); // Coordinates of box on layer
                                indLayer = this->locateLayerGDSII((cellIMP.boxes[indCond]).getLayer()); // Index of layer in structure field
                                portSupRet = { boxCoord[0], boxCoord[1], (this->layers)[indLayer].getZStart(), boxCoord[0], boxCoord[1], 0. }; // xsup = LR xcoord, ysup = LR ycoord, zsup = layer zStart, xret = xsup, yret = ysup, zret = 0.
                            }
                            ports.emplace_back(Port(groupName, 'B', Z_near, portSupRet, (this->layers)[indLayer].getGDSIINum()));
                        }
                        // Keep moving down the port table
                        getline(impFile, fileLine);
                    }

                    // Propagate port information to Solver Database now
                    this->para = Parasitics(ports, spMat(), spMat());
                }

                // Keep reading new lines in file
                getline(impFile, fileLine);
            }

            // Close file
            impFile.close();

            // Update ASCII database
            adbIMP.setLibName(this->designName);
            adbIMP.appendCell(cellIMP);
            adbIMP.setdbUnits(adbIMP.getdbUnits() * 1.e-3); // Rescale .imp file units 0.001x to allow integer representation in GDSII

            // Print the ASCII database
            adbIMP.print({ 0 });
            (this->settings).print();

            // Write GDSII file to hard drive
            bool dumpPassed = adbIMP.dump();
            return dumpPassed;
        }
        else
        {
            // File could not be opened
            return false;
        }
    }

    // Load simulation input file 
    bool readSimInput(std::string inputFileName)
    {
        // Attempt to open simulation input file
        ifstream inputFile(inputFileName.c_str());
        if (inputFile.is_open())
        {
            // File is readable line-by-line
            std::string fileLine;
            getline(inputFile, fileLine);

            // Read rest of file line-by-line
            while (!inputFile.eof())
            {
                // Handle total size
                if (fileLine.compare(0, 10, "TOTAL SIZE") == 0)
                {
                    // Move down one line
                    getline(inputFile, fileLine);
                    
                    // Obtain limits of IC design size
                    double xmin, xmax, ymin, ymax, zmin, zmax;
                    size_t indCoordStart = 0;
                    size_t indCoordEnd = fileLine.find(" ", indCoordStart);
                    xmin = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                    indCoordStart = indCoordEnd + 1;
                    indCoordEnd = fileLine.find(" ", indCoordStart);
                    xmax = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                    indCoordStart = indCoordEnd + 1;
                    indCoordEnd = fileLine.find(" ", indCoordStart);
                    ymin = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                    indCoordStart = indCoordEnd + 1;
                    indCoordEnd = fileLine.find(" ", indCoordStart);
                    ymax = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                    indCoordStart = indCoordEnd + 1;
                    indCoordEnd = fileLine.find(" ", indCoordStart);
                    zmin = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                    indCoordStart = indCoordEnd + 1;
                    zmax = stod(fileLine.substr(indCoordStart));
                    
                    // Move down one line
                    getline(inputFile, fileLine);
                    
                    // Obtain length unit
                    double lengthUnit = stod(fileLine.substr(fileLine.find("lengthUnit = ") + 13, fileLine.find(" #")));
                    
                    // Create simulation settings
                    this->settings = SimSettings(lengthUnit, { xmin * lengthUnit, xmax * lengthUnit, ymin * lengthUnit, ymax * lengthUnit, zmin * lengthUnit, zmax * lengthUnit }, 1., 0., {});
                }
                // Handle frequency sweep parameters
                else if (fileLine.compare(0, 9, "FREQUENCY") == 0)
                {
                    // Move down one line
                    getline(inputFile, fileLine);
                    
                    // Obtain frequency unit
                    double freqUnit = stod(fileLine.substr(fileLine.find("freqUnit = ") + 11, fileLine.find(" #")));
                    
                    // Move down one line
                    getline(inputFile, fileLine);
                    
                    // Obtain first frequency
                    double freqStart = stod(fileLine.substr(fileLine.find("freqStart = ") + 12, fileLine.find(" #")));
                    
                    // Move down one line
                    getline(inputFile, fileLine);
                    
                    // Obtain last frequency
                    double freqEnd = stod(fileLine.substr(fileLine.find("freqEnd = ") + 10, fileLine.find(" #")));
                    
                    // Move down one line
                    getline(inputFile, fileLine);
                    
                    // Obtain number of frequencies
                    size_t nFreq = stoi(fileLine.substr(fileLine.find("nfreq = ") + 8, fileLine.find(" #")));
                    
                    // Move down one line
                    getline(inputFile, fileLine);
                    
                    // Obtain frequency scaling
                    double freqScale = stod(fileLine.substr(fileLine.find("freqScale = ") + 12, fileLine.find(" #")));

                    // Establish frequency list
                    vector<double> freqList;
                    if (nFreq == 1)
                    {
                        freqList.push_back(freqStart);
                    }
                    else if (nFreq == 2)
                    {
                        freqList.push_back(freqStart);
                        freqList.push_back(freqEnd);
                    }
                    else
                    {
                        double exp10Step = log10(freqEnd / freqStart) / (nFreq - 1);
                        freqList.push_back(freqStart);
                        for (size_t indi = 1; indi < nFreq - 1; indi++)
                        {
                            freqList.push_back(freqList.back() * pow(10, exp10Step));
                        }
                        freqList.push_back(freqEnd); // Ensure last frequency is exact
                    }
                    
                    // Update simulation settings
                    (this->settings).setFreqUnit(freqUnit);
                    (this->settings).setFreqScale(freqScale);
                    (this->settings).setFreqs(freqList);
                }
                // Handle dielectric stack-up
                else if (fileLine.compare(0, 16, "DIELECTRIC STACK") == 0)
                {
                    // Move down one line
                    getline(inputFile, fileLine);

                    // Obtain number of dieletric layers in stack-up
                    size_t numStack = stoi(fileLine.substr(fileLine.find("numStack = ") + 11, fileLine.find(" #")));

                    // Move down one line
                    getline(inputFile, fileLine);

                    // Read each line in dielectric stack
                    for (size_t indStack = 0; indStack < numStack; indStack++)
                    {
                        // Get dielectric stack delimiters
                        size_t indNameEnd = fileLine.find(" ");
                        size_t indHeight = fileLine.find("h = ") + 4;
                        size_t indRelPermit = fileLine.find("e = ") + 4;
                        size_t indComment = fileLine.find(" #");

                        // Save dielectric stack information to variables
                        std::string layerName = fileLine.substr(0, indNameEnd);
                        double layerHeight = stod(fileLine.substr(indHeight, indRelPermit - indHeight - 5));
                        double layerEpsilonR = stod(fileLine.substr(indRelPermit, indComment - indRelPermit));

                        // Register a new layer in layers field
                        (this->layers).emplace_back(Layer(layerName, -1, -1.0, layerHeight, layerEpsilonR, 0.0, 0.0));

                        // Keep moving down the dielectric stack
                        getline(inputFile, fileLine);
                    }
                }
                // Handle port list if not already populated by GDSII file
                else if ((fileLine.compare(0, 4, "PORT") == 0) && ((this->para).getPorts().size() == 0))
                {
                    // Move down one line
                    getline(inputFile, fileLine);

                    // Obtain number of ports in list
                    size_t numPort = 0;
                    if (fileLine.find("numPorts = ") < string::npos)
                    {
                        numPort = stoi(fileLine.substr(fileLine.find("numPorts = ") + 11, fileLine.find(" #")));
                    }
                    else if (fileLine.find("numPort = ") < string::npos)
                    {
                        numPort = stoi(fileLine.substr(fileLine.find("numPort = ") + 10, fileLine.find(" #")));
                    }

                    // Move down one line
                    getline(inputFile, fileLine);

                    // Read each line in the port list
                    vector<Port> ports;
                    for (size_t indPort = 0; indPort < numPort; indPort++)
                    {
                        // Port information recorded differently based on how many numbers appear
                        size_t nSpace = count(fileLine.begin(), fileLine.end(), ' ');

                        // Obtain limits of IC design size
                        double xsup, ysup, zsup, xret, yret, zret;
                        int sourceDir;
                        char portDir;
                        int portLayer = -1;
                        if (nSpace == 4)
                        {
                            size_t indCoordStart = 0;
                            size_t indCoordEnd = fileLine.find(" ", indCoordStart);
                            xsup = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            ysup = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            xret = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            yret = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            portLayer = stoi(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart)); // Neglect comments
                            sourceDir = 0; // No soure direction information included

                            zsup = (this->layers)[this->locateLayerGDSII(portLayer)].getZStart();
                            zret = zsup + (this->layers)[this->locateLayerGDSII(portLayer)].getZHeight(); // Return assumed to be on same layer but on upper end
                        }
                        else if (nSpace >= 6)
                        {
                            size_t indCoordStart = 0;
                            size_t indCoordEnd = fileLine.find(" ", indCoordStart);
                            xsup = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            ysup = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            zsup = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            xret = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            yret = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            zret = stod(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart));
                            indCoordStart = indCoordEnd + 1;
                            indCoordEnd = fileLine.find(" ", indCoordStart);
                            sourceDir = stoi(fileLine.substr(indCoordStart, indCoordEnd - indCoordStart)); // Neglect comments
                            portLayer = (this->layers)[this->locateLayerZStart(zsup)].getGDSIINum();
                        }

                        // Append port information to vectors (all ports officially unnamed, so use number)
                        switch (sourceDir) // Assign port direction
                        {
                        case -1:
                            portDir = 'O';
                            break;
                        case +1:
                            portDir = 'I';
                            break;
                        case 0:
                            portDir = 'B';
                            break;
                        default:
                            portDir = 'B'; // Treat as bidirectional if direction unclear
                        }
                        ports.emplace_back(Port("port" + to_string(indPort + 1), portDir, 50., { xsup, ysup, zsup, xret, yret, zret }, portLayer));

                        // Move down one line
                        getline(inputFile, fileLine);
                    }

                    // Propagate port list to parasitics data structure now
                    this->para = Parasitics(ports, spMat(numPort, numPort), spMat(numPort, numPort));
                }

                // Keep reading new lines in file
                getline(inputFile, fileLine);
            }

            // Close file
            inputFile.close();
            return true;
        }
        else
        {
            // File could not be opened
            return false;
        }
    }

    // Convert to fdtdMesh
    fdtdMesh convertToFDTDMesh(int numCdtRow)
    {
        // Initialize the output
        unordered_set<double> portCoorx, portCoory;
        fdtdMesh data;
        data.numCdtRow = numCdtRow; // Use the sole argument

        // Use simulation settings to set fields
        data.lengthUnit = (this->settings).getLengthUnit();
        data.freqUnit = (this->settings).getFreqUnit();
        data.nfreq = (this->settings).getNFreq();
        data.freqStart = (this->settings).getFreqs().front();
        data.freqEnd = (this->settings).getFreqs().back();

        // Use layer stack-up information to set fields
        data.numStack = this->getNumLayer();
        data.stackEps = (double*)malloc(sizeof(double) * data.numStack);
        data.stackBegCoor = (double*)malloc(sizeof(double) * data.numStack);
        data.stackEndCoor = (double*)malloc(sizeof(double) * data.numStack);
        //data.stackEpsn; // Is this needed from SolverDataBase?
        for (size_t indi = 0; indi < data.numStack; indi++)
        {
            Layer thisLayer = this->getLayer(indi); // Get copy of layer for this iteration
            data.stackEps[indi] = thisLayer.getEpsilonR(); // See if relative or absolute permittivity!!!
            data.stackBegCoor[indi] = thisLayer.getZStart();
            data.stackEndCoor[indi] = thisLayer.getZStart() + thisLayer.getZHeight();
            data.stackName.push_back(thisLayer.getLayerName());
        }

        // Set conductor information from saved stack information
        data.stackCdtMark = (double*)malloc(sizeof(double) * data.numStack);
        for (size_t indi = 0; indi < data.numCdtRow; indi++)
        {
            for (size_t indj = 0; indj < data.numStack; indj++)
            {
                // See if conductor layer numbers being built correspond to existing layer name
                /*if (to_string(data.conductorIn[indi].layer).compare(data.stackName[indj]) == 0)
                {
                    data.conductorIn[indi].zmin = data.stackBegCoor[indj];
                    data.conductorIn[indi].zmax = data.stackEndCoor[indj];
                    data.stackCdtMark[indj] = 1;
                }*/

                // See if conductor layer numbers being built correspond to GDSII layer number
                if (data.conductorIn[indi].layer == this->getLayer(indj).getGDSIINum())
                {
                    data.conductorIn[indi].zmin = data.stackBegCoor[indj];
                    data.conductorIn[indi].zmax = data.stackEndCoor[indj];
                    data.stackCdtMark[indj] = 1;
                }
            }
        }

        // Use port information in parasitics data member to set fields
        data.numPorts = (this->para).getNPort();
        data.portCoor = (fdtdPort*)malloc(sizeof(fdtdPort) * data.numPorts);
        //data.portArea; // Is this needed at all?
        //data.portEdge; // Is this needed from SolverDataBase?
        //data.portNno; // Is this needed at all?
        for (size_t indi = 0; indi < data.numPorts; indi++)
        {
            Port thisPort = (this->para).getPort(indi); // Get copy of port information for this interation

            // Save coordinates of the port (class automatically ensures x1 < x2, y1 < y2, and z1 < z2)
            data.portCoor[indi].x1 = thisPort.getCoord()[0];
            data.portCoor[indi].y1 = thisPort.getCoord()[1];
            data.portCoor[indi].z1 = thisPort.getCoord()[2];
            data.portCoor[indi].x2 = thisPort.getCoord()[3];
            data.portCoor[indi].y2 = thisPort.getCoord()[4];
            data.portCoor[indi].z2 = thisPort.getCoord()[5];

            // Add coordinates of port to unordered maps if need be
            if (thisPort.getCoord()[0] == thisPort.getCoord()[3])
            {
                // Record this x-coordinate along the x-axis
                portCoorx.insert(thisPort.getCoord()[0]);
            }
            if (thisPort.getCoord()[1] == thisPort.getCoord()[4])
            {
                // Record this y-coordinate along the y-axis
                portCoory.insert(thisPort.getCoord()[1]);
            }
        }

        // Return fdtdMesh
        return data;
    }

    // Print the solver database and dump to SPEF file
    bool printDump()
    {
        // Print
        int numLayer = this->getNumLayer();
        cout << "Solver Database of IC Design, " << this->designName << ":" << endl;
        cout << " Settings for the simulation:" << endl;
        (this->settings).print(); // Print the simulation settings
        cout << " Details of the " << numLayer << " layers:" << endl;
        for (size_t indi = 0; indi < numLayer; indi++)
        {
            (this->layers)[indi].print();
        }
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
            // Write to file
            spef::Spef designPara = (this->para).toSPEF(this->designName, 1.0e-5);
            designPara.dump(spefFile);
            
            // Close file
            spefFile.close();
            return true;
        }
        else
        {
            // File could not be opened
            return false;
        }
    }

    // Destructor
    ~SolverDataBase()
    {
        // Nothing
    }
};

#endif
