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
#include <parser-spef/parser-spef.hpp>
#include <Eigen/Sparse>
#include "limboint.h"
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

// Custom classes for SPEF writer
class SimSettings
{
private:
    double lengthUnit;     // Units for lengths (m)
    vector<double> limits; // xmin, xmax, ymin, ymax, zmin, zmax (m)
    double freqUnit;       // Units for frequency (Hz)
    size_t nFreq;          // Number of frequencies in simulation
    double freqScale;      // Frequency scaling?
    vector<double> freqs;  // List of logarithmically spaced frequencies
public:
    // Default constructor
    SimSettings()
    {
        this->lengthUnit = 1.;
        this->limits = {0., 0., 0., 0., 0., 0.};
        this->freqUnit = 1.;
        this->nFreq = 0;
        this->freqScale = 0.;
        this->freqs = {};
    }

    // Parametrized constructor
    SimSettings(double lengthUnit, vector<double> limits, double freqUnit, size_t nFreq, double freqScale, vector<double> freqs)
    {
        this->lengthUnit = lengthUnit;
        this->limits = limits;
        this->freqUnit = freqUnit;
        this->nFreq = nFreq;
        this->freqScale = freqScale;
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

    // Get number of frequencies in simulatoin
    size_t getNFreq() const
    {
        return this->nFreq;
    }

    // Get frequency scaling
    double getFreqScale() const
    {
        return this->freqScale;
    }

    // Get list of frequencies (logarithmic spacing)
    vector<double> getFreqs() const
    {
        return this->freqs;
    }

    // Print the simulation settings
    void print() const
    {
        size_t numFreq = this->getNFreq();
        cout << " ------" << endl;
        cout << " Simulation Settings:" << endl;
        cout << "  Limits in x-direction: " << (this->limits)[0] << " m to " << (this->limits)[1] << "m" << endl;
        cout << "  Limits in y-direction: " << (this->limits)[2] << " m to " << (this->limits)[3] << "m" << endl;
        cout << "  Limits in z-direction: " << (this->limits)[4] << " m to " << (this->limits)[5] << "m" << endl;
        cout << "  List of " << numFreq << " frequencies to simulate with " << this->freqScale << " scaling:" << endl;
        for (size_t indi = 0; indi < numFreq; indi++)
        {
            if (numFreq - indi == 1) // Single frequency left to print
            {
                cout << "   #" << indi + 1 << " is " << (this->freqs)[indi] << " Hz" << endl;
            }
            else // Two frequencies per line
            {
                cout << "   #" << indi + 1 << " is " << (this->freqs)[indi++] << " Hz, and #" << indi + 2 << " is " << (this->freqs)[indi + 1] << " Hz" << endl;
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
    // Metallic layers are nonnegative, and -1 is used for dielectric, substrates, and undescribed planes
    int getGDSIINum() const
    {
        return this->gdsiiNum;
    }

    // Get layer bottom z-coordinate
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
    size_t nPorts;                             // Number of ports
    vector<std::string> ports;                 // Name of each port
    vector<char> portDir;                      // Direction of each port
    vector<double> Z_port_source;              // Impedance of source attached to port (ohm)
    spMat matG;                                // Conductance matrix (S)
    spMat matC;                                // Capacitance matrix (F)
public:
    // Default constructor
    Parasitics()
    {
        vector<std::string> ports = {};
        vector<char> portDir = {};
        vector<double> Z_port_source = {};
        spMat emptMat;
        this->nPorts = 0;
        this->ports = ports;
        this->portDir = portDir;
        this->matG = emptMat;
        this->matC = emptMat;
    }

    // Parametrized constructor
    Parasitics(size_t nPorts, vector<std::string> ports, vector<char> portDir, vector<double> Z_port_source, spMat matG, spMat matC)
    {
        this->nPorts = nPorts;
        this->ports = ports;
        this->portDir = portDir;
        this->Z_port_source = Z_port_source;
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

    // Get impedance of source attached to port (ohm)
    vector<double> getZPortSource() const
    {
        return this->Z_port_source;
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

    // Print the parasitics information
    void print() const
    {
        int numPort = getNPort();
        cout << " ------" << endl;
        cout << " Parasitics Details:" << endl;
        cout << "  Port Names:" << endl;
        for (size_t indi = 0; indi < numPort; indi++) // Handle each port name
        {
            cout << "   #" << indi + 1 << ": " << (this->ports)[indi] << ", direction " << (this->portDir)[indi] << endl;
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
            para.name_map.emplace(indi + 1, (this->ports)[indi]); // Create name map for each port
            para.ports.emplace_back("*" + to_string(indi + 1)); // Instantiate and push new port entry by name map
            switch ((this->portDir)[indi]) // Assign port direction
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
        for (size_t indi = 0; indi < numPort; indi++)
        {
            para.nets.back().connections.emplace_back(spef::Connection());
            para.nets.back().connections.back().name = (this->ports)[indi];
            para.nets.back().connections.back().type = spef::ConnectionType::EXTERNAL; // All ports are external connections, not internal to cells
            para.nets.back().connections.back().direction = (para.ports[indi]).direction; // Same as port direction
            for (spMat::InnerIterator it(this->matC, indi); it; ++it)
            {
                if ((it.row() == it.col()) && (abs(this->getCNodeGround(indi)) >= saveThresh * capTot)) // Diagonal element, so sufficient admittance to ground
                {
                    para.nets.back().caps.emplace_back(forward_as_tuple((para.ports)[indi].name, "", this->getCNodeGround(indi)));
                }
                else if (it.col() > it.row()) // Upper triangular element, so negate node-to-node admittance
                {
                    para.nets.back().caps.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[it.col()].name, -it.value()));
                } // Skip lower triangular element, which is repeated
            }
            for (spMat::InnerIterator it(this->matG, indi); it; ++it)
            {
                if ((it.row() == it.col()) && (abs(this->getGNodeGround(indi)) >= saveThresh * condTot)) // Diagonal element, so sufficient admittance to ground
                {
                    para.nets.back().ress.emplace_back(forward_as_tuple((para.ports)[indi].name, "", 1.0 / this->getGNodeGround(indi)));
                }
                else if (it.col() > it.row()) // Upper triangular element, so negate node-to-node admittance
                {
                    para.nets.back().ress.emplace_back(forward_as_tuple((para.ports)[indi].name, (para.ports)[it.col()].name, -1.0 / it.value()));
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
        /*vector<vector<double>> emptMat = { {} };*/
        spMat emptMat;
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

    // Set SPEF output file name
    void setOutSPEF(std::string fileName)
    {
        this->outSPEF = fileName;
    }

    // Find index of layer by name
    // Returns index past number of layers if not found
    size_t locateLayer(std::string name) const
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

    // Read IMAP 3D file data into the solver database
    bool readIMAPwriteGDSII(std::string imapFileName, std::string gdsiiFileName)
    {
        // Attempt to open IMAP 3D file
        ifstream imapFile(imapFileName.c_str());
        if (imapFile.is_open())
        {
            // Time manipulations
            char timeStr[80]; // Character array to hold formatted time
            time_t rawtime; // Timer variable
            strftime(timeStr, sizeof(timeStr), "%d-%m-%Y %H:%M:%S", localtime(&rawtime)); // Use local timezone to format string
            std::string time(timeStr); // Formatted time to string parametrized constructor

            // Build single geometric cell from limboint.h to store information
            GeoCell cellIMAP;
            cellIMAP.cellName = imapFileName.substr(0, imapFileName.find(".", 1));
            cellIMAP.dateCreate = time;
            cellIMAP.dateMod = time;

            // Build ASCII database from limboint.h
            AsciiDataBase adbIMAP;
            adbIMAP.setFileName(gdsiiFileName);
            adbIMAP.setDateMod(time);
            adbIMAP.setDateAccess(time);

            // File is readable line-by-line
            std::string fileLine;
            getline(imapFile, fileLine);

            // Save opening lines as metadata (mostly ignored)
            if (fileLine.compare(0, 7, "IMAP3D ") == 0)
            {
                // Save IMAP 3D version number
                std::string version = fileLine.substr(7, fileLine.length() - 7); // 7 characters in header syntax
            }

            // Read rest of file line-by-line
            while (!imapFile.eof())
            {
                // Handle units
                double stripLength = 2000; // Hard-coded (for now)
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
                    adbIMAP.setdbUserUnits(1.);
                    adbIMAP.setdbUnits(multSI);
                }
                // Handle design name
                else if (fileLine.compare(0, 5, "NAME ") == 0)
                {
                    this->designName = fileLine.substr(5, fileLine.length() - 6);
                }
                // Handle layer stack
                if (fileLine.compare(0, 5, "STACK") == 0)
                {
                    // Move down one line
                    getline(imapFile, fileLine);

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
                            }
                            double zStart = 0.;
                            if (indZStart != string::npos)
                            {
                                zStart = stod(fileLine.substr(indZStart + 2, indHeight - indZStart - 3)) * adbIMAP.getdbUnits();
                            }
                            double zHeight = stod(fileLine.substr(indHeight + 2, fileLine.find(" ", indHeight) - indHeight - 2)) * adbIMAP.getdbUnits();
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
                        getline(imapFile, fileLine);
                    }
                }
                // Handle conductors
                if (fileLine.compare(0, 10, "CONDUCTORS") == 0)
                {
                    // Move down one line
                    getline(imapFile, fileLine);
                    stripLength *= adbIMAP.getdbUnits(); // Hard-coded rescaling of length (for now)

                    // Keep reading until end of conductor list
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
                            double x1 = stod(fileLine.substr(indX1 + 3, indY1 - indX1 - 4)) * adbIMAP.getdbUnits(); // Length is index difference minus space
                            double y1 = stod(fileLine.substr(indY1 + 3, indZ1 - indY1 - 4)) * adbIMAP.getdbUnits();
                            double x2 = stod(fileLine.substr(indX2 + 3, indY2 - indX2 - 4)) * adbIMAP.getdbUnits();
                            double y2 = stod(fileLine.substr(indY2 + 3, indZ2 - indY2 - 4)) * adbIMAP.getdbUnits();
                            std::string sigma = fileLine.substr(indSigma, indLayer - indSigma - 1);
                            std::string group;
                            if (indGroup != string::npos)
                            {
                                group = fileLine.substr(indGroup + 6);
                            }

                            // Assign layer number based on layer stack-up
                            int gdsiiNum = 1;
                            size_t indThisLayer = this->locateLayer(fileLine.substr(indLayer + 6, fileLine.find(" ", indLayer)));
                            if (indThisLayer < (this->layers).size())
                            {
                                // Layer name for conductor matched name in layer stack-up
                                gdsiiNum = ((this->layers)[indThisLayer]).getGDSIINum();
                            }

                            // Push new box to the geometric cell
                            cellIMAP.boxes.emplace_back(box({ x2, y1, x2, y2 + stripLength, x1, y2 + stripLength, x1, y1, x2, y1 }, gdsiiNum, {sigma, condName, category, group}, 0));
                        }
                        // Keep moving down the conductor list
                        getline(imapFile, fileLine);
                    }
                }
                // Handle port table
                if (fileLine.compare(0, 9, "PORTTABLE") == 0)
                {
                    // Move down one line
                    getline(imapFile, fileLine);

                    // Keep reading until end of port table
                    vector<std::string> ports;
                    vector<char> portDir;
                    vector<double> Z_source;
                    while (fileLine.compare(0, 8, "ANALYSIS") != 0)
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
                            ports.push_back(groupName);
                            portDir.push_back('B');
                            Z_source.push_back(Z_near);
                        }
                        // Keep moving down the port table
                        getline(imapFile, fileLine);
                    }

                    // Propagate port information to Solver Database now
                    this->para = Parasitics(ports.size(), ports, portDir, Z_source, spMat(), spMat());
                }
                // Handle analysis parameters
                if (fileLine.compare(0, 8, "ANALYSIS") == 0)
                {
                    // Move down one line
                    getline(imapFile, fileLine);
                }
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
                    int numFreqPts = stoi(fileLine.substr(indNoP + 15));
                    vector<double> freqList;
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

                    // Save simulation settings (hard-coded limits for now)
                    this->settings = SimSettings(adbIMAP.getdbUnits(), {0.0, 3.000e-4, 0.0, 2.00005e-3, (this->layers).back().getZStart(), (this->layers).front().getZStart() + (this->layers).front().getZHeight()}, 1.0, (size_t) numFreqPts, 0.0, freqList);
                }
                // Handle length in this weird spot
                if (fileLine.compare(0, 6, "Length") == 0)
                {
                    //double stripLength = stod(fileLine.substr(7));
                }

                // Keep reading new lines in file
                getline(imapFile, fileLine);
            }

            // Close file
            imapFile.close();

            // Update ASCII database
            adbIMAP.setLibName(this->designName);
            adbIMAP.appendCell(cellIMAP);
            adbIMAP.setdbUnits(adbIMAP.getdbUnits() * 1.e-3); // Rescale IMAP 0.001x to allow integer representation in GDSII

            // Print the ASCII database
            adbIMAP.print({ 0 });

            // Write GDSII file to hard drive
            bool dumpPassed = adbIMAP.dump();
            return dumpPassed;
        }
        else
        {
            // File could not be opened
            return false;
        }
    }

    // Print the solver database and dump to SPEF file
    bool printDump()
    {
        // Print
        int numLayer = this->getNumLayer();
        cout << "Solver Database of IC Design, " << this->designName << ":" << endl;
        cout << " Settings for the simulation:" << endl;
        (this->settings).print(); // Print the simulation settings
        cout << " Layers:" << endl;
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