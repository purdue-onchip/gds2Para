// Make sure this file is included only once
#ifndef Interface_h
#define Interface_h

#include <iostream>
#include <fstream>
#include <vector>

// Definition of Geometry class to store GDT file information
// Must first use GDS2GDT in shell
using namespace std;

class GeoCell
{
  private:
    std::string cellName;                    // Cell name
	std::string metadata;                    // Included cell metadata
	vector<unsigned long> boundLayers;       // Boundary layer number vector
	vector<vector<double>*> bounds;          // Boundaries in cell coordinates vector
	vector<vector<std::string>*> boundProps; // Boundary properties vector
	vector<unsigned long> pathLayers;        // Path layer number vector
	vector<unsigned long> pathTypes;         // Path type vector
	vector<double> pathWidths;               // Path width vector
	vector<vector<double>*> paths;           // Paths in cell coordinates vector
	vector<vector<std::string>*> pathProps;  // Path properties vector
	vector<std::string> srefNames;           // Structure reference name vector
	vector<vector<double>> srefs;            // Structure references in cell coordinate pairs
	vector<vector<std::string>*> srefProps;  // Structure reference properties vector
	vector<unsigned long> nodeLayers;        // Node layer number vector
	vector<unsigned long> nodeTypes;         // Node type vector
	vector<vector<double>*> nodes;           // Nodes in cell coordinates vector
	vector<vector<std::string>*> nodeProps;  // Node properties vector
	vector<unsigned long> textLayers;        // Text layer number vector
	vector<double> textWidths;               // Text width vector
	vector<vector<unsigned long>> textJusts; // Text justification vector
	vector<vector<double>> texts;            // Text boxes in cell coordinate pairs
	vector<std::string> textStr;             // Text in text boxes
  public:
    // Constructor
	GeoCell(std::string, std::string, vector<unsigned long>, vector<vector<double>*>, vector<vector<std::string>*>, vector<unsigned long>, vector<unsigned long>, vector<double>, vector<vector<double>*>, vector<vector<std::string>*>, vector<std::string>, vector<vector<double>>, vector<vector<std::string>*>, vector<unsigned long>, vector<unsigned long>, vector<vector<double>*>, vector<vector<std::string>*>, vector<unsigned long>, vector<double>, vector<vector<unsigned long>>, vector<vector<double>>, vector<std::string>);
	
	// Modify a geometric cell
	bool modCell();

    // Get the name of the cell
    std::string getCellName() const;

    // Get the metadata of the cell
    std::string getCellMeta() const;

    // Get number of boundaries
    int getNumBound() const;

    // Get number of paths
    int getNumPath() const;

    // Get number of structure references
    int getNumSRef() const;

	// Get number of nodes
	int getNumNode() const;

    // Get number of text boxes
    int getNumText() const;

    // Return boundary layer numbers
    vector<unsigned long> getBoundLayers() const;

    // Return boundaries in cell coordinates
    vector<vector<double>*> getBounds() const;

    // Return boundary properties
    vector<vector<std::string>*> getBoundProps() const;

    // Return path layer numbers
    vector<unsigned long> getPathLayers() const;

    // Return path types
	// 0 = square ends at vertices, 1 = round ends, 2 = square ends overshoot vertices by half width
    vector<unsigned long> getPathTypes() const;

    // Return path widths
	// Negative values mean independent of structure scaling
    vector<double> getPathWidths() const;

    // Return paths in cell coordinates
    vector<vector<double>*> getPaths() const;

    // Return path properties
    vector<vector<std::string>*> getPathProps() const;

    // Return structure reference names
    vector<std::string> getSRefNames() const;

    // Return coordinate pairs of placed structure references
    vector<vector<double>> getSRefs() const;

    // Return structure reference properties
    vector<vector<std::string>*> getSRefProps() const;

	// Return node layer numbers
	vector<unsigned long> getNodeLayers() const;

	// Return node types
	vector<unsigned long> getNodeTypes() const;

	// Return nodes in cell coordinates
	vector<vector<double>*> getNodes() const;

	// Return node properties
	vector<vector<std::string>*> getNodeProps() const;

    // Return text layer numbers
    vector<unsigned long> getTextLayers() const;

    // Return text widths
	// Negative values mean independent of structure scaling
    vector<double> getTextWidths() const;

    // Return text justifications
    vector<vector<unsigned long>> getTextJusts() const;

    // Return coordinate pairs of placed text boxes
    vector<vector<double>> getTexts() const;

    // Return text in text boxes
    vector<std::string> getTextStr() const;

    // Print the geometric cell
	void print() const;

	// Destructor
	~GeoCell();
};

class Geometry
{
  private:
    std::string fileName;   // File name
    std::string metadata;   // Included metadata
    long fileSize;          // Number of bytes in file besides metadata
	vector<GeoCell> cells;  // Vector of cells in design
  public:
    // Constructor
    Geometry();

    // Load file contents
    // Return false if unable to load.
    bool load(std::string filename);

    // Get file size
    long getSize() const;

	// Get number of cells
	int getNumCell() const;

    // Find index of cell by name
    size_t locateCell(std::string name) const;

    // Return a cell
    GeoCell getCell(size_t indCell) const;

    // Return all cell names
    vector<string> findNames() const;

    // Return all included layers in order
    vector<unsigned long> findLayers() const;

	// Return via list
	vector<vector<double>> findVias(size_t indCell, vector<double> center, std::string viaName);

    // Print the geometry
    void print() const;

    // Write gnuplot file for layer of a cell
    vector<int> writeCellPlotter(ofstream& gnuFile, size_t indCell, vector<double> center, unsigned long layer, int bNum, int pNum, int nNum, int tNum);

	// Write gnuplot environment
	// Return false if unable to write a gnuplot file
	bool writePlotEnv(std::string filePrefix, std::string masterCell, unsigned long pixels, vector<double> layerBounds);

    // Destructor
    ~Geometry();
};

#endif
