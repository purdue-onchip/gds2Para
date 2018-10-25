/**
* @file   limboint.h
* @author Michael R. Hayashi
* @date   18 October 2018
* @brief  Header for Limbo/Geometry Handling Interface for GDSII Files
*/

// Make sure this file is included only once
#ifndef limboint_h
#define limboint_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
using std::cout;
using std::endl;

class boundary
{
private:
	vector<double>* bounds;        // Coordinates of boundary
	int layer;                     // Layer of boundary
	vector<std::string>* props;    // Properties of boundary
public:
	// Default constructor
	boundary()
	{
		this->bounds = new vector<double>;
		this->layer = 0;
		this->props = NULL;
	}

	// Parametrized constructor
	boundary(vector<double>* bounds, int layer, vector<std::string>* props)
	{
		this->bounds = bounds;
		this->layer = layer;
		this->props = props;
	}

	// Get boundary points in cell coordinates
	vector<double>* getBounds() const
	{
		return this->bounds;
	}

	// Get boundary layer number
	int getLayer() const
	{
		return this->layer;
	}

	// Get boundary properties
	vector<std::string>* getProps() const
	{
		return this->props;
	}

	// Destructor
	~boundary()
	{
		//this->bounds = NULL;
		delete bounds;
		//(this->bounds).clear();
		this->layer = 0;
		//this->props = NULL;
		delete props;
	}
};

class path
{
private:
	vector<double>* paths;         // Coordinates of path
	int layer;                     // Layer of path
	vector<std::string>* props;    // Properties of path
	int type;                      // Endcap type of path
	double width;                  // Width of path (m)
public:
	// Default constructor
	path()
	{
		this->paths = new vector<double>;
		this->layer = 0;
		this->props = NULL;
		this->type = 2;
		this->width = 0.;
	}

	// Parametrized constructor
	path(vector<double>* paths, int layer, vector<std::string>* props, int type, double width)
	{
		this->paths = paths;
		this->layer = layer;
		this->props = props;
		if ((type >= 0) || (type <= 2))
		{
			this->type = type;
		}
		else
		{
			cout << "Path type must be 0, 1, or 2. Defaulting to type 2." << endl;
			this->type = 2; // Default to type 2 if invalid
		}
		this->width = width;
	}

	// Get path vertices in cell coordinates
	vector<double>* getPaths() const
	{
		return this->paths;
	}

	// Get path layer number
	int getLayer() const
	{
		return this->layer;
	}

	// Get path properties
	vector<std::string>* getProps() const
	{
		return this->props;
	}

	// Get path type
	// 0 = square ends at vertices, 1 = round ends, 2 = square ends overshoot vertices by half width
	int getType() const
	{
		return this->type;
	}

	// Get path width
	// Negative values mean independent of structure scaling
	double getWidth() const
	{
		return this->width;
	}

	// Destructor
	~path()
	{
		//this->paths = NULL;
		delete paths;
		this->layer = 0;
		delete props;
		this->type = 2;
		this->width = 0.;
	}
};

class node
{
private:
	vector<double>* nodes;         // Coordinates of electrical node
	int layer;                     // Layer of electrical node
	vector<std::string>* props;    // Properties of electrical node
	int type;                      // Type of electrical node
public:
	// Default constructor
	node()
	{
		this->nodes = new vector<double>;
		this->layer = 0;
		this->props = NULL;
		this->type = 0;
	}

	// Parametrized constructor
	node(vector<double>* nodes, int layer, vector<std::string>* props, int type)
	{
		this->nodes = nodes;
		this->layer = layer;
		this->props = props;
		this->type = 0;
	}

	// Get node vertices in cell coordinates
	vector<double>* getNodes() const
	{
		return this->nodes;
	}

	// Get node layer number
	int getLayer() const
	{
		return this->layer;
	}

	// Get node properties
	vector<std::string>* getProps() const
	{
		return this->props;
	}

	// Get node type
	int getType() const
	{
		return this->type;
	}

	// Destructor
	~node()
	{
		//this->nodes = NULL;
		delete nodes;
		this->layer = 0;
		delete props;
		this->type = 0;
	}
};

class box
{
private:
	vector<double>* boxes;         // Coordinates of box outline
	int layer;                     // Layer of box outline
	vector<std::string>* props;    // Properties of box outline
	int type;                      // Type of box outline
public:
	// Default constructor
	box()
	{
		this->boxes = new vector<double>;
		this->layer = 0;
		this->props = NULL;
		this->type = 0;
	}

	// Parametrized constructor
	box(vector<double>* boxes, int layer, vector<std::string>* props, int type)
	{
		this->boxes = boxes;
		this->layer = layer;
		this->props = props;
		this->type = 0;
	}

	// Get box outline vertices in cell coordinates
	vector<double>* getBoxes() const
	{
		return this->boxes;
	}

	// Get box outline layer number
	int getLayer() const
	{
		return this->layer;
	}

	// Get box outline properties
	vector<std::string>* getProps() const
	{
		return this->props;
	}

	// Get box outline type
	int getType() const
	{
		return this->type;
	}

	// Destructor
	~box()
	{
		//this->boxes = NULL;
		delete boxes;
		this->layer = 0;
		delete props;
		this->type = 0;
	}
};

class textbox
{
private:
	vector<double> texts;          // Coordinate pair of text box center
	int layer;                     // Layer of text box
	vector<std::string>* props;    // Properties of text box
	int type;                      // Type of text box
	int fontID;                    // Font identifier for text box
	vector<int> justs;             // Justifications of text box
	double width;                  // Width of text (m)
	std::string textStr;           // String in text box
public:
	// Default constructor
	textbox()
	{
		vector<double> texts = { 0., 0. };
		vector<int> justs = { 0, 0 };
		this->texts = texts;
		this->layer = 0;
		this->props = NULL;
		this->type = 0;
		this->fontID = 0;
		this->justs = justs;
		this->width = 0.;
		this->textStr = "";
	}

	// Parametrized constructor
	textbox(vector<double> texts, int layer, vector<std::string>* props, int type, int fontID, vector<int> justs, double width, std::string textStr)
	{
		this->texts = texts;
		this->layer = layer;
		this->props = props;
		this->type = type;
		this->fontID = fontID;
		this->justs = justs;
		this->width = width;
		this->textStr = textStr;
	}

	// Get text box center in cell coordinates
	vector<double> getTexts() const
	{
		return this->texts;
	}

	// Get text box layer number
	int getLayer() const
	{
		return this->layer;
	}

	// Get text box properties
	vector<std::string>* getProps() const
	{
		return this->props;
	}

	// Get text box type
	int getType() const
	{
		return this->type;
	}

	// Get text box font ID (0-3)
	int getFontID() const
	{
		return this->fontID;
	}

	// Get text box justifications
	// First (vertical): 0 = top, 1 = middle, 2 = bottom
	// Second (horizontal): 0 = left, 1 = center, 2 = right
	vector<int> getJusts() const
	{
		return this->justs;
	}

	// Get text width
	// Negative values mean independent of structure scaling
	double getWidth() const
	{
		return this->width;
	}

	// Get string in text box
	std::string getTextStr() const
	{
		return this->textStr;
	}

	// Destructor
	~textbox()
	{
		vector<double> texts = { 0., 0. };
		vector<int> justs = { 0, 0 };
		this->texts = texts;
		this->layer = 0;
		delete props;
		this->type = 2;
		this->fontID = 0;
		this->justs = justs;
		this->width = 0.;
		this->textStr = "";
	}
};

class sref
{
private:
	vector<double> srefs;          // Coordinate pair of placed structure reference
	vector<std::string>* props;    // Properties of structure reference
	std::string srefName;          // Name of structure reference
public:
	// Default constructor
	sref()
	{
		vector<double> srefs = { 0., 0. };
		this->srefs = srefs;
		this->props = NULL;
		this->srefName = "";
	}

	// Parametrized constructor
	sref(vector<double> srefs, vector<std::string>* props, std::string srefName)
	{
		this->srefs = srefs;
		this->props = props;
		this->srefName = srefName;
	}

	// Get placed structure reference center in cell coordinates
	vector<double> getSRefs() const
	{
		return this->srefs;
	}

	// Get structure reference properties
	vector<std::string>* getProps() const
	{
		return this->props;
	}

	// Get string in text box
	std::string getSRefName() const
	{
		return this->srefName;
	}

	// Destructor
	~sref()
	{
		vector<double> srefs = { 0., 0. };
		this->srefs = srefs;
		delete props;
		this->srefName = "";
	}
};

class GeoCell
{
public:
	std::string cellName;                        // Cell name
	std::string dateCreate;                      // Date of cell creation
	std::string dateMod;                         // Date of last modification
	vector<boundary> boundaries;                 // All boundaries in geometric cell
	vector<path> paths;                          // All paths in geometric cell
	vector<node> nodes;                          // All nodes in geometric cell
	vector<box> boxes;                           // All box outlines in geometric cell
	vector<textbox> textboxes;                   // All text boxes in geometric cell
	vector<sref> sreferences;                    // All structure references in geometric cell
public:
	// Default constructor
	GeoCell()
	{
		std::string cellName, dateCreate, dateMod;
		vector<boundary> boundaries;
		vector<path> paths;
		vector<node> nodes;
		vector<box> boxes;
		vector<textbox> textboxes;
		vector<sref> sreferences;
		this->cellName = cellName;
		this->dateCreate = dateCreate;
		this->dateMod = dateMod;
		this->boundaries = boundaries;
		this->paths = paths;
		this->nodes = nodes;
		this->boxes = boxes;
		this->textboxes = textboxes;
		this->sreferences = sreferences;
	}

	// Get the name of the cell
	std::string getCellName() const
	{
		return this->cellName;
	}

	// Get the cell creation date
	std::string getDateCreate() const
	{
		return this->dateCreate;
	}

	// Get the last modification date of cell
	std::string getDateMod() const
	{
		return this->dateMod;
	}

	// Get vector of boundaries
	vector<boundary> getBoundVec() const
	{
		return this->boundaries;
	}

	// Get vector of paths
	vector<path> getPathVec() const
	{
		return this->paths;
	}

	// Get vector of nodes
	vector<node> getNodeVec() const
	{
		return this->nodes;
	}

	// Get vector of box outlines
	vector<box> getBoxVec() const
	{
		return this->boxes;
	}

	// Get vector of text boxes
	vector<textbox> getTextVec() const
	{
		return this->textboxes;
	}

	// Get vector of structure references
	vector<sref> getSRefVec() const
	{
		return this->sreferences;
	}

	// Return number of boundaries
	int getNumBound() const
	{
		return (this->boundaries).size();
	}

	// Return number of paths
	int getNumPath() const
	{
		return (this->paths).size();
	}

	// Return number of nodes
	int getNumNode() const
	{
		return (this->nodes).size();
	}

	// Return number of box outlines
	int getNumBox() const
	{
		return (this->boxes).size();
	}

	// Return number of text boxes
	int getNumText() const
	{
		return (this->textboxes).size();
	}

	// Return number of structure references
	int getNumSRef() const
	{
		return (this->sreferences).size();
	}

	// Print the geometric cell
	void print() const
	{
		int numBound = GeoCell::getNumBound();
		int numPath = GeoCell::getNumPath();
		int numNode = GeoCell::getNumNode();
		int numBox = GeoCell::getNumBox();
		int numText = GeoCell::getNumText();
		int numSRef = GeoCell::getNumSRef();

		cout << " ------" << endl;
		cout << " Cell Details:" << endl;
		cout << "  Cell Name: " << this->cellName << endl;
		cout << "  Created: " << this->dateCreate << endl;
		cout << "  Last Modified: " << this->dateMod << endl;
		cout << "  List of " << numBound << " boundaries:" << endl;
		for (size_t indi = 0; indi < numBound; indi++) // Handle each boundary
		{
			if (((this->boundaries)[indi]).getProps() != NULL) // Check if NULL pointer
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->boundaries)[indi]).getLayer() << " with " << (*(((this->boundaries)[indi]).getProps()))[0] << " as Property " << 1 << endl;
			}
			else
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->boundaries)[indi]).getLayer() << endl;
			}
			char point[128];
			string strPoints;
			vector<double> boundCoord = *(((this->boundaries)[indi]).getBounds());
			for (size_t indj = 0; indj < boundCoord.size(); indj++) // C string for each ordered pair
			{
				sprintf(point, "(%1.4g, %1.4g) ", boundCoord[indj++], boundCoord[indj + 1]);
				strPoints.append(point);
			}
			cout << "    " << strPoints << endl; // Print line of ordered pairs
			//delete[] point; // Free array pointer
		}
		cout << "  List of " << numPath << " paths:" << endl;
		for (size_t indi = 0; indi < numPath; indi++) // Handle each path
		{
			if (((this->paths)[indi]).getProps() != NULL) // Check if NULL pointer
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->paths)[indi]).getLayer() << " with path type #" << ((this->paths)[indi]).getType() << " and width " << ((this->paths)[indi]).getWidth() << " with " << (*(((this->paths)[indi]).getProps()))[0] << " as Property " << 1 << endl;
			}
			else
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->paths)[indi]).getLayer() << " with path type #" << ((this->paths)[indi]).getType() << " and width " << ((this->paths)[indi]).getWidth() << endl;
			}
			char point[128];
			string strPoints;
			vector<double> pathCoord = *(((this->paths)[indi]).getPaths());
			for (size_t indj = 0; indj < pathCoord.size(); indj++) // C string for each ordered pair
			{
				sprintf(point, "(%1.4g, %1.4g) ", pathCoord[indj++], pathCoord[indj + 1]);
				strPoints.append(point);
			}
			cout << "    " << strPoints << endl; // Print line of ordered pairs
			//delete[] point; // Free array pointer
		}
		cout << "  List of " << numNode << " nodes:" << endl;
		for (size_t indi = 0; indi < numNode; indi++) // Handle each node
		{
			if (((this->nodes)[indi]).getProps() != NULL) // Check if NULL pointer
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->nodes)[indi]).getLayer() << " with node type #" << ((this->nodes)[indi]).getType() << " with " << (*(((this->nodes)[indi]).getProps()))[0] << " as Property " << 1 << endl;
			}
			else
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->nodes)[indi]).getLayer() << " with node type #" << ((this->nodes)[indi]).getType() << endl;
			}
			char point[128];
			string strPoints;
			vector<double> nodeCoord = *(((this->nodes)[indi]).getNodes());
			for (size_t indj = 0; indj < nodeCoord.size(); indj++) // C string for each ordered pair
			{
				sprintf(point, "(%1.4g, %1.4g) ", nodeCoord[indj++], nodeCoord[indj + 1]);
				strPoints.append(point);
			}
			cout << "    " << strPoints << endl; // Print line of ordered pairs
			//delete[] point; // Free array pointer
		}
		cout << "  List of " << numBox << " box outlines:" << endl;
		for (size_t indi = 0; indi < numBox; indi++) // Handle each box outline
		{
			if (((this->boxes)[indi]).getProps() != NULL) // Check if NULL pointer
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->boxes)[indi]).getLayer() << " with box outline type #" << ((this->boxes)[indi]).getType() << " with " << (*(((this->boxes)[indi]).getProps()))[0] << " as Property " << 1 << endl;
			}
			else
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->boxes)[indi]).getLayer() << " with box outline type #" << ((this->boxes)[indi]).getType() << endl;
			}
			char point[128];
			string strPoints;
			vector<double> boxCoord = *(((this->boxes)[indi]).getBoxes());
			for (size_t indj = 0; indj < boxCoord.size(); indj++) // C string for each ordered pair (should be 5, sometimes 4)
			{
				sprintf(point, "(%1.4g, %1.4g) ", boxCoord[indj++], boxCoord[indj + 1]);
				strPoints.append(point);
			}
			cout << "    " << strPoints << endl; // Print line of ordered pairs
			//delete[] point; // Free array pointer
		}
		cout << "  List of " << numText << " text boxes:" << endl;
		for (size_t indi = 0; indi < numText; indi++) // Handle each text box
		{
			if (((this->textboxes)[indi]).getProps() != NULL)  // Check if NULL pointer
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->textboxes)[indi]).getLayer() << " with text box type #" << ((this->textboxes)[indi]).getType() << " and message \"" << ((this->textboxes)[indi]).getTextStr() << "\" with " << (*(((this->textboxes)[indi]).getProps()))[0] << " as Property " << 1 << endl;
			}
			else
			{
				cout << "   #" << indi + 1 << " on Layer " << ((this->textboxes)[indi]).getLayer() << " with text box type #" << ((this->textboxes)[indi]).getType() << " and message \"" << ((this->textboxes)[indi]).getTextStr() << "\"" << endl;
			}
			char point[128];
			sprintf(point, "(%1.4g, %1.4g) ", (((this->textboxes)[indi]).getTexts())[0], (((this->textboxes)[indi]).getTexts())[1]);
			cout << "    " << point << endl; // Print ordered pair of center
			//delete[] point; // Free array pointer
		}
		cout << "  List of " << numSRef << " structure references:" << endl;
		for (size_t indi = 0; indi < numSRef; indi++) // Handle each structure reference
		{
			if (((this->sreferences)[indi]).getProps() != NULL) // Check if NULL pointer
			{
				cout << "   #" << indi + 1 << " named " << ((this->sreferences)[indi]).getSRefName() << " with " << (*(((this->sreferences)[indi]).getProps()))[0] << " as Property " << 1 << endl;
			}
			else
			{
				cout << "   #" << indi + 1 << " named " << ((this->sreferences)[indi]).getSRefName() << endl;
			}
			char point[128];
			sprintf(point, "(%1.4g, %1.4g)", (((this->sreferences)[indi]).getSRefs())[0], (((this->sreferences)[indi]).getSRefs())[1]);
			cout << "    " << point << endl; // Print ordered pair of center
			//delete[] point; // Free array pointer
		}
		cout << " ------" << endl;
	}

	// Destructor
	~GeoCell()
	{
		this->cellName = "";
		this->dateCreate = "";
		this->dateMod = "";
		/*~(this->boundaries);
		~(this->paths);
		~(this->nodes);
		~(this->boxes);
		~(this->textboxes);
		~(this->sreferences);*/
	}
};

/// @brief test ascii callbacks 
struct AsciiDataBase : public GdsParser::GdsDataBase
{
private:
	std::string fileName;                        // File name
	std::string version;                         // GDSII version number
	std::string dateMod;                         // Date of last modification
	std::string dateAccess;                      // Date of last access
	std::string libName;                         // Name of this library
	double dbUserUnits;                          // Database units measured in user units
	double dbUnits;                              // Database units measured in SI (m)
	size_t numCell;                              // Index of present working geometric cell
	char element;                                // Present working element read
	size_t numProp;                              // Present working property number
	vector<GeoCell> cells;                       // Vector of cells in design
public:
	/// @brief constructor
	AsciiDataBase()
	{
		std::string fileName, version, dateMod, dateAccess, libName;
		double databaseUserUnits = 1.;
		double databaseUnits = 1.;
		size_t numCell = 0;
		char element;
		size_t numProp = 0;
		vector<GeoCell> cells;

		this->fileName = fileName;
		this->version = version;
		this->dateMod = dateMod;
		this->dateAccess = dateAccess;
		this->libName = libName;
		this->dbUserUnits = databaseUserUnits;
		this->dbUnits = databaseUnits;
		this->numCell = numCell;
		this->element = element;
		this->numProp = 0;
		this->cells = cells;
	}

	// Get file name
	std::string getFileName() const
	{
		return this->fileName;
	}

	// Get version
	std::string getVersion() const
	{
		return this->version;
	}

	// Get date of last modification
	std::string getDateMod() const
	{
		return this->dateMod;
	}

	// Get date of last access
	std::string getDateAccess() const
	{
		return this->dateAccess;
	}

	// Get library name
	std::string getLibName() const
	{
		return this->libName;
	}

	// Get database user units
	double getdbUserUnits() const
	{
		return this->dbUserUnits;
	}

	// Get database units in SI
	double getdbUnits() const
	{
		return this->dbUnits;
	}

	// Get total number of cells (current cell index when parsing)
	size_t getNumCell() const
	{
		return this->numCell;
	}

	// Get current element designator when parsing
	char getElement() const
	{
		return this->element;
	}

	// Get current property number
	size_t getNumProp()
	{
		return this->numProp;
	}

	// Set file name
	void setFileName(std::string fileName)
	{
		this->fileName = fileName;
	}

	// Find index of cell by name
	// Returns index past number of cells if not found
	size_t locateCell(std::string name) const
	{
		size_t indCell;
		for (indCell = 0; indCell < this->getNumCell(); indCell++)
		{
			if (name.compare((this->cells[indCell]).getCellName()) == 0)
			{
				return indCell;
			}
		}
		return indCell;
	}

	// Return a cell
	GeoCell getCell(size_t indCell) const
	{
		return (this->cells)[indCell];
	}

	// Return all cell names
	vector<string> findNames() const
	{
		vector<string> names;
		for (size_t indi = 0; indi < this->getNumCell(); indi++)
		{
			GeoCell cell = this->getCell(indi);
			names.push_back(cell.getCellName());
		}
		return names;
	}

	// Return all included layers in order
	vector<int> findLayers() const
	{
		vector<int> layers;
		for (size_t indi = 0; indi < this->getNumCell(); indi++)
		{
			GeoCell cell = this->getCell(indi);
			for (size_t indj = 0; indj < cell.getNumBound(); indj++)
			{
				layers.push_back(((cell.boundaries)[indj]).getLayer());
			}
			for (size_t indj = 0; indj < cell.getNumPath(); indj++)
			{
				layers.push_back(((cell.paths)[indj]).getLayer());
			}
			for (size_t indj = 0; indj < cell.getNumNode(); indj++)
			{
				layers.push_back(((cell.nodes)[indj]).getLayer());
			}
			for (size_t indj = 0; indj < cell.getNumBox(); indj++)
			{
				layers.push_back(((cell.boxes)[indj]).getLayer());
			}
			for (size_t indj = 0; indj < cell.getNumText(); indj++)
			{
				layers.push_back(((cell.textboxes)[indj]).getLayer());
			}
		}
		sort(layers.begin(), layers.end());
		layers.erase(unique(layers.begin(), layers.end()), layers.end());
		return layers;
	}

	// Return via list
	vector<vector<double>> findVias(size_t indCell, vector<double> center, std::string viaName)
	{
		GeoCell thisCell = this->getCell(indCell);
		size_t indVia = this->locateCell(viaName);
		vector<vector<double>> viaList;
		for (size_t indi = 0; indi < thisCell.getNumSRef(); indi++)
		{
			string srefName = (thisCell.sreferences)[indi].getSRefName();
			size_t indNextCell = this->locateCell(srefName);
			if (indNextCell == indVia) // Cell is a via ...
			{
				viaList.push_back((thisCell.sreferences)[indi].getSRefs());
			}
			else if (indNextCell < this->getNumCell()) // Non-via cell name was found ...
			{
				vector<double> thisSRef = (thisCell.sreferences)[indi].getSRefs();
				vector<vector<double>> newVias = this->findVias(indNextCell, { (center[0] + thisSRef[0]), (center[1] + thisSRef[1]) }, viaName); // Recursion step
				viaList.insert(viaList.end(), newVias.begin(), newVias.end());
			}
		}
		return viaList;
	}

	// Print the ASCII database with the design geometry
	void print(vector<size_t> indCellPrint) const
	{
		int numCell = getNumCell();

		cout << "ASCII Database of IC Design:" << endl;
		cout << " File Name: " << this->fileName << endl;
		cout << " Version: " << this->getVersion() << endl;
		cout << " Date of last modification: " << this->getDateMod() << endl;
		cout << " Date of last access: " << this->getDateAccess() << endl;
		cout << " Name of this library: " << this->getLibName() << endl;
		cout << " Database units: " << this->getdbUnits() << " m" << endl;
		cout << " Database units: " << this->getdbUserUnits() << " user units" << endl;
		cout << " List of " << numCell << " cells:" << endl;
		for (size_t indi = 0; indi < numCell; indi++)
		{
			cout << "  " << indi + 1 << ". " << ((this->cells)[indi]).getCellName() << endl;
			cout << "   Counts: " << (this->cells)[indi].getNumBound() << " boundaries, " << (this->cells)[indi].getNumPath() << " paths, " << (this->cells)[indi].getNumNode() << " nodes, " << (this->cells)[indi].getNumBox() << " boxes," << endl << "     " << (this->cells)[indi].getNumText() << " text boxes, and " << (this->cells)[indi].getNumSRef() << " structure references" << endl;
		}
		for (size_t indi = 0; indi < indCellPrint.size(); indi++)
		{
			(this->cells)[indCellPrint[indi]].print();
		}
		cout << "------" << endl;
	}

	///////////////////// required callbacks /////////////////////
	/// @brief bit array callback 
	/// @param ascii_record_type record 
	/// @param ascii_data_type data type 
	/// @param vBitArray data array  
	virtual void bit_array_cbk(const char* ascii_record_type, const char* ascii_data_type, vector<int> const& vBitArray)
	{
		//cout << __func__ << endl;
		this->general_cbk(ascii_record_type, ascii_data_type, vBitArray);
	}
	/// @brief 2-byte integer callback 
	/// @param ascii_record_type record 
	/// @param ascii_data_type data type 
	/// @param vInteger data array  
	virtual void integer_2_cbk(const char* ascii_record_type, const char* ascii_data_type, vector<int> const& vInteger)
	{
		//cout << __func__ << endl;
		this->general_cbk(ascii_record_type, ascii_data_type, vInteger);
	}
	/// @brief 4-byte integer callback 
	/// @param ascii_record_type record 
	/// @param ascii_data_type data type 
	/// @param vInteger data array  
	virtual void integer_4_cbk(const char* ascii_record_type, const char* ascii_data_type, vector<int> const& vInteger)
	{
		//cout << __func__ << endl;
		this->general_cbk(ascii_record_type, ascii_data_type, vInteger);
	}
	/// @brief 4-byte floating point number callback 
	/// @param ascii_record_type record 
	/// @param ascii_data_type data type 
	/// @param vFloat data array  
	virtual void real_4_cbk(const char* ascii_record_type, const char* ascii_data_type, vector<double> const& vFloat)
	{
		//cout << __func__ << endl;
		this->general_cbk(ascii_record_type, ascii_data_type, vFloat);
	}
	/// @brief 8-byte floating point number callback 
	/// @param ascii_record_type record 
	/// @param ascii_data_type data type 
	/// @param vFloat data array  
	virtual void real_8_cbk(const char* ascii_record_type, const char* ascii_data_type, vector<double> const& vFloat)
	{
		//cout << __func__ << endl;
		this->general_cbk(ascii_record_type, ascii_data_type, vFloat);
	}
	/// @brief string callback 
	/// @param ascii_record_type record 
	/// @param ascii_data_type data type 
	/// @param str data 
	virtual void string_cbk(const char* ascii_record_type, const char* ascii_data_type, string const& str)
	{
		//cout << __func__ << endl;
		this->general_cbk(ascii_record_type, ascii_data_type, str);
	}
	/// @brief begin or end indicator of a block 
	/// @param ascii_record_type record 
	virtual void begin_end_cbk(const char* ascii_record_type)
	{
		//cout << __func__ << endl;
		this->general_cbk(ascii_record_type, "", vector<int>(0));
	}

	/// @brief A generic callback function handles all other callback functions. 
	/// It is not efficient but concise as a demo. 
	/// @tparam ContainerType container type 
	/// @param ascii_record_type record 
	/// @param ascii_data_type data type 
	/// @param data data values 
	template <typename ContainerType>
	void general_cbk(string const& ascii_record_type, string const& ascii_data_type, ContainerType const& data)
	{
		// Data Printing
		/*if (this->getElement() == 's')
		{
		cout << "ascii_record_type: " << ascii_record_type << endl
		<< "ascii_data_type: " << ascii_data_type << endl
		<< "data size: " << data.size() << endl;
		}*/

		// Data handling
		if (ascii_record_type == "HEADER")
		{
			this->version = std::to_string(data[0]);
		}
		else if (ascii_record_type == "BGNLIB")
		{
			this->dateMod = std::to_string(data[0]) + "-" + std::to_string(data[1]) + "-" + std::to_string(data[2]) + " " + std::to_string(data[3]) + ":" + std::to_string(data[4]) + ":" + std::to_string(data[5]);
			this->dateAccess = std::to_string(data[6]) + "-" + std::to_string(data[7]) + "-" + std::to_string(data[8]) + " " + std::to_string(data[9]) + ":" + std::to_string(data[10]) + ":" + std::to_string(data[11]);
		}
		else if (ascii_record_type == "LIBNAME")
		{
			// Fix and print library name
			char label[64];
			for (size_t indj = 0; indj < data.size(); indj++) // Only store printable characters
			{
				if (((int)data[indj] < 32) || ((int)data[indj] > 128))
				{
					label[indj] = '\0';
					break; // Break if unprintable character
				}
				label[indj] = data[indj];
			}
			//cout << "Name of this library: " << label << endl;

			// Store library name
			this->libName = label;
		}
		else if (ascii_record_type == "UNITS")
		{
			this->dbUserUnits = data[0];
			this->dbUnits = data[1];
		}
		else if (ascii_record_type == "BGNSTR")
		{
			(this->cells).push_back(GeoCell()); // Create new geometric cell in vector
			((this->cells)[this->numCell]).dateCreate = std::to_string(data[0]) + "-" + std::to_string(data[1]) + "-" + std::to_string(data[2]) + " " + std::to_string(data[3]) + ":" + std::to_string(data[4]) + ":" + std::to_string(data[5]); // Save cell creation date
			((this->cells)[this->numCell]).dateMod = std::to_string(data[6]) + "-" + std::to_string(data[7]) + "-" + std::to_string(data[8]) + " " + std::to_string(data[9]) + ":" + std::to_string(data[10]) + ":" + std::to_string(data[11]); // Save last modification date of cell
		}
		else if (ascii_record_type == "STRNAME")
		{
			// Print and store name
			char label[64];
			for (size_t indj = 0; indj < data.size(); indj++) // Only store printable characters
			{
				if (((int)data[indj] < 32) || ((int)data[indj] > 128))
				{
					label[indj] = '\0';
					break; // Break if unprintable character
				}
				label[indj] = data[indj];
			}
			//cout << "Geometric cell name: " << label << endl;
			((this->cells)[this->numCell]).cellName = label; // Save cell name
		}
		else if (ascii_record_type == "BOUNDARY")
		{
			(this->element) = 'b'; // Present working element designated as boundary
			((this->cells)[this->numCell]).boundaries.emplace_back(boundary()); // Record an empty boundary for now
		}
		else if (ascii_record_type == "PATH")
		{
			(this->element) = 'p';
			((this->cells)[this->numCell]).paths.emplace_back(path());
		}
		else if (ascii_record_type == "NODE")
		{
			(this->element) = 'n';
			((this->cells)[this->numCell]).nodes.emplace_back(node());
		}
		else if (ascii_record_type == "BOX")
		{
			(this->element) = 'x';
			((this->cells)[this->numCell]).boxes.emplace_back(box());
		}
		else if (ascii_record_type == "TEXT")
		{
			(this->element) = 't';
			((this->cells)[this->numCell]).textboxes.emplace_back(textbox());
		}
		else if (ascii_record_type == "SREF")
		{
			(this->element) = 's';
			((this->cells)[this->numCell]).sreferences.emplace_back(sref());
		}
		else if (ascii_record_type == "LAYER")
		{
			if (this->getElement() == 'b')
			{
				boundary modBound = getCell(this->numCell).boundaries.back(); // Get copy of boundary
				((this->cells)[this->numCell]).boundaries.pop_back(); // Remove last boundary vector entry
				((this->cells)[this->numCell]).boundaries.emplace_back(boundary(modBound.getBounds(), data[0], modBound.getProps())); // Put boundary back with layer update
			}
			else if (this->getElement() == 'p')
			{
				path modPath = getCell(this->numCell).paths.back();
				((this->cells)[this->numCell]).paths.pop_back();
				((this->cells)[this->numCell]).paths.emplace_back(path(modPath.getPaths(), data[0], modPath.getProps(), modPath.getType(), modPath.getWidth()));
			}
			else if (this->getElement() == 'n')
			{
				node modNode = getCell(this->numCell).nodes.back();
				((this->cells)[this->numCell]).nodes.pop_back();
				((this->cells)[this->numCell]).nodes.emplace_back(node(modNode.getNodes(), data[0], modNode.getProps(), modNode.getType()));
			}
			else if (this->getElement() == 'x')
			{
				box modBox = getCell(this->numCell).boxes.back();
				((this->cells)[this->numCell]).boxes.pop_back();
				((this->cells)[this->numCell]).boxes.emplace_back(box(modBox.getBoxes(), data[0], modBox.getProps(), modBox.getType()));
			}
			else if (this->getElement() == 't')
			{
				textbox modText = getCell(this->numCell).textboxes.back();
				((this->cells)[this->numCell]).textboxes.pop_back();
				((this->cells)[this->numCell]).textboxes.emplace_back(textbox(modText.getTexts(), data[0], modText.getProps(), modText.getType(), modText.getFontID(), modText.getJusts(), modText.getWidth(), modText.getTextStr()));
			}
		}
		else if (ascii_record_type == "XY")
		{
			// Scaling for coordinates
			double unitFactor = this->getdbUnits();

			// Print ordered pairs
			/*cout << "    " << std::scientific;
			cout.precision(4);
			for (size_t indj = 0; indj < data.size(); indj++) // Direct stream write for each ordered pair
			{
			cout << "(" << unitFactor * data[indj++] << ", " << unitFactor * data[indj + 1] << ") ";
			}
			cout << endl;*/

			// Calculate coordinates
			vector<double>* coord = new vector<double>;
			for (size_t indi = 0; indi < data.size(); indi++)
			{
				(*coord).push_back(unitFactor * data[indi]);
			}

			// Store coordinates
			if (this->getElement() == 'b')
			{
				boundary modBound = getCell(this->numCell).boundaries.back(); // Get copy of boundary
				((this->cells)[this->numCell]).boundaries.pop_back(); // Remove last boundary vector entry
				((this->cells)[this->numCell]).boundaries.emplace_back(boundary(coord, modBound.getLayer(), modBound.getProps())); // Put boundary back with coordinate update
			}
			else if (this->getElement() == 'p')
			{
				path modPath = getCell(this->numCell).paths.back();
				((this->cells)[this->numCell]).paths.pop_back();
				((this->cells)[this->numCell]).paths.emplace_back(path(coord, modPath.getLayer(), modPath.getProps(), modPath.getType(), modPath.getWidth()));
			}
			else if (this->getElement() == 'n')
			{
				node modNode = getCell(this->numCell).nodes.back();
				((this->cells)[this->numCell]).nodes.pop_back();
				((this->cells)[this->numCell]).nodes.emplace_back(node(coord, modNode.getLayer(), modNode.getProps(), modNode.getType()));
			}
			else if (this->getElement() == 'x')
			{
				box modBox = getCell(this->numCell).boxes.back();
				((this->cells)[this->numCell]).boxes.pop_back();
				((this->cells)[this->numCell]).boxes.emplace_back(box(coord, modBox.getLayer(), modBox.getProps(), modBox.getType()));
			}
			else if (this->getElement() == 't')
			{
				textbox modText = getCell(this->numCell).textboxes.back();
				((this->cells)[this->numCell]).textboxes.pop_back();
				((this->cells)[this->numCell]).textboxes.emplace_back(textbox(*coord, modText.getLayer(), modText.getProps(), modText.getType(), modText.getFontID(), modText.getJusts(), modText.getWidth(), modText.getTextStr()));
			}
			else if (this->getElement() == 's')
			{
				sref modSRef = getCell(this->numCell).sreferences.back();
				((this->cells)[this->numCell]).sreferences.pop_back();
				((this->cells)[this->numCell]).sreferences.emplace_back(sref(*coord, modSRef.getProps(), modSRef.getSRefName()));
			}

			// Delete coordinates variable
			delete coord;
		}
		else if (ascii_record_type == "DATATYPE") // Unimplemented in the GDSII standard
		{
			//cout << "Datatype for " << (this->element) << ": " << data[0] << endl;
		}
		else if (ascii_record_type == "PATHTYPE")
		{
			if (this->getElement() == 'p')
			{
				path modPath = getCell(this->numCell).paths.back();
				((this->cells)[this->numCell]).paths.pop_back();
				((this->cells)[this->numCell]).paths.emplace_back(path(modPath.getPaths(), modPath.getLayer(), modPath.getProps(), data[0], modPath.getWidth()));
			}
		}
		else if (ascii_record_type == "NODETYPE")
		{
			if (this->getElement() == 'n')
			{
				node modNode = getCell(this->numCell).nodes.back();
				((this->cells)[this->numCell]).nodes.pop_back();
				((this->cells)[this->numCell]).nodes.emplace_back(node(modNode.getNodes(), modNode.getLayer(), modNode.getProps(), data[0]));
			}
		}
		else if (ascii_record_type == "BOXTYPE")
		{
			if (this->getElement() == 'x')
			{
				box modBox = getCell(this->numCell).boxes.back();
				((this->cells)[this->numCell]).boxes.pop_back();
				((this->cells)[this->numCell]).boxes.emplace_back(box(modBox.getBoxes(), modBox.getLayer(), modBox.getProps(), data[0]));
			}
		}
		else if (ascii_record_type == "TEXTTYPE")
		{
			if (this->getElement() == 't')
			{
				textbox modText = getCell(this->numCell).textboxes.back();
				((this->cells)[this->numCell]).textboxes.pop_back();
				((this->cells)[this->numCell]).textboxes.emplace_back(textbox(modText.getTexts(), modText.getLayer(), modText.getProps(), data[0], modText.getFontID(), modText.getJusts(), modText.getWidth(), modText.getTextStr()));
			}
		}
		else if (ascii_record_type == "WIDTH")
		{
			// Scaling for widths
			double unitFactor = this->getdbUnits();

			// Store widths
			if (this->getElement() == 'p')
			{
				path modPath = getCell(this->numCell).paths.back();
				((this->cells)[this->numCell]).paths.pop_back();
				((this->cells)[this->numCell]).paths.emplace_back(path(modPath.getPaths(), modPath.getLayer(), modPath.getProps(), modPath.getType(), unitFactor * data[0]));
			}
		}
		else if (ascii_record_type == "PRESENTATION")
		{
			// Declare binary masks
			int horizJustMask = 0b0000000000000011;
			int vertJustMask = horizJustMask << 2;
			int fontIDMask = vertJustMask << 2;

			// Store text presentation information
			if (this->getElement() == 't')
			{
				// Get existing text box record
				textbox modText = getCell(this->numCell).textboxes.back();
				((this->cells)[this->numCell]).textboxes.pop_back();

				// Find font identifier based on bit pattern
				int fontID = ((int)data[0] & fontIDMask) >> 4;

				// Find vertical and horizontal justification based on bit pattern
				int vertJust = ((int)data[0] & vertJustMask) >> 2;
				int horizJust = ((int)data[0] & horizJustMask);
				if (vertJust == 4)
				{
					cout << "Illegal text box vertical justification reset to middle" << endl;
					vertJust = 1;
				}
				if (horizJust == 4)
				{
					cout << "Illegal text box horizontal justification reset to center" << endl;
					vertJust = 1;
				}

				// Put text box back with presentation update
				((this->cells)[this->numCell]).textboxes.emplace_back(textbox(modText.getTexts(), modText.getLayer(), modText.getProps(), modText.getType(), fontID, { vertJust, horizJust }, modText.getWidth(), modText.getTextStr()));
			}
		}
		else if (ascii_record_type == "STRANS") // Unimplemented by user until needed
		{
			//cout << "Linear transformation of element: " << data[0] << endl;
		}
		else if (ascii_record_type == "STRING")
		{
			// Print and store text box string
			char label[64];
			for (size_t indj = 0; indj < data.size(); indj++) // Only store printable characters
			{
				if (((int)data[indj] < 32) || ((int)data[indj] > 128))
				{
					label[indj] = '\0';
					break; // Break if unprintable character
				}
				label[indj] = data[indj];
			}
			//cout << "Text box string: " << label << endl;
			if (this->getElement() == 't')
			{
				textbox modText = getCell(this->numCell).textboxes.back();
				((this->cells)[this->numCell]).textboxes.pop_back();
				((this->cells)[this->numCell]).textboxes.emplace_back(textbox(modText.getTexts(), modText.getLayer(), modText.getProps(), modText.getType(), modText.getFontID(), modText.getJusts(), modText.getWidth(), label));
			}
		}
		else if (ascii_record_type == "SNAME")
		{
			// Print and store name
			char label[64];
			for (size_t indj = 0; indj < data.size(); indj++) // Only store printable characters
			{
				if (((int)data[indj] < 32) || ((int)data[indj] > 128))
				{
					label[indj] = '\0';
					break; // Break if unprintable character
				}
				label[indj] = data[indj];
			}
			//cout << "Structure reference name: " << label << endl;
			if (this->getElement() == 's')
			{
				sref modSRef = getCell(this->numCell).sreferences.back();
				((this->cells)[this->numCell]).sreferences.pop_back();
				((this->cells)[this->numCell]).sreferences.emplace_back(sref(modSRef.getSRefs(), modSRef.getProps(), label));
			}
		}
		else if (ascii_record_type == "PROPATTR")
		{
			this->numProp = data[0]; // Store position for attribute value to be saved
		}
		else if (ascii_record_type == "PROPVALUE")
		{
			// Fix and print property
			char label[64];
			for (size_t indj = 0; indj < data.size(); indj++) // Only store printable characters
			{
				if (((int)data[indj] < 32) || ((int)data[indj] > 128))
				{
					label[indj] = '\0';
					break; // Break if unprintable character
				}
				label[indj] = data[indj];
			}
			//cout << "Element property value: " << label << endl;

			// Store new property (might only handle a single one)
			if (this->getElement() == 'b')
			{
				boundary modBound = getCell(this->numCell).boundaries.back(); // Get copy of boundary
				if (modBound.getProps() != NULL)
				{
					vector<std::string> oldProps = *(modBound.getProps()); // Get properties of this boundary (to be used)
				}
				else
				{
					vector<std::string> oldProps = {}; // Use empty string
				}
				vector<std::string>* modProps = new vector<std::string>; // Create dynamic string vector for changed property
				(*modProps).assign(this->numProp, ""); // Ensure enough elements available for property to go in correct location
				(*modProps)[this->numProp - 1] = label; // Put new property value in place
				((this->cells)[this->numCell]).boundaries.pop_back(); // Remove last boundary vector entry
				((this->cells)[this->numCell]).boundaries.emplace_back(boundary(modBound.getBounds(), modBound.getLayer(), modProps)); // Put boundary back with properties update
				delete modProps; // Free pointer to vector
			}
			else if (this->getElement() == 'p')
			{
				path modPath = getCell(this->numCell).paths.back();
				if (modPath.getProps() != NULL)
				{
					vector<std::string> oldProps = *(modPath.getProps());
				}
				else
				{
					vector<std::string> oldProps = {};
				}
				vector<std::string>* modProps = new vector<std::string>;
				(*modProps).assign(this->numProp, "");
				(*modProps)[this->numProp - 1] = label;
				((this->cells)[this->numCell]).paths.pop_back();
				((this->cells)[this->numCell]).paths.emplace_back(path(modPath.getPaths(), modPath.getLayer(), modProps, modPath.getType(), modPath.getWidth()));
				delete modProps; // Free pointer to vector
			}
			else if (this->getElement() == 'n')
			{
				node modNode = getCell(this->numCell).nodes.back();
				if (modNode.getProps() != NULL)
				{
					vector<std::string> oldProps = *(modNode.getProps());
				}
				else
				{
					vector<std::string> oldProps = {};
				}
				vector<std::string>* modProps = new vector<std::string>;
				(*modProps).assign(this->numProp, "");
				(*modProps)[this->numProp - 1] = label;
				((this->cells)[this->numCell]).nodes.pop_back();
				((this->cells)[this->numCell]).nodes.emplace_back(node(modNode.getNodes(), modNode.getLayer(), modProps, modNode.getType()));
				delete modProps; // Free pointer to vector
			}
			else if (this->getElement() == 'x')
			{
				box modBox = getCell(this->numCell).boxes.back();
				if (modBox.getProps() != NULL)
				{
					vector<std::string> oldProps = *(modBox.getProps());
				}
				else
				{
					vector<std::string> oldProps = {};
				}
				vector<std::string>* modProps = new vector<std::string>;
				(*modProps).assign(this->numProp, "");
				(*modProps)[this->numProp - 1] = label;
				((this->cells)[this->numCell]).boxes.pop_back();
				((this->cells)[this->numCell]).boxes.emplace_back(box(modBox.getBoxes(), modBox.getLayer(), modProps, modBox.getType()));
				delete modProps; // Free pointer to vector
			}
			else if (this->getElement() == 't')
			{
				textbox modText = getCell(this->numCell).textboxes.back();
				if (modText.getProps() != NULL)
				{
					vector<std::string> oldProps = *(modText.getProps());
				}
				else
				{
					vector<std::string> oldProps = {};
				}
				vector<std::string>* modProps = new vector<std::string>;
				(*modProps).assign(this->numProp, "");
				(*modProps)[this->numProp - 1] = label;
				((this->cells)[this->numCell]).textboxes.pop_back();
				((this->cells)[this->numCell]).textboxes.emplace_back(textbox(modText.getTexts(), modText.getLayer(), modProps, modText.getType(), modText.getFontID(), modText.getJusts(), modText.getWidth(), modText.getTextStr()));
				delete modProps; // Free pointer to vector
			}
			else if (this->getElement() == 's')
			{
				sref modSRef = getCell(this->numCell).sreferences.back();
				if (modSRef.getProps() != NULL)
				{
					vector<std::string> oldProps = *(modSRef.getProps());
				}
				else
				{
					vector<std::string> oldProps = {};
				}
				vector<std::string>* modProps = new vector<std::string>;
				(*modProps).assign(this->numProp, "");
				(*modProps)[this->numProp - 1] = label;
				((this->cells)[this->numCell]).sreferences.pop_back();
				((this->cells)[this->numCell]).sreferences.emplace_back(sref(modSRef.getSRefs(), modProps, modSRef.getSRefName()));
				delete modProps; // Free pointer to vector
			}
		}
		else if (ascii_record_type == "ENDEL")
		{
			(this->element) = '\0'; // Reset current element
		}
		else if (ascii_record_type == "ENDSTR")
		{
			(this->numCell)++; // Increment number of cells to use as index
		}
		else if (ascii_record_type == "ENDLIB")
		{
			//cout << "Entire GDSII file read" << endl; // Nothing more to do
		}
		else // Unhandled
		{
			cout << "ascii_record_type: " << ascii_record_type << endl
				<< "ascii_data_type: " << ascii_data_type << endl
				<< "data size: " << data.size() << endl;
		}
	}

	/// @brief destructor
	~AsciiDataBase()
	{
		//~(this->cells);
	}
};

/// @brief test enum callbacks
struct EnumDataBase : public GdsParser::GdsDataBaseKernel
{
	/// @brief constructor 
	EnumDataBase()
	{
		//cout << "constructing EnumDataBase" << endl;
	}
	///////////////////// required callbacks /////////////////////
	/// @brief bit array callback 
	/// @param record_type record 
	/// @param data_type data type 
	/// @param vBitArray data array  
	virtual void bit_array_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, vector<int> const& vBitArray)
	{
		//cout << __func__ << endl;
		this->general_cbk(record_type, data_type, vBitArray);
	}
	/// @brief 2-byte integer callback 
	/// @param record_type record 
	/// @param data_type data type 
	/// @param vInteger data array  
	virtual void integer_2_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, vector<int> const& vInteger)
	{
		//cout << __func__ << endl;
		this->general_cbk(record_type, data_type, vInteger);
	}
	/// @brief 4-byte integer callback 
	/// @param record_type record 
	/// @param data_type data type 
	/// @param vInteger data array  
	virtual void integer_4_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, vector<int> const& vInteger)
	{
		//cout << __func__ << endl;
		this->general_cbk(record_type, data_type, vInteger);
	}
	/// @brief 4-byte floating point number callback 
	/// @param record_type record 
	/// @param data_type data type 
	/// @param vFloat data array  
	virtual void real_4_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, vector<double> const& vFloat)
	{
		//cout << __func__ << endl;
		this->general_cbk(record_type, data_type, vFloat);
	}
	/// @brief 8-byte floating point number callback 
	/// @param record_type record 
	/// @param data_type data type 
	/// @param vFloat data array  
	virtual void real_8_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, vector<double> const& vFloat)
	{
		//cout << __func__ << endl;
		this->general_cbk(record_type, data_type, vFloat);
	}
	/// @brief string callback 
	/// @param record_type record 
	/// @param data_type data type 
	/// @param str data 
	virtual void string_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, string const& str)
	{
		//cout << __func__ << endl;
		this->general_cbk(record_type, data_type, str);
	}
	/// @brief begin or end indicator of a block 
	/// @param record_type record 
	virtual void begin_end_cbk(GdsParser::GdsRecords::EnumType record_type)
	{
		//cout << __func__ << endl;
		this->general_cbk(record_type, GdsParser::GdsData::NO_DATA, vector<int>(0));
	}

	/// @brief A generic callback function handles all other callback functions. 
	/// It is not efficient but concise as a demo. 
	/// @tparam ContainerType container type 
	/// @param record_type record 
	/// @param data_type data type 
	/// @param data data values 
	template <typename ContainerType>
	void general_cbk(GdsParser::GdsRecords::EnumType record_type, GdsParser::GdsData::EnumType data_type, ContainerType const& data)
	{
		/*cout << "ascii_record_type: " << GdsParser::gds_record_ascii(record_type) << endl
		<< "ascii_data_type: " << GdsParser::gds_data_ascii(data_type) << endl
		<< "data size: " << data.size() << endl;*/
		switch (record_type)
		{
		case GdsParser::GdsRecords::UNITS:
			break;
		case GdsParser::GdsRecords::BOUNDARY:
			break;
		case GdsParser::GdsRecords::LAYER:
			/*cout << "LAYER = " << data[0] <<  endl;*/
			break;
		case GdsParser::GdsRecords::XY:
			/*for (typename ContainerType::const_iterator it = data.begin(); it != data.end(); ++it)
			cout << *it << " ";
			cout << endl;
			cout << data.size() << endl;*/
			break;
		case GdsParser::GdsRecords::ENDEL:
			break;
		default:
			break;
		}
	}
};
#endif