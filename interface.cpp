// Geometry Implementation

// Used for cout << "string"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "interface.h"
using namespace std;

// Constructor
GeoCell::GeoCell(string cellName, string metadata, vector<unsigned long> boundLayers, vector<vector<double>*> bounds, vector<vector<string>*> boundProps, vector<unsigned long> pathLayers, vector<unsigned long> pathTypes, vector<double> pathWidths, vector<vector<double>*> paths, vector<vector<string>*> pathProps, vector<string> srefNames, vector<vector<double>> srefs, vector<vector<string>*> srefProps, vector<unsigned long> nodeLayers, vector<unsigned long> nodeTypes, vector<vector<double>*> nodes, vector<vector<string>*> nodeProps, vector<unsigned long> textLayers, vector<double> textWidths, vector<vector<unsigned long>> textJusts, vector<vector<double>> texts, vector<string> textStr)
{
	this->cellName = cellName;
	this->metadata = metadata;
	this->boundLayers = boundLayers;
	this->bounds = bounds;
	this->boundProps = boundProps;
	this->pathLayers = pathLayers;
	this->pathTypes = pathTypes;
	this->pathWidths = pathWidths;
	this->paths = paths;
	this->pathProps = pathProps;
	this->srefNames = srefNames;
	this->srefs = srefs;
	this->srefProps = srefProps;
	this->nodeLayers = nodeLayers;
	this->nodeTypes = nodeTypes;
	this->nodes = nodes;
	this->nodeProps = nodeProps;
	this->textLayers = textLayers;
	this->textWidths = textWidths;
	this->textJusts = textJusts;
	this->texts = texts;
	this->textStr = textStr;
}

// Modify a geometric cell
bool GeoCell::modCell()
{
	// To be implemted...
    return true;
}

// Get the name of the cell
string GeoCell::getCellName() const
{
    return this->cellName;
}

// Get the metadata of the cell
string GeoCell::getCellMeta() const
{
    return this->metadata;
}

// Get number of boundaries
int GeoCell::getNumBound() const
{
    return bounds.size();
}

// Get number of paths
int GeoCell::getNumPath() const
{
    return paths.size();
}

// Get number of structure references
int GeoCell::getNumSRef() const
{
    return srefs.size();
}

// Get number of nodes
int GeoCell::getNumNode() const
{
	return nodes.size();
}

// Get number of text boxes
int GeoCell::getNumText() const
{
	return texts.size();
}

// Return boundary layer numbers
vector<unsigned long> GeoCell::getBoundLayers() const
{
    return this->boundLayers;
}

// Return boundaries in cell coordinates
vector<vector<double>*> GeoCell::getBounds() const
{
    return this->bounds;
}

// Return boundary properties
vector<vector<std::string>*> GeoCell::getBoundProps() const
{
    return this->boundProps;
}

// Return path layer numbers
vector<unsigned long> GeoCell::getPathLayers() const
{
    return this->pathLayers;
}

// Return path types
// 0 = square ends at vertices, 1 = round ends, 2 = square ends overshoot vertices by half width
vector<unsigned long> GeoCell::getPathTypes() const
{
    return this->pathTypes;
}

// Return path widths
// Negative values mean independent of structure scaling
vector<double> GeoCell::getPathWidths() const
{
    return this->pathWidths;
}

// Return paths in cell coordinates
vector<vector<double>*> GeoCell::getPaths() const
{
    return this->paths;
}

// Return path properties
vector<vector<std::string>*> GeoCell::getPathProps() const
{
    return this->pathProps;
}

// Return structure reference names
vector<std::string> GeoCell::getSRefNames() const
{
    return this->srefNames;
}

// Return coordinate pairs of placed structure references
vector<vector<double>> GeoCell::getSRefs() const
{
    return this->srefs;
}

// Return structure reference properties
vector<vector<std::string>*> GeoCell::getSRefProps() const
{
    return this->srefProps;
}

// Return node layer numbers
vector<unsigned long> GeoCell::getNodeLayers() const
{
	return this->nodeLayers;
}

// Return node types
vector<unsigned long> GeoCell::getNodeTypes() const
{
	return this->nodeTypes;
}

// Return nodes in cell coordinates
vector<vector<double>*> GeoCell::getNodes() const
{
	return this->nodes;
}

// Return node properties
vector<vector<std::string>*> GeoCell::getNodeProps() const
{
	return this->nodeProps;
}

// Return text layer numbers
vector<unsigned long> GeoCell::getTextLayers() const
{
    return this->textLayers;
}

// Return text widths
// Negative values mean independent of structure scaling
vector<double> GeoCell::getTextWidths() const
{
    return this->textWidths;
}

// Return text justifications
vector<vector<unsigned long>> GeoCell::getTextJusts() const
{
    return this->textJusts;
}

// Return coordinate pairs of placed text boxes
vector<vector<double>> GeoCell::getTexts() const
{
    return this->texts;
}

// Return text in text boxes
vector<std::string> GeoCell::getTextStr() const
{
    return this->textStr;
}

// Print the cell contents
void GeoCell::print() const
{
    int numBound = GeoCell::getNumBound();
    int numPath = GeoCell::getNumPath();
    int numSRef = GeoCell::getNumSRef();
	int numNode = GeoCell::getNumNode();
	int numText = GeoCell::getNumText();

    cout << " ------" << endl;
    cout << " Cell Details:" << endl;
    cout << "  Cell Name: " << this->cellName << endl;
    cout << "  Metadata: " << this->metadata << endl;
    cout << "  List of " << numBound << " boundaries:" << endl;
    for (size_t indi = 0; indi < numBound; indi++) // Handle each boundary
    {
		if (boundProps[indi]) // Check if NULL pointer
		{
            cout << "   #" << indi + 1 << " on Layer " << boundLayers[indi] << " with " << (*(boundProps[indi]))[0] << " as Property " << 1 << endl;
		}
        else
        {
            cout << "   #" << indi + 1 << " on Layer " << boundLayers[indi] << endl;
        }
        char point[128];
        string strPoints;
        for (size_t indj = 0; indj < (*(bounds[indi])).size(); indj++) // C string for each ordered pair
        {
            sprintf(point, "(%1.4g, %1.4g) ", (*(bounds[indi]))[indj++], (*(bounds[indi]))[indj + 1]);
            strPoints.append(point);
        }
        cout << "    " << strPoints << endl; // Print line of ordered pairs
    }
	cout << "  List of " << numPath << " paths:" << endl;
	for (size_t indi = 0; indi < numPath; indi++) // Handle each path
	{
        if (pathProps[indi]) // Check if NULL pointer
        {
            cout << "   #" << indi + 1 << " on Layer " << pathLayers[indi] << " with path type #" << pathTypes[indi] << " and width " << pathWidths[indi] << " with " << (*(pathProps[indi]))[0] << " as Property " << 1 << endl;
        }
        else
        {
            cout << "   #" << indi + 1 << " on Layer " << pathLayers[indi] << " with path type #" << pathTypes[indi] << " and width " << pathWidths[indi] << endl;
        }
        char point[128];
        string strPoints;
        for (size_t indj = 0; indj < (*(paths[indi])).size(); indj++) // C string for each ordered pair
        {
            sprintf(point, "(%1.4g, %1.4g) ", (*(paths[indi]))[indj++], (*(paths[indi]))[indj + 1]);
            strPoints.append(point);
        }
        cout << "    " << strPoints << endl; // Print line of ordered pairs
	}
	cout << "  List of " << numSRef << " structure references:" << endl;
	for (size_t indi = 0; indi < numSRef; indi++) // Handle each structure
	{
        if (srefProps[indi]) // Check if NULL pointer
        {
            cout << "   #" << indi + 1 << " named " << srefNames[indi] << " with " << (*(srefProps[indi]))[0] << " as Property " << 1 << endl;
        }
        else
        {
            cout << "   #" << indi + 1 << " named " << srefNames[indi] << endl;
        }
        char point[128];
        sprintf(point, "(%1.4g, %1.4g)", srefs[indi][0], srefs[indi][1]);
        cout << "    " << point << endl;
	}
	cout << "  List of " << numNode << " nodes:" << endl;
	for (size_t indi = 0; indi < numNode; indi++) // Handle each node
	{
		if (nodeProps[indi]) // Check if NULL pointer
		{
			cout << "   #" << indi + 1 << " on Layer " << nodeLayers[indi] << " with node type #" << nodeTypes[indi] << " with " << (*(nodeProps[indi]))[0] << " as Property " << 1 << endl;
		}
		else
		{
			cout << "   #" << indi + 1 << " on Layer " << nodeLayers[indi] << " with node type #" << nodeTypes[indi] << endl;
		}
		char point[128];
		string strPoints;
		for (size_t indj = 0; indj < (*(nodes[indi])).size(); indj++) // C string for each ordered pair
		{
			sprintf(point, "(%1.4g, %1.4g) ", (*(nodes[indi]))[indj++], (*(nodes[indi]))[indj + 1]);
			strPoints.append(point);
		}
		cout << "    " << strPoints << endl; // Print line of ordered pairs
	}
    cout << " ------" << endl;
}

// Destructor
GeoCell::~GeoCell()
{
    /*for (size_t indi = 0; indi < this->getNumBound(); indi++)
    {
        cout << this->cellName << ": " << this->bounds[indi] << " (" << (*(bounds[indi]))[0] << ", " << (*(bounds[indi]))[1] << endl;
        cout << indi << " out of " << this->getNumBound() << endl;
        //delete this->bounds[indi];
        //delete this->boundProps[indi];
        delete this->bounds.at(indi);
        delete this->boundProps.at(indi);
    }*/
    /*for (size_t indi = 0; indi < this->getNumPath(); indi++)
    {
        delete this->paths[indi];
        delete this->pathProps[indi];
    }
    for (size_t indi = 0; indi < this->getNumSRef(); indi++)
    {
        delete this->srefProps[indi];
    }*/
}


// Constructor
Geometry::Geometry()
{
    this->fileName = "\0";
    this->metadata = "\0";
    this->fileSize = 0;
	//this->cells = GeoCell::GeoCell(string(""), string(""), vector<vector<double>*>(), NULL);
}

// Load file contents
// Return false if unable to load.
bool Geometry::load(string fileName)
{
	// Get relative path to file
    this->fileName = fileName;
    string fullName = new char[fileName.length() + 2];
    fullName = "./" + fileName;

	// Attempt to open GDT file
	ifstream gdtFile (fullName.c_str());
    if (gdtFile.is_open())
	{
		// File is readable line-by-line
	    string fileLine;
		getline(gdtFile, fileLine);

		// Save first and subsequent lines before comments as metadata
		if (fileLine[0] != '#')
		{
		    this->metadata = fileLine;
			getline(gdtFile, fileLine);
		}
		while (fileLine[0] != '#')
		{
		    this->metadata.append("\n" + fileLine);
			getline(gdtFile, fileLine);
		}

		// Go past the block of comments
		while (fileLine[0] == '#')
		{
		    getline(gdtFile, fileLine);
		}

		// Read file cell-by-cell
		while (!gdtFile.eof())
		{
			// Identify cell
			if (fileLine.compare(0, 5, "cell{") == 0)
			{
				// Describe the cell
				size_t indID = fileLine.find("'");
				string cellMetadata = fileLine.substr(5, indID - 5);
				string cellName = fileLine.substr(indID + 1, fileLine.length() - cellMetadata.length() - 7); // 5 characters for cell syntax and 2 to remove trailing quote

				// Read contents to end of cell
				vector<unsigned long> boundLayers;
				vector<vector<double>*> boundsAll;
				vector<vector<string>*> boundProps;
				vector<unsigned long> pathLayers;
				vector<unsigned long> pathTypes;
				vector<double> pathWidths;
				vector<vector<double>*> pathsAll;
				vector<vector<string>*> pathProps;
				vector<string> srefNames;
				vector<vector<double>> srefsAll;
				vector<vector<string>*> srefProps;
				vector<unsigned long> nodeLayers;
				vector<unsigned long> nodeTypes;
				vector<vector<double>*> nodesAll;
				vector<vector<string>*> nodeProps;
				vector<unsigned long> textLayers;
				vector<double> textWidths;
				vector<vector<unsigned long>> textJusts;
				vector<vector<double>> textsAll;
				vector<string> textStr;
				while (fileLine[0] != '}')
				{
                    getline(gdtFile, fileLine);
				    // Handle boundaries
				    if (fileLine.compare(0, 2, "b{") == 0)
				    {
				        // Extract layer number
				        size_t indLayer = fileLine.find(" ") - 1;
				        boundLayers.push_back(std::stoul(fileLine.substr(2, indLayer)));
				        
				        // Split string at spaces between successive coordinates
				        size_t indCoordStart = fileLine.find("xy(") + 3;
				        size_t indCoordStop = fileLine.find(")") - 1;
				        vector<double>* coord = new vector<double>;
				        size_t coordBegin = indCoordStart;
				        while (coordBegin <= indCoordStop)
				        {
				            size_t coordEnd = fileLine.find(" ", coordBegin);
				            if (coordEnd == string::npos)
				            {
				                coordEnd = indCoordStop; // Final coordinate
				            }
				            else
				            {
				                coordEnd -= 1; // Remove trailing whitespace
				            }
				            (*coord).push_back(std::stod(fileLine.substr(coordBegin, coordEnd)));
				            coordBegin = coordEnd + 2;
				        }
				        boundsAll.push_back(coord);
						//delete coord;

				        // Search for boundary properties (might only handle one...)
				        size_t indProp = fileLine.find("pr{");
				        if (indProp == string::npos)
				        {
				            boundProps.push_back(NULL);
				        }
				        else
				        {
				            size_t indPropSpace = fileLine.find(" ", indProp);
				            unsigned long propAttr = std::stoul(fileLine.substr(indProp + 3, indPropSpace - indProp));
				            vector<string>* property = new vector<string>;
				            (*property).assign(propAttr, ""); // Make sure enough elements available for property to go to correct location
				            (*property)[propAttr - 1] = fileLine.substr(indPropSpace + 2, fileLine.find("}", indProp) - indPropSpace - 3); // Remove surrounding single quotation marks
				            boundProps.push_back(property);
				        }
				    }
				    // Handle paths
				    else if (fileLine.compare(0, 2, "p{") == 0)
				    {
				        // Extract layer number
				        size_t indLayer = fileLine.find(" ") - 1;
				        pathLayers.push_back(std::stoul(fileLine.substr(2, indLayer)));
				        
				        // Extract path types
				        size_t indType = fileLine.find(" pt");
				        if (indType == string::npos)
				        {
				            pathTypes.push_back(0); // Default is square ends at vertices
				        }
				        else
				        {
				            pathTypes.push_back(std::stoul(fileLine.substr(indType + 3, fileLine.find(" ", indType + 1) - indType - 3)));
				        }

				        // Extract path widths
				        size_t indWidth = fileLine.find(" w");
				        if (indWidth == string::npos)
				        {
				            pathWidths.push_back(0);
				        }
				        else
				        {
				            pathWidths.push_back(std::stod(fileLine.substr(indWidth + 2, fileLine.find(" ", indWidth + 1) - indWidth - 2)));
				        }
				        
				        // Split string at spaces between successive coordinates
				        size_t indCoordStart = fileLine.find("xy(") + 3;
				        size_t indCoordStop = fileLine.find(")") - 1;
				        vector<double>* coord = new vector<double>;
				        size_t coordBegin = indCoordStart;
				        while (coordBegin <= indCoordStop)
				        {
				            size_t coordEnd = fileLine.find(" ", coordBegin);
				            if (coordEnd == string::npos)
				            {
				                coordEnd = indCoordStop; // Final coordinate
				            }
				            else
				            {
				                coordEnd -= 1; // Remove trailing whitespace
				            }
				            (*coord).push_back(std::stod(fileLine.substr(coordBegin, coordEnd)));
				            coordBegin = coordEnd + 2;
				        }
				        pathsAll.push_back(coord);

				        // Search for path properties (might only handle one...)
				        size_t indProp = fileLine.find("pr{");
				        if (indProp == string::npos)
				        {
				            pathProps.push_back(NULL);
				        }
				        else
				        {
				            size_t indPropSpace = fileLine.find(" ", indProp);
				            unsigned long propAttr = std::stoul(fileLine.substr(indProp + 3, indPropSpace - indProp));
                            vector<string>* property = new vector<string>;
                            (*property).assign(propAttr, ""); // Make sure enough elements available for property to go to correct location
                            (*property)[propAttr - 1] = fileLine.substr(indPropSpace + 2, fileLine.find("}", indProp) - indPropSpace - 3); // Remove surrounding single quotation marks
				            pathProps.push_back(property);
				        }
				    }
				    // Handle structure references
				    else if (fileLine.compare(0, 3, "s{'") == 0)
				    {
				        // Extract structure reference names
				        size_t indSName = fileLine.find("'", 3) - 1;
				        srefNames.push_back(fileLine.substr(3, indSName - 2));
				        
				        // Split string at space between two coordinates
				        size_t indCoordStart = fileLine.find("xy(") + 3;
				        size_t indCoordSep = fileLine.find(" ", indCoordStart);
				        size_t indCoordStop = fileLine.find(")") - 1;
				        vector<double> coord;
				        coord.assign(2, 0); // Coordinate pair initialized to origin
				        coord[0] = std::stod(fileLine.substr(indCoordStart, indCoordSep - indCoordStart));
				        coord[1] = std::stod(fileLine.substr(indCoordSep + 1, indCoordStop - indCoordSep));
				        srefsAll.push_back(coord);

				        // Search for structure properties (might only handle one...)
				        size_t indProp = fileLine.find("pr{");
				        if (indProp == string::npos)
				        {
				            srefProps.push_back(NULL);
				        }
				        else
				        {
				            size_t indPropSpace = fileLine.find(" ", indProp);
				            unsigned long propAttr = std::stoul(fileLine.substr(indProp + 3, indPropSpace - indProp));
                            vector<string>* property = new vector<string>;
                            (*property).assign(propAttr, ""); // Make sure enough elements available for property to go to correct location
                            (*property)[propAttr - 1] = fileLine.substr(indPropSpace + 2, fileLine.find("}", indProp) - indPropSpace - 3); // Remove surrounding single quotation marks
				            srefProps.push_back(property);
				        }
				    }
					// Handle nodes
					else if (fileLine.compare(0, 2, "n{") == 0)
					{
						// Extract layer number
						size_t indLayer = fileLine.find(" ") - 1;
						nodeLayers.push_back(std::stoul(fileLine.substr(2, indLayer)));

						// Extract node types
						size_t indType = fileLine.find(" nt");
						if (indType == string::npos)
						{
							indType = fileLine.find(" xt"); // Specified as box type instead
							if (indType == string::npos)
							{
								nodeTypes.push_back(0);
							}
							else
							{
								nodeTypes.push_back(std::stoul(fileLine.substr(indType + 3, fileLine.find(" ", indType + 1) - indType - 3)));
							}
						}
						else
						{
							nodeTypes.push_back(std::stoul(fileLine.substr(indType + 3, fileLine.find(" ", indType + 1) - indType - 3)));
						}

						// Split string at spaces between successive coordinates
						size_t indCoordStart = fileLine.find("xy(") + 3;
						size_t indCoordStop = fileLine.find(")") - 1;
						vector<double>* coord = new vector<double>;
						size_t coordBegin = indCoordStart;
						while (coordBegin <= indCoordStop)
						{
							size_t coordEnd = fileLine.find(" ", coordBegin);
							if (coordEnd == string::npos)
							{
								coordEnd = indCoordStop; // Final coordinate
							}
							else
							{
								coordEnd -= 1; // Remove trailing whitespace
							}
							(*coord).push_back(std::stod(fileLine.substr(coordBegin, coordEnd)));
							coordBegin = coordEnd + 2;
						}
						nodesAll.push_back(coord);

						// Search for node properties (might only handle one...)
						size_t indProp = fileLine.find("pr{");
						if (indProp == string::npos)
						{
							nodeProps.push_back(NULL);
						}
						else
						{
							size_t indPropSpace = fileLine.find(" ", indProp);
							unsigned long propAttr = std::stoul(fileLine.substr(indProp + 3, indPropSpace - indProp));
							vector<string>* property = new vector<string>;
							(*property).assign(propAttr, ""); // Make sure enough elements available for property to go to correct location
							(*property)[propAttr - 1] = fileLine.substr(indPropSpace + 2, fileLine.find("}", indProp) - indPropSpace - 3); // Remove surrounding single quotation marks
							nodeProps.push_back(property);
						}
					}
					// Handle text boxes
					else if (fileLine.compare(0, 2, "t{") == 0)
					{
						// Extract layer number
						size_t indLayer = fileLine.find(" ") - 1;
						textLayers.push_back(std::stoul(fileLine.substr(2, indLayer)));

				        // Extract text box widths
				        size_t indWidth = fileLine.find(" w");
				        if (indWidth == string::npos)
				        {
				            textWidths.push_back(0);
				        }
				        else
				        {
				            textWidths.push_back(std::stod(fileLine.substr(indWidth + 2, fileLine.find(" ", indWidth + 1) - indWidth - 2)));
				        }

                        // Extract text justification
                        size_t indType = fileLine.find(" tl ");
                        if (indType == string::npos)
                        {
                            indType = fileLine.find(" tc ");
                            if (indType == string::npos)
                            {
                                indType = fileLine.find(" tr ");
                                if (indType == string::npos)
                                {
                                    indType = fileLine.find(" ml ");
                                    if (indType == string::npos)
                                    {
                                        indType = fileLine.find(" mc ");
                                        if (indType == string::npos)
                                        {
indType = fileLine.find(" mr ");
if (indType == string::npos)
{
    indType = fileLine.find(" bl ");
    if (indType == string::npos)
    {
        indType = fileLine.find(" bc ");
        if (indType == string::npos)
        {
            indType = fileLine.find(" br ");
            if (indType == string::npos)
            {
                textJusts.push_back({2, 2}); // Bottom right
            }
            else
            {
                textJusts.push_back({1, 1}); // Default: middle center
            }
        }
        else
        {
            textJusts.push_back({2, 1}); // Bottom center
        }
    }
    else
    {
        textJusts.push_back({2, 0}); // Bottom left
    }
}
else
{
    textJusts.push_back({1, 2}); // Middle right
}
                                        }
                                        else
                                        {
                                            textJusts.push_back({1, 1}); // Middle center
                                        }
                                    }
                                    else
                                    {
                                        textJusts.push_back({1, 0}); // Middle left
                                    }
                                }
                                else
                                {
                                    textJusts.push_back({0, 2}); // Top right
                                }
                            }
                            else
                            {
                                textJusts.push_back({0, 1}); // Top center
                            }
                        }
                        else
                        {
                            textJusts.push_back({0, 0}); // Top left
                        }

						// Split string at space between two coordinates
				        size_t indCoordStart = fileLine.find("xy(") + 3;
				        size_t indCoordSep = fileLine.find(" ", indCoordStart);
				        size_t indCoordStop = fileLine.find(")") - 1;
				        vector<double> coord;
				        coord.assign(2, 0); // Coordinate pair initialized to origin
				        coord[0] = std::stod(fileLine.substr(indCoordStart, indCoordSep - indCoordStart));
				        coord[1] = std::stod(fileLine.substr(indCoordSep + 1, indCoordStop - indCoordSep));
				        textsAll.push_back(coord);

                        // Extract text box string
                        size_t indTStrStart = fileLine.find("'", 3) + 1;
                        size_t indTStrEnd = fileLine.find("'", indTStrStart);
				        textStr.push_back(fileLine.substr(indTStrStart, indTStrEnd));
					}
				}

                // Create full structure of geometric cell
                cells.push_back(GeoCell(cellName, cellMetadata, boundLayers, boundsAll, boundProps, pathLayers, pathTypes, pathWidths, pathsAll, pathProps, srefNames, srefsAll, srefProps, nodeLayers, nodeTypes, nodesAll, nodeProps, textLayers, textWidths, textJusts, textsAll, textStr));
			}
			// Skip additional comments below comment block
            else if (fileLine[0] == '#')
            {
                continue;
            }
			// Update amount of file read
			this->fileSize = gdtFile.tellg();
			getline(gdtFile, fileLine);
		}

		// Close GDT file
	    gdtFile.close();
		return true;
	}
	else
	{
		// File could not be opened
	    return false;
	}
}

// Get file size
long Geometry::getSize() const
{
    return this->fileSize;
}

// Get number of cells
int Geometry::getNumCell() const
{
    return cells.size();
}

// Find index of cell by name
size_t Geometry::locateCell(string name) const
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
GeoCell Geometry::getCell(size_t indCell) const
{
    return cells[indCell];
}

// Return all cell names
vector<string> Geometry::findNames() const
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
vector<unsigned long> Geometry::findLayers() const
{
    vector<unsigned long> layers;
    //size_t allEntities = 0;
    for (size_t indi = 0; indi < this->getNumCell(); indi++)
    {
        GeoCell cell = this->getCell(indi);
        //allEntities += (cell.getNumBound() + cell.getNumPath());
        vector<unsigned long> bLayers = cell.getBoundLayers();
        vector<unsigned long> pLayers = cell.getPathLayers();
		vector<unsigned long> nLayers = cell.getNodeLayers();
		vector<unsigned long> tLayers = cell.getTextLayers();
        layers.insert(layers.end(), bLayers.begin(), bLayers.end());
        layers.insert(layers.end(), pLayers.begin(), pLayers.end());
		layers.insert(layers.end(), nLayers.begin(), nLayers.end());
		layers.insert(layers.end(), tLayers.begin(), tLayers.end());
    }
    sort(layers.begin(), layers.end());
    layers.erase(unique(layers.begin(), layers.end()), layers.end());
    return layers;
}

// Return via list
vector<vector<double>> Geometry::findVias(size_t indCell, vector<double> center, string viaName)
{
	GeoCell thisCell = this->getCell(indCell);
	size_t indVia = this->locateCell(viaName);
	vector<string> srefNames = thisCell.getSRefNames();
	vector<vector<double>> srefs = thisCell.getSRefs();
	vector<vector<double>> viaList;
	for (size_t indj = 0; indj < srefNames.size(); indj++) // Check each structure reference ...
	{
		size_t indNextCell = this->locateCell(srefNames[indj]);
		if (indNextCell == indVia) // Cell is a via ...
		{
			viaList.push_back(srefs[indj]);
		}
		else if (indNextCell < this->getNumCell()) // Non-via cell name was found ...
		{
			vector<double> thisSRef = srefs[indj];
			vector<vector<double>> newVias = this->findVias(indNextCell, { (center[0] + thisSRef[0]), (center[1] + thisSRef[1]) }, viaName); // Recursion step
			viaList.insert(viaList.end(), newVias.begin(), newVias.end());
		}
	}
	return viaList;
}

// Print the geometry
void Geometry::print() const
{
	int numCell = Geometry::getNumCell();

    cout << "Geometry:" << endl;
    cout << " File Name: " << this->fileName << endl;
    cout << " Metadata: " << this->metadata << endl;
	cout << " List of " << numCell << " cells:" << endl;
	for (size_t indi = 0; indi < numCell; indi++)
	{
		cout << "  " << indi + 1 << ". " << cells[indi].getCellName() << endl;
		cout << "   Metadata: " << cells[indi].getCellMeta() << endl;
		cout << "   Counts: " << cells[indi].getNumBound() << " boundaries, " << cells[indi].getNumPath() << " paths, " << cells[indi].getNumSRef() << " structures, " << cells[indi].getNumNode() << " nodes, and " << cells[indi].getNumText() << " text boxes" << endl;
	}
    cout << " File size: " << this->fileSize << endl;
    cout << "------" << endl;
}

// Write gnuplot file for layer of a cell
vector<int> Geometry::writeCellPlotter(ofstream& gnuFile, size_t indCell, vector<double> center, unsigned long layer, int bNum, int pNum, int nNum, int tNum)
{
    GeoCell thisCell = this->getCell(indCell);
    vector<unsigned long> bLayers = thisCell.getBoundLayers();
    vector<vector<double>*> bounds = thisCell.getBounds();
    vector<unsigned long> pLayers = thisCell.getPathLayers();
	vector<unsigned long> pTypes = thisCell.getPathTypes();
    vector<double> pWidths = thisCell.getPathWidths();
    vector<vector<double>*> paths = thisCell.getPaths();
    vector<string> srefNames = thisCell.getSRefNames();
    vector<vector<double>> srefs = thisCell.getSRefs();
	vector<unsigned long> nLayers = thisCell.getNodeLayers();
	vector<vector<double>*> nodes = thisCell.getNodes();
	vector<unsigned long> tLayers = thisCell.getTextLayers();
	vector<vector<unsigned long>> tJusts = thisCell.getTextJusts();
	vector<vector<double>> texts = thisCell.getTexts();
	vector<string> textStr = thisCell.getTextStr();
    for (size_t indj = 0; indj < bLayers.size(); indj++) // Check each boundary ...
    {
        if (bLayers[indj] == layer) // To draw boundaries on that layer
        {
            vector<double> thisBound = *(bounds[indj]);
            gnuFile << "set obj " << bNum << " polygon from ";
            gnuFile << (thisBound[0] + center[0]) << ", " << (thisBound[1] + center[1]) << ", 0";
            for (size_t indk = 2; indk < thisBound.size(); indk++)
            {
                gnuFile << " to " << (thisBound[indk++] + center[0]) << ", " << (thisBound[indk + 1] + center[1]) << ", 0";
            }
            gnuFile << " to " << (thisBound[0] + center[0]) << ", " << (thisBound[1] + center[1]) << ", 0" << endl; // Boundary polygon ends on the first point
            gnuFile << "set obj " << bNum << " fs solid 0.7" << endl;
            bNum++;
        }
    }
    for (size_t indj = 0; indj < pLayers.size(); indj++) // Check each path ...
    {
        if (pLayers[indj] == layer) // To draw paths on that layer
        {
            vector<double> thisPath = *(paths[indj]);
            for (size_t indk = 0; indk < thisPath.size() - 2; indk++)
            {
                //gnuFile << "set arrow nohead linewidth " << pWidths[indj] << " from " << (thisPath[indk++] + center[0]) << ", " << (thisPath[indk + 1] + center[1]) << ", 0" << " to " << (thisPath[indk + 2] + center[0]) << ", " << (thisPath[indk + 3] + center[1]) << ", 0" << endl;
                if (thisPath[indk] == thisPath[indk + 2]) // Mahattan move with same x-coordinate
                {
					gnuFile << "set obj " << pNum << " polygon from ";
                    gnuFile << (thisPath[indk] + center[0] - 0.5*pWidths[indj]) << ", " << (thisPath[indk + 1] + center[1]) << ", 0"; // First x diminished
                    gnuFile << " to " << (thisPath[indk + 2] + center[0] - 0.5*pWidths[indj]) << ", " << (thisPath[indk + 3] + center[1]) << ", 0"; // Second x diminished
                    gnuFile << " to " << (thisPath[indk + 2] + center[0] + 0.5*pWidths[indj]) << ", " << (thisPath[indk + 3] + center[1]) << ", 0"; // Second x augemented
                    gnuFile << " to " << (thisPath[indk] + center[0] + 0.5*pWidths[indj]) << ", " << (thisPath[indk + 1] + center[1]) << ", 0"; // First x augmented
                    gnuFile << " to " << (thisPath[indk] + center[0] - 0.5*pWidths[indj]) << ", " << (thisPath[indk + 1] + center[1]) << ", 0" << endl; // Path segment polygon ends on the first point
					gnuFile << "set obj " << pNum << " fs solid 0.9" << endl;
					pNum++;
					if ((indk != 0) || (pTypes[indj] == 2))
					{
						gnuFile << "set obj " << pNum << "rect center " << (thisPath[indk] + center[0]) << ", " << (thisPath[indk + 1] + center[1]) << ", 0 size " << pWidths[indj] << "," << pWidths[indj] << endl; // Square covering first point
						gnuFile << "set obj " << pNum << " fc black fs solid 0.9" << endl;
						pNum++;
					}
                    indk++;
                }
                else if (thisPath[indk + 1] == thisPath[indk + 3]) // Manhattan move with same y-coordinate
                {
					gnuFile << "set obj " << pNum << " polygon from ";
                    gnuFile << (thisPath[indk] + center[0]) << ", " << (thisPath[indk + 1] + center[1] - 0.5*pWidths[indj]) << ", 0"; // First y diminished
                    gnuFile << " to " << (thisPath[indk + 2] + center[0]) << ", " << (thisPath[indk + 3] + center[1] - 0.5*pWidths[indj]) << ", 0"; // Second y diminished
                    gnuFile << " to " << (thisPath[indk + 2] + center[0]) << ", " << (thisPath[indk + 3] + center[1] + 0.5*pWidths[indj]) << ", 0"; // Second y augemented
                    gnuFile << " to " << (thisPath[indk] + center[0]) << ", " << (thisPath[indk + 1] + center[1] + 0.5*pWidths[indj]) << ", 0"; // First y augmented
                    gnuFile << " to " << (thisPath[indk] + center[0]) << ", " << (thisPath[indk + 1] + center[1] - 0.5*pWidths[indj]) << ", 0" << endl; // Path segment polygon ends on the first point
					gnuFile << "set obj " << pNum << " fs solid 0.9" << endl;
					pNum++;
					if ((indk != 0) || (pTypes[indj] == 2))
					{
						gnuFile << "set obj " << pNum << "rect center " << (thisPath[indk] + center[0]) << ", " << (thisPath[indk + 1] + center[1]) << ", 0 size " << pWidths[indj] << "," << pWidths[indj] << endl; // Square covering first point in segment
						gnuFile << "set obj " << pNum << " fc black fs solid 0.9" << endl;
						pNum++;
					}
                    indk++;
                }
                else // Non-Manhattan move
                {
                    gnuFile << "set arrow nohead linewidth " << pWidths[indj] << " from " << (thisPath[indk++] + center[0]) << ", " << (thisPath[indk + 1] + center[1]) << ", 0" << " to " << (thisPath[indk + 2] + center[0]) << ", " << (thisPath[indk + 3] + center[1]) << ", 0" << endl;
                }
            }
			if (pTypes[indj] == 2) // Type 2 overshoot end style
			{
				gnuFile << "set obj " << pNum << "rect center " << (thisPath[thisPath.size() - 2] + center[0]) << ", " << (thisPath[thisPath.size() - 1] + center[1]) << ", 0 size " << pWidths[indj] << "," << pWidths[indj] << endl; // Square covering last point (very first point handled in loop)
				gnuFile << "set obj " << pNum << " fc black fs solid 0.9" << endl;
				pNum++;
			}
			else if (pTypes[indj] == 1) // Type 1 round end style
			{
				gnuFile << "set obj " << pNum << "circ at " << (thisPath[0] + center[0]) << ", " << (thisPath[1] + center[1]) << ", 0 size " << (0.5*pWidths[indj]) << endl; // Circle covering very first point
				gnuFile << "set obj " << pNum << " fc black fs solid 0.9" << endl;
				pNum++;
				gnuFile << "set obj " << pNum << "circ at " << (thisPath[thisPath.size() - 2] + center[0]) << ", " << (thisPath[thisPath.size() - 1] + center[1]) << ", 0 size " << (0.5*pWidths[indj]) << endl; // Circle covering last point
				gnuFile << "set obj " << pNum << " fc black fs solid 0.9" << endl;
				pNum++;
			}
        }
    }
	for (size_t indj = 0; indj < nLayers.size(); indj++) // Check each node ...
	{
		if (nLayers[indj] == layer) // To draw nodes on that layer
		{
			vector<double> thisNode = *(nodes[indj]);
			gnuFile << "set obj " << nNum << " polygon from ";
			gnuFile << (thisNode[0] + center[0]) << ", " << (thisNode[1] + center[1]) << ", 0";
			for (size_t indk = 2; indk < thisNode.size(); indk++)
			{
				gnuFile << " to " << (thisNode[indk++] + center[0]) << ", " << (thisNode[indk + 1] + center[1]) << ", 0";
			}
			gnuFile << endl; // Node polygons inherently end on the first point
			gnuFile << "set obj " << nNum << " fs transparent solid 0.3 dt \"..\"" << endl;
			nNum++;
		}
	}
	for (size_t indj = 0; indj < tLayers.size(); indj++) // Check each text box ...
	{
		if (tLayers[indj] == layer) // To write text boxes on that layer
		{
			vector<double> thisText = texts[indj];
			vector<unsigned long> thisJust = tJusts[indj];
			gnuFile << "set label " << tNum << " norotate front at ";
			gnuFile << (thisText[0] + center[0]) << ", " << (thisText[1] + center[1]) << ", 0 ";
			gnuFile << "\"" << textStr[indj] << "\" ";
			if (thisJust[1] == 0) // gnuplot only supports horizontal alignment
			{
			    gnuFile << " left";
			}
			else if (thisJust[1] == 2)
			{
			    gnuFile << " right";
			}
			else
			{
			    gnuFile << " center"; // Center by default
			}
			gnuFile << endl;
			tNum++;
		}
	}
    int newBNum = bNum;
	int newPNum = pNum;
	int newNNum = nNum;
	int newTNum = tNum;
	vector<int> newObjNums = {bNum, pNum, nNum, tNum};
    for (size_t indj = 0; indj < srefNames.size(); indj++) // Check each structure reference ...
    {
        size_t indNextCell = this->locateCell(srefNames[indj]);
        if (indNextCell < this->getNumCell()) // Cell name was found ...
        {
            vector<double> thisSRef = srefs[indj];
            newObjNums = this->writeCellPlotter(gnuFile, indNextCell, {(center[0] + thisSRef[0]), (center[1] + thisSRef[1])}, layer, newBNum, newPNum, newNNum, newTNum); // Recursion step
			newBNum = newObjNums[0];
			newPNum = newObjNums[1];
			newNNum = newObjNums[2];
			newTNum = newObjNums[3];
        }
    }
    return newObjNums;
}

// Write gnuplot environment
// Return false if unable to write a gnuplot file
bool Geometry::writePlotEnv(std::string filePrefix, std::string masterCell, unsigned long pixels, vector<double> layerBounds)
{
	// Attempt to open file
	string gnuName = "./" + filePrefix + "_recreate.gnu";
	ofstream gnuFile(gnuName.c_str());
	if (gnuFile.is_open())
	{
		// Set gnuplot environment
		gnuFile << "set terminal pngcairo transparent enhanced font \"arial, 10\" fontscale 1.0 size " << pixels << ", " << pixels << endl;

		// Plot cell layer by layer
		size_t indCell = this->locateCell(masterCell); // The master design
		vector<unsigned long> layers = this->findLayers();
		int bNum = 1;
		int pNum = 10001;
		int nNum = 20001;
		int tNum = 30001;
		for (size_t indi = 0; indi < layers.size(); indi++) // On each layer...
		{
			gnuFile << "# Layer " << layers[indi] << endl;
			gnuFile << "reset" << endl;
			gnuFile << "set output '" << filePrefix << "_recreate/" << filePrefix << "_recreate" << layers[indi] << ".png'" << endl;
			gnuFile << "set key off" << endl;
			gnuFile << "set grid" << endl;
			gnuFile << "set xrange [" << layerBounds[0] << ":" << layerBounds[1] << "] noreverse nowriteback" << endl;
			gnuFile << "set yrange [" << layerBounds[2] << ":" << layerBounds[3] << "] noreverse nowriteback" << endl;
			//gnuFile << "set size square 1,1" << endl;
			gnuFile << "set view equal xy" << endl;
			vector<int> objNums = this->writeCellPlotter(gnuFile, indCell, { 0, 0 }, layers[indi], bNum, pNum, nNum, tNum); // Call recursive function for first time
			bNum = objNums[0];
			pNum = objNums[1];
			nNum = objNums[2];
			tNum = objNums[3];
			gnuFile << "plot -1000" << endl;
		}

		// Close file
		gnuFile.close();
		return true;
	}
	else
	{
		// File could not be opened
		return false;
	}
}

// Destructor
Geometry::~Geometry()
{
	/*for (size_t indi = 0; indi < this->getNumCell(); indi++)
	{
	    cells[indi].~GeoCell();
	}*/
}
