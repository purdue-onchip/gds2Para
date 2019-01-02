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
#include <cstring>
#include <string>
#include <algorithm>
#include <stack>
#include <tr1/unordered_map>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
#include <limbo/parsers/gdsii/stream/GdsWriter.h>
#include "fdtd.h"
using std::cerr;
using std::cout;
using std::endl;

class boundary
{
private:
    vector<double> bounds;         // Coordinates of boundary
    int layer;                     // Layer of boundary
    vector<std::string> props;     // Properties of boundary
public:
    // Default constructor
    boundary()
    {
        this->bounds = {};
        this->layer = 0;
        this->props = {};
    }

    // Parametrized constructor
    boundary(vector<double> bounds, int layer, vector<std::string> props)
    {
        this->bounds = bounds;
        this->layer = layer;
        this->props = props;
    }

    // Get boundary points in cell coordinates
    vector<double> getBounds() const
    {
        return this->bounds;
    }

    // Get boundary layer number
    int getLayer() const
    {
        return this->layer;
    }

    // Get boundary properties
    vector<std::string> getProps() const
    {
        return this->props;
    }

    // Return number of boundary points
    size_t getNBoundPt() const
    {
        size_t nBoundCoord = (this->bounds).size(); // Number of coordinates stored
        if ((nBoundCoord % 2) == 0)
        {
            return nBoundCoord / 2;
        }
        else
        {
            cerr << "Should be an even number of coordinates stored for boundary points" << endl;
            return 0;
        }
    }

    // Reorder points with lower-right first, going counterclockwise
    void reorder()
    {
        // Remove last point if equal to first point
        size_t nBoundPt = this->getNBoundPt(); // Number of boundary points (coordinate pairs)
        if (((this->bounds)[2 * 0] == (this->bounds)[2 * (nBoundPt - 1)]) && ((this->bounds)[2 * 0 + 1] == (this->bounds)[2 * (nBoundPt - 1) + 1]))
        {
            (this->bounds).pop_back(); // Remove repeated y-coordinate
            (this->bounds).pop_back(); // Remove repeated x-coordinate
        }

        // Find lower-right most point
        size_t ptLR = 0; // Index of lower-right point
        for (size_t indi = 1; indi < nBoundPt; indi++)
        {
            if ((this->bounds)[2 * indi] >= (this->bounds)[2 * ptLR]) // x-coordinate right of rightmost
            {
                if ((this->bounds)[2 * indi + 1] <= (this->bounds)[2 * ptLR + 1]) // y-coordinate below lowest
                {
                    ptLR = indi;
                }
            }
        }

        // Circularly rotate the vector to lower-right most point first
        rotate((this->bounds).begin(), (this->bounds).begin() + 2 * ptLR, (this->bounds).end());

        // Correct direction if clockwise
        if (((this->bounds)[2 * 0] > (this->bounds)[2 * 1]) && ((this->bounds)[2 * 0 + 1] < (this->bounds)[2 * 1 + 1]))
        {
            vector<double> revBound;
            revBound.push_back((this->bounds)[0]); // Same first x-coordinate
            revBound.push_back((this->bounds)[1]); // Same first y-coordinate
            for (size_t indi = nBoundPt - 1; indi > 0; indi++)
            {
                revBound.push_back((this->bounds)[2 * indi]); // x-coordinate
                revBound.push_back((this->bounds)[2 * indi + 1]); // y-coordinate
            }
            this->bounds = revBound; // New value for bounds
        }

        // Ensure last point is equal to first in lower-right
        (this->bounds).push_back((this->bounds)[0]); // Repeat the first x-coordinate
        (this->bounds).push_back((this->bounds)[1]); // Repeat the first y-coordinate
    }

    // Destructor
    ~boundary()
    {
        //this->bounds = {};
        //(this->bounds).clear();
        this->layer = 0;
        //this->props = {};
    }
};

class path
{
private:
    vector<double> paths;          // Coordinates of path
    int layer;                     // Layer of path
    vector<std::string> props;     // Properties of path
    int type;                      // Endcap type of path
    double width;                  // Width of path (m)
public:
    // Default constructor
    path()
    {
        this->paths = {};
        this->layer = 0;
        this->props = {};
        this->type = 2;
        this->width = 0.;
    }

    // Parametrized constructor
    path(vector<double> paths, int layer, vector<std::string> props, int type, double width)
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
    vector<double> getPaths() const
    {
        return this->paths;
    }

    // Get path layer number
    int getLayer() const
    {
        return this->layer;
    }

    // Get path properties
    vector<std::string> getProps() const
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

    // Return number of path points
    size_t getNPathPt() const
    {
        size_t nPathCoord = (this->paths).size(); // Number of coordinates stored
        if ((nPathCoord % 2) == 0)
        {
            return nPathCoord / 2;
        }
        else
        {
            cerr << "Should be an even number of coordinates stored for path vertex points" << endl;
            return 0;
        }
    }

    // Destructor
    ~path()
    {
        //this->paths = {};
        this->layer = 0;
        this->type = 2;
        this->width = 0.;
    }
};

class node
{
private:
    vector<double> nodes;          // Coordinates of electrical node
    int layer;                     // Layer of electrical node
    vector<std::string> props;     // Properties of electrical node
    int type;                      // Type of electrical node
public:
    // Default constructor
    node()
    {
        this->nodes = {};
        this->layer = 0;
        this->props = {};
        this->type = 0;
    }

    // Parametrized constructor
    node(vector<double> nodes, int layer, vector<std::string> props, int type)
    {
        this->nodes = nodes;
        this->layer = layer;
        this->props = props;
        this->type = type;
    }

    // Get node vertices in cell coordinates
    vector<double> getNodes() const
    {
        return this->nodes;
    }

    // Get node layer number
    int getLayer() const
    {
        return this->layer;
    }

    // Get node properties
    vector<std::string> getProps() const
    {
        return this->props;
    }

    // Get node type
    int getType() const
    {
        return this->type;
    }

    // Return number of node points
    size_t getNNodePt() const
    {
        size_t nNodeCoord = (this->nodes).size(); // Number of coordinates stored
        if ((nNodeCoord % 2) == 0)
        {
            return nNodeCoord / 2;
        }
        else
        {
            cerr << "Should be an even number of coordinates stored for node points" << endl;
            return 0;
        }
    }

    // Reorder points with lower-right first, going counterclockwise
    void reorder()
    {
        // Remove last point if equal to first point
        size_t nNodePt = this->getNNodePt(); // Number of node points (coordinate pairs)
        if (((this->nodes)[2 * 0] == (this->nodes)[2 * (nNodePt - 1)]) && ((this->nodes)[2 * 0 + 1] == (this->nodes)[2 * (nNodePt - 1) + 1]))
        {
            (this->nodes).pop_back(); // Remove repeated y-coordinate
            (this->nodes).pop_back(); // Remove repeated x-coordinate
        }

        // Find lower-right most point
        size_t ptLR = 0; // Index of lower-right point
        for (size_t indi = 1; indi < nNodePt; indi++)
        {
            if ((this->nodes)[2 * indi] >= (this->nodes)[2 * ptLR]) // x-coordinate right of rightmost
            {
                if ((this->nodes)[2 * indi + 1] <= (this->nodes)[2 * ptLR + 1]) // y-coordinate below lowest
                {
                    ptLR = indi;
                }
            }
        }

        // Circularly rotate the vector to lower-right most point first
        rotate((this->nodes).begin(), (this->nodes).begin() + 2 * ptLR, (this->nodes).end());

        // Correct direction if clockwise
        if (((this->nodes)[2 * 0] > (this->nodes)[2 * 1]) && ((this->nodes)[2 * 0 + 1] < (this->nodes)[2 * 1 + 1]))
        {
            vector<double> revNode;
            revNode.push_back((this->nodes)[0]); // Same first x-coordinate
            revNode.push_back((this->nodes)[1]); // Same first y-coordinate
            for (size_t indi = nNodePt - 1; indi > 0; indi++)
            {
                revNode.push_back((this->nodes)[2 * indi]); // x-coordinate
                revNode.push_back((this->nodes)[2 * indi + 1]); // y-coordinate
            }
            this->nodes = revNode; // New value for bounds
        }

        // Ensure last point is equal to first in lower-right
        (this->nodes).push_back((this->nodes)[0]); // Repeat the first x-coordinate
        (this->nodes).push_back((this->nodes)[1]); // Repeat the first y-coordinate
    }

    // Destructor
    ~node()
    {
        //this->nodes = {};
        this->layer = 0;
        this->type = 0;
    }
};

class box
{
private:
    vector<double> boxes;          // Coordinates of box outline
    int layer;                     // Layer of box outline
    vector<std::string> props;     // Properties of box outline
    int type;                      // Type of box outline
public:
    // Default constructor
    box()
    {
        this->boxes = {};
        this->layer = 0;
        this->props = {};
        this->type = 0;
    }

    // Parametrized constructor
    box(vector<double> boxes, int layer, vector<std::string> props, int type)
    {
        this->boxes = boxes;
        this->layer = layer;
        this->props = props;
        this->type = type;
    }

    // Get box outline vertices in cell coordinates
    vector<double> getBoxes() const
    {
        return this->boxes;
    }

    // Get box outline layer number
    int getLayer() const
    {
        return this->layer;
    }

    // Get box outline properties
    vector<std::string> getProps() const
    {
        return this->props;
    }

    // Get box outline type
    int getType() const
    {
        return this->type;
    }

    // Return number of box outline points
    size_t getNBoxPt() const
    {
        size_t nBoxCoord = (this->boxes).size(); // Number of coordinates stored
        if ((nBoxCoord % 2) == 0)
        {
            return nBoxCoord / 2;
        }
        else
        {
            cerr << "Should be an even number of coordinates stored for box outline points" << endl;
            return 0;
        }
    }

    // Reorder points with lower-right first, going counterclockwise
    void reorder()
    {
        // Remove last point if equal to first point
        size_t nBoxPt = this->getNBoxPt(); // Number of node points (coordinate pairs)
        if (((this->boxes)[2 * 0] == (this->boxes)[2 * (nBoxPt - 1)]) && ((this->boxes)[2 * 0 + 1] == (this->boxes)[2 * (nBoxPt - 1) + 1]))
        {
            (this->boxes).pop_back(); // Remove repeated y-coordinate
            (this->boxes).pop_back(); // Remove repeated x-coordinate
        }

        // Find lower-right most point
        size_t ptLR = 0; // Index of lower-right point
        for (size_t indi = 1; indi < nBoxPt; indi++)
        {
            if ((this->boxes)[2 * indi] >= (this->boxes)[2 * ptLR]) // x-coordinate right of rightmost
            {
                if ((this->boxes)[2 * indi + 1] <= (this->boxes)[2 * ptLR + 1]) // y-coordinate below lowest
                {
                    ptLR = indi;
                }
            }
        }

        // Circularly rotate the vector to lower-right most point first
        rotate((this->boxes).begin(), (this->boxes).begin() + 2 * ptLR, (this->boxes).end());

        // Correct direction if clockwise
        if (((this->boxes)[2 * 0] > (this->boxes)[2 * 1]) && ((this->boxes)[2 * 0 + 1] < (this->boxes)[2 * 1 + 1]))
        {
            vector<double> revBox;
            revBox.push_back((this->boxes)[0]); // Same first x-coordinate
            revBox.push_back((this->boxes)[1]); // Same first y-coordinate
            for (size_t indi = nBoxPt - 1; indi > 0; indi++)
            {
                revBox.push_back((this->boxes)[2 * indi]); // x-coordinate
                revBox.push_back((this->boxes)[2 * indi + 1]); // y-coordinate
            }
            this->boxes = revBox; // New value for bounds
        }

        // Ensure last point is equal to first in lower-right
        (this->boxes).push_back((this->boxes)[0]); // Repeat the first x-coordinate
        (this->boxes).push_back((this->boxes)[1]); // Repeat the first y-coordinate
    }

    // Destructor
    ~box()
    {
        //this->boxes = {};
        this->layer = 0;
        this->type = 0;
    }
};

class textbox
{
private:
    vector<double> texts;          // Coordinate pair of text box center
    int layer;                     // Layer of text box
    vector<std::string> props;     // Properties of text box
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
        this->props = {};
        this->type = 0;
        this->fontID = 0;
        this->justs = justs;
        this->width = 0.;
        this->textStr = "";
    }

    // Parametrized constructor
    textbox(vector<double> texts, int layer, vector<std::string> props, int type, int fontID, vector<int> justs, double width, std::string textStr)
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
    vector<std::string> getProps() const
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
    vector<std::string> props;     // Properties of structure reference
    std::string srefName;          // Name of structure reference
public:
    // Default constructor
    sref()
    {
        vector<double> srefs = { 0., 0. };
        this->srefs = srefs;
        this->props = {};
        this->srefName = "";
    }

    // Parametrized constructor
    sref(vector<double> srefs, vector<std::string> props, std::string srefName)
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
    vector<std::string> getProps() const
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
            if (((this->boundaries)[indi]).getProps().size() != 0) // Check if property vector has nonzero length
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->boundaries)[indi]).getLayer() << " with " << (((this->boundaries)[indi]).getProps())[0] << " as Property " << 1 << endl;
            }
            else
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->boundaries)[indi]).getLayer() << endl;
            }
            char point[128];
            string strPoints;
            vector<double> boundCoord = ((this->boundaries)[indi]).getBounds();
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
            if (((this->paths)[indi]).getProps().size() != 0) // Check if property vector has nonzero length
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->paths)[indi]).getLayer() << " with path type #" << ((this->paths)[indi]).getType() << " and width " << ((this->paths)[indi]).getWidth() << " with " << (((this->paths)[indi]).getProps())[0] << " as Property " << 1 << endl;
            }
            else
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->paths)[indi]).getLayer() << " with path type #" << ((this->paths)[indi]).getType() << " and width " << ((this->paths)[indi]).getWidth() << endl;
            }
            char point[128];
            string strPoints;
            vector<double> pathCoord = ((this->paths)[indi]).getPaths();
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
            if (((this->nodes)[indi]).getProps().size() != 0) // Check if property vector has nonzero length
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->nodes)[indi]).getLayer() << " with node type #" << ((this->nodes)[indi]).getType() << " with " << (((this->nodes)[indi]).getProps())[0] << " as Property " << 1 << endl;
            }
            else
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->nodes)[indi]).getLayer() << " with node type #" << ((this->nodes)[indi]).getType() << endl;
            }
            char point[128];
            string strPoints;
            vector<double> nodeCoord = ((this->nodes)[indi]).getNodes();
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
            if (((this->boxes)[indi]).getProps().size() != 0) // Check if property vector has nonzero length
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->boxes)[indi]).getLayer() << " with box outline type #" << ((this->boxes)[indi]).getType() << " with " << (((this->boxes)[indi]).getProps())[0] << " as Property " << 1 << endl;
            }
            else
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->boxes)[indi]).getLayer() << " with box outline type #" << ((this->boxes)[indi]).getType() << endl;
            }
            char point[128];
            string strPoints;
            vector<double> boxCoord = ((this->boxes)[indi]).getBoxes();
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
            if (((this->textboxes)[indi]).getProps().size() != 0) // Check if property vector has nonzero length
            {
                cout << "   #" << indi + 1 << " on Layer " << ((this->textboxes)[indi]).getLayer() << " with text box type #" << ((this->textboxes)[indi]).getType() << " and message \"" << ((this->textboxes)[indi]).getTextStr() << "\" with " << (((this->textboxes)[indi]).getProps())[0] << " as Property " << 1 << endl;
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
            if (((this->sreferences)[indi]).getProps().size() != 0) // Check if property vector has nonzero length
            {
                cout << "   #" << indi + 1 << " named " << ((this->sreferences)[indi]).getSRefName() << " with " << (((this->sreferences)[indi]).getProps())[0] << " as Property " << 1 << endl;
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

    // Terse alternative to print the geometric cell
    void printAlt() const
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
        cout << "  List of " << numBound << " boundaries [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numBound; indi++) // Handle each boundary
        {
            cout << "   " << ((this->boundaries)[indi]).getLayer() << " ";
            vector<double> boundCoord = ((this->boundaries)[indi]).getBounds();
            for (size_t indj = 0; indj < boundCoord.size(); indj++)
            {
                cout << boundCoord[indj++] << " " << boundCoord[indj + 1] << " ";
            }
            cout << ((this->boundaries)[indi]).getNBoundPt() << endl;
        }
        cout << "  List of " << numPath << " paths [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numPath; indi++) // Handle each path
        {
            cout << "   " << ((this->paths)[indi]).getLayer() << " ";
            vector<double> pathCoord = ((this->paths)[indi]).getPaths();
            for (size_t indj = 0; indj < pathCoord.size(); indj++)
            {
                cout << pathCoord[indj++] << " " << pathCoord[indj + 1] << " ";
            }
            cout << ((this->paths)[indi]).getNPathPt() << endl;
        }
        cout << "  List of " << numNode << " nodes [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numNode; indi++) // Handle each node
        {
            cout << "   " << ((this->nodes)[indi]).getLayer() << " ";
            vector<double> nodeCoord = ((this->nodes)[indi]).getNodes();
            for (size_t indj = 0; indj < nodeCoord.size(); indj++) // C string for each ordered pair
            {
                cout << nodeCoord[indj++] << " " << nodeCoord[indj + 1] << " ";
            }
            cout << ((this->nodes)[indi]).getNNodePt() << endl;
        }
        cout << "  List of " << numBox << " box outlines [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numBox; indi++) // Handle each box outline
        {
            cout << "   " << ((this->boxes)[indi]).getLayer() << " ";
            vector<double> boxCoord = ((this->boxes)[indi]).getBoxes();
            for (size_t indj = 0; indj < boxCoord.size(); indj++) // C string for each ordered pair (should be 5, sometimes 4)
            {
                cout << boxCoord[indj++] << " " << boxCoord[indj + 1] << " ";
            }
            cout << ((this->boxes)[indi]).getNBoxPt() << endl;
        }
        cout << "  List of " << numText << " text boxes [layer (x-coord, y-coord)]:" << endl;
        for (size_t indi = 0; indi < numText; indi++) // Handle each text box
        {
            cout << "   " << ((this->textboxes)[indi]).getLayer() << " ";
            char point[128];
            sprintf(point, "(%1.4g, %1.4g) ", (((this->textboxes)[indi]).getTexts())[0], (((this->textboxes)[indi]).getTexts())[1]);
            cout << point << endl; // Print ordered pair of center
        }
        cout << "  List of " << numSRef << " structure references [name (x-coord, y-coord)]:" << endl;
        for (size_t indi = 0; indi < numSRef; indi++) // Handle each structure reference
        {
            cout << "   " << ((this->sreferences)[indi]).getSRefName() << " ";
            char point[128];
            sprintf(point, "(%1.4g, %1.4g)", (((this->sreferences)[indi]).getSRefs())[0], (((this->sreferences)[indi]).getSRefs())[1]);
            cout << point << endl; // Print ordered pair of center
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
    int numCdtIn;                                // Number of conductor rows
    std::string strPoints;                       // Polygon output
    fdtdOneCondct a;                             // Template for the conductor vector
public:
    /// @brief constructor
    AsciiDataBase()
    {
        std::string fileName, version, dateMod, dateAccess, libName;
        double databaseUserUnits = 1.;
        double databaseUnits = 1.;
        char element;
        vector<GeoCell> cells;

        this->fileName = fileName;
        this->version = version;
        this->dateMod = dateMod;
        this->dateAccess = dateAccess;
        this->libName = libName;
        this->dbUserUnits = databaseUserUnits;
        this->dbUnits = databaseUnits;
        this->numCell = 0;
        this->element = element;
        this->numProp = 0;
        this->cells = cells;
        this->numCdtIn = 0;
        this->strPoints = "";
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
    size_t getNumProp() const
    {
        return this->numProp;
    }

    // Get number of conductor rows
    int getNumCdtIn() const
    {
        return this->numCdtIn;
    }

    // Get polygon output string
    std::string getPolygon()
    {
        return this->strPoints;
    }

    // Set file name
    void setFileName(std::string fileName)
    {
        this->fileName = fileName;
    }

    // Set date of last modification
    void setDateMod(std::string dateMod)
    {
        this->dateMod = dateMod;
    }

    // Set date of last access
    void setDateAccess(std::string dateAccess)
    {
        this->dateAccess = dateAccess;
    }

    // Set library name
    void setLibName(std::string libName)
    {
        this->libName = libName;
    }

    // Set database user units
    void setdbUserUnits(double dbUserUnits)
    {
        this->dbUserUnits = dbUserUnits;
    }

    // Set database units in SI
    void setdbUnits(double dbUnits)
    {
        this->dbUnits = dbUnits;
    }

    // Set number of conductor rows
    void setNumCdtIn(int numCdtIn)
    {
        this->numCdtIn = numCdtIn;
    }

    // Set the string for the polygon
    void setPolygon(std::string strPoints)
    {
        this->strPoints = strPoints;
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

    // Append a cell
    void appendCell(GeoCell extraCell)
    {
        (this->cells).push_back(extraCell);
        (this->numCell)++;
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

    // Print all the conductor information (fdtdMesh style)
    void printall(std::string name, double xo, double yo, fdtdMesh *sys)
    {
        int si, condj;    // Size of the vector sys->conductorIn

        // Open output file
        /*std::ofstream outfile;*/
        size_t indExtension = this->getFileName().find(".", 1);
        /*std::string polyFileName = this->getFileName().substr(0, indExtension) + "_polygon.txt";
        outfile.open(polyFileName, std::ofstream::out | std::ofstream::app);
        if (!outfile.is_open()) // Failed to open the polygon file to write
        {
            return;
        }*/

        // Get information about this cell in ASCII database
        const GeoCell cell = this->cells[this->locateCell(name)];
        int numBound = cell.getNumBound();
        int numPath = cell.getNumPath();
        int numNode = cell.getNumNode();
        int numBox = cell.getNumBox();
        int numText = cell.getNumText();
        int numSRef = cell.getNumSRef();

        // Print cell information
        cout << "  List of " << numBound << " boundaries:" << endl;
        for (size_t indi = 0; indi < numBound; indi++) // Handle each boundary
        {
            (this->numCdtIn)++;

            vector<double> boundCoord = ((cell.boundaries)[indi]).getBounds();
            sys->conductorIn.push_back(this->a);
            si = sys->conductorIn.size() - 1;
            sys->conductorIn[si].numVert = boundCoord.size() / 2 - 1;
            sys->conductorIn[si].xmax = DOUBLEMIN;
            sys->conductorIn[si].xmin = DOUBLEMAX;
            sys->conductorIn[si].ymax = DOUBLEMIN;
            sys->conductorIn[si].ymin = DOUBLEMAX;
            sys->conductorIn[si].x = (double*)calloc(sys->conductorIn[si].numVert, sizeof(double));
            sys->conductorIn[si].y = (double*)calloc(sys->conductorIn[si].numVert, sizeof(double));
            sys->conductorIn[si].layer = ((cell.boundaries)[indi]).getLayer();
            condj = 0;
            for (size_t indj = 0; indj < boundCoord.size() - 2; indj++) // C string for each ordered pair, -2 because the last point is the starting point
            {
                sys->conductorIn[si].x[condj] = boundCoord[indj] + xo;
                indj++;
                sys->conductorIn[si].y[condj] = boundCoord[indj] + yo;
                if (sys->conductorIn[si].x[condj] > sys->conductorIn[si].xmax){
                    sys->conductorIn[si].xmax = sys->conductorIn[si].x[condj];
                }
                if (sys->conductorIn[si].x[condj] < sys->conductorIn[si].xmin){
                    sys->conductorIn[si].xmin = sys->conductorIn[si].x[condj];
                }
                if (sys->conductorIn[si].y[condj] > sys->conductorIn[si].ymax){
                    sys->conductorIn[si].ymax = sys->conductorIn[si].y[condj];
                }
                if (sys->conductorIn[si].y[condj] < sys->conductorIn[si].ymin){
                    sys->conductorIn[si].ymin = sys->conductorIn[si].y[condj];
                }

                condj++;
            }
        }
        cout << "  List of " << numPath << " paths:" << endl;
        for (size_t indi = 0; indi < numPath; indi++) // Handle each path
        {
            vector<double> pathCoord = ((cell.paths)[indi]).getPaths();
            double width = ((cell.paths)[indi]).getWidth();
            sys->conductorIn.push_back(this->a);
            si = sys->conductorIn.size() - 1;
            sys->conductorIn[si].numVert = 4;
            sys->conductorIn[si].xmax = DOUBLEMIN;
            sys->conductorIn[si].xmin = DOUBLEMAX;
            sys->conductorIn[si].ymax = DOUBLEMIN;
            sys->conductorIn[si].ymin = DOUBLEMAX;
            sys->conductorIn[si].x = (double*)calloc(sys->conductorIn[si].numVert, sizeof(double));
            sys->conductorIn[si].y = (double*)calloc(sys->conductorIn[si].numVert, sizeof(double));
            sys->conductorIn[si].layer = ((cell.paths)[indi]).getLayer();
            condj = 0;
            for (size_t indj = 0; indj < pathCoord.size()-2; indj++) // C string for each ordered pair
            {
                if (pathCoord[indj] == pathCoord[indj + 2]) // along y axis
                {
                    if (pathCoord[indj + 1] > pathCoord[indj + 3]) // first point is on top of the second point
                    {
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        sys->conductorIn[si].xmin = pathCoord[indj] - width / 2 + xo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        sys->conductorIn[si].xmax = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].ymax = pathCoord[indj + 1] + width / 2 + yo;
                        indj++;
                        indj++;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        sys->conductorIn[si].ymin = pathCoord[indj + 1] - width / 2 + yo;
                        condj++;
                        this->numCdtIn++;
                        indj--;
                    }
                    else // second point is on the top of the first point
                    {
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        sys->conductorIn[si].xmin = pathCoord[indj] - width / 2 + xo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        sys->conductorIn[si].xmax = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].ymin = pathCoord[indj + 1] - width / 2 + yo;
                        indj++;
                        indj++;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        sys->conductorIn[si].ymax = pathCoord[indj + 1] + width / 2 + yo;
                        condj++;
                        this->numCdtIn++;
                        indj--;
                    }
                }
                else // along x axis
                {
                    if (pathCoord[indj] > pathCoord[indj + 2]) // first point is on the right of the second point
                    {
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        sys->conductorIn[si].xmax = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].ymax = pathCoord[indj + 1] + width / 2 + yo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        sys->conductorIn[si].ymin = pathCoord[indj + 1] - width / 2 + yo;
                        indj++;
                        indj++;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        sys->conductorIn[si].xmin = pathCoord[indj] - width / 2 + xo;
                        condj++;
                        this->numCdtIn++;
                        indj--;
                    }
                    else // second point is on the right of the first point
                    {
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        sys->conductorIn[si].xmin = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].ymax = pathCoord[indj + 1] + width / 2 + yo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] - width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        sys->conductorIn[si].ymin = pathCoord[indj + 1] - width / 2 + yo;
                        indj++;
                        indj++;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] - width / 2 + yo;
                        condj++;
                        sys->conductorIn[si].x[condj] = pathCoord[indj] + width / 2 + xo;
                        sys->conductorIn[si].y[condj] = pathCoord[indj + 1] + width / 2 + yo;
                        sys->conductorIn[si].xmax = pathCoord[indj] + width / 2 + xo;
                        condj++;
                        this->numCdtIn++;
                        indj--;
                    }
                }
            }

        }
        cout << "  List of " << numNode << " nodes:" << endl;
        cout << "  List of " << numBox << " box outlines:" << endl;
        for (size_t indi = 0; indi < numBox; indi++) // Handle each box outline
        {
            this->numCdtIn++;

            sys->conductorIn.push_back(this->a);
            vector<double> boxCoord = ((cell.boxes)[indi]).getBoxes();
            si = sys->conductorIn.size() - 1;
            sys->conductorIn[si].numVert = ((cell.boxes)[indi]).getNBoxPt() - 1;
            sys->conductorIn[si].layer = ((cell.boxes)[indi]).getLayer();
            sys->conductorIn[si].xmax = DOUBLEMIN;
            sys->conductorIn[si].xmin = DOUBLEMAX;
            sys->conductorIn[si].ymax = DOUBLEMIN;
            sys->conductorIn[si].ymin = DOUBLEMAX;
            sys->conductorIn[si].x = (double*)calloc(sys->conductorIn[si].numVert, sizeof(double));
            sys->conductorIn[si].y = (double*)calloc(sys->conductorIn[si].numVert, sizeof(double));
            condj = 0;
            for (size_t indj = 0; indj < boxCoord.size() - 2; indj++) {
                sys->conductorIn[si].x[condj] = boxCoord[indj] + xo;
                sys->conductorIn[si].y[condj] = boxCoord[indj + 1] + yo;
                if (sys->conductorIn[si].x[condj] > sys->conductorIn[si].xmax){
                    sys->conductorIn[si].xmax = sys->conductorIn[si].x[condj];
                }
                if (sys->conductorIn[si].x[condj] < sys->conductorIn[si].xmin){
                    sys->conductorIn[si].xmin = sys->conductorIn[si].x[condj];
                }
                if (sys->conductorIn[si].y[condj] > sys->conductorIn[si].ymax){
                    sys->conductorIn[si].ymax = sys->conductorIn[si].y[condj];
                }
                if (sys->conductorIn[si].y[condj] < sys->conductorIn[si].ymin){
                    sys->conductorIn[si].ymin = sys->conductorIn[si].y[condj];
                }

                indj++;
                condj++;
            }
        }
        cout << "  List of " << numText << " text boxes:" << endl;
        cout << "  List of " << numSRef << " structure references:" << endl;
        for (size_t indi = 0; indi < numSRef; indi++) // Handle each structure reference
        {
            printall((cell.sreferences)[indi].getSRefName(), (((cell.sreferences)[indi]).getSRefs())[0] + xo, (((cell.sreferences)[indi]).getSRefs())[1] + yo, sys);
        }

        // Close output file
        /*outfile.close();*/
    }

    // Print all the conductor information (no fdtdMesh)
    void printall(std::string name, double xo, double yo)
    {
        // Open output file
        /*std::ofstream outfile;*/
        size_t indExtension = this->getFileName().find(".", 1);
        /*std::string polyFileName = this->getFileName().substr(0, indExtension) + "_polygon.txt";
        outfile.open(polyFileName, std::ofstream::out | std::ofstream::app);
        if (!outfile.is_open()) // Failed to open the polygon file to write
        {
        return;
        }*/

        // Get information about this cell in ASCII database
        const GeoCell cell = this->cells[this->locateCell(name)];
        int numBound = cell.getNumBound();
        int numPath = cell.getNumPath();
        int numNode = cell.getNumNode();
        int numBox = cell.getNumBox();
        int numText = cell.getNumText();
        int numSRef = cell.getNumSRef();

        // Print cell information
        cout << "  List of " << numBound << " boundaries:" << endl;
        for (size_t indi = 0; indi < numBound; indi++) // Handle each boundary
        { 
            char point[128];
            vector<double> boundCoord = ((cell.boundaries)[indi]).getBounds();
            this->strPoints.append("    " + to_string(boundCoord.size() / 2 - 1) + " " + to_string(((cell.boundaries)[indi]).getLayer()) + " "); // Number of nodes, then layer number
            int lay = ((cell.boundaries)[indi]).getLayer();
            for (size_t indj = 0; indj < boundCoord.size() - 2; indj++) // C string for each ordered pair, -2 because the last point is the starting point
            {
                sprintf(point, "%1.4g, %1.4g ", boundCoord[indj++] + xo, boundCoord[indj + 1] + yo);
                this->strPoints.append(point);
            }
            this->strPoints.append("\n");
            /*outfile << strPoints;*/
            (this->numCdtIn)++;
        }
        cout << "  List of " << numPath << " paths:" << endl;
        for (size_t indi = 0; indi < numPath; indi++) // Handle each path
        {
            char point[128];
            this->strPoints.append("    4" + to_string(((cell.paths)[indi]).getLayer()) + " ");
            vector<double> pathCoord = ((cell.paths)[indi]).getPaths();
            double width = ((cell.paths)[indi]).getWidth();
            for (size_t indj = 0; indj < pathCoord.size(); indj++) // C string for each ordered pair
            {
                if (indj < pathCoord.size() - 2)
                {
                    if (pathCoord[indj] == pathCoord[indj + 2]) // along y axis
                    {
                        if (pathCoord[indj + 1] > pathCoord[indj + 3]) // first point is on top of the second point
                        {
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            this->strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            this->strPoints.append(point);
                            indj++;
                            indj++;
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            this->strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            this->strPoints.append(point);
                            this->strPoints.append("\n");
                            this->numCdtIn++;
                            indj--;
                        }
                        else // second point is on the top of the first point
                        {
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            this->strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            this->strPoints.append(point);
                            indj++;
                            indj++;
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            this->strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            this->strPoints.append(point);
                            this->strPoints.append("\n");
                            this->numCdtIn++;
                            indj--;
                        }
                    }
                    else // along x axis
                    {
                        if (pathCoord[indj] > pathCoord[indj + 2]) // first point is on the right of the second point
                        {
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            this->strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            this->strPoints.append(point);
                            indj++;
                            indj++;
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            this->strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            this->strPoints.append(point);
                            this->strPoints.append("\n");
                            this->numCdtIn++;
                            indj--;
                        }
                        else // second point is on the right of the first point
                        {
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] - width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            strPoints.append(point);
                            indj++;
                            indj++;
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] - width / 2 + yo);
                            strPoints.append(point);
                            sprintf(point, "%1.4g %1.4g ", pathCoord[indj] + width / 2 + xo, pathCoord[indj + 1] + width / 2 + yo);
                            strPoints.append(point);
                            strPoints.append("\n");
                            this->numCdtIn++;
                            indj--;
                        }
                    }
                }
            }
            /*outfile << strPoints;*/
        }
        cout << "  List of " << numNode << " nodes:" << endl;
        cout << "  List of " << numBox << " box outlines:" << endl;
        for (size_t indi = 0; indi < numBox; indi++) // Handle each box outline
        {
            /*outfile << "    " << ((cell.boxes)[indi]).getNBoxPt() - 1 << " " << ((cell.boxes)[indi]).getLayer() << " ";*/
            char point[128];
            this->strPoints.append("    " + to_string(((cell.boxes)[indi]).getNBoxPt() - 1) + " " + to_string(((cell.boxes)[indi]).getLayer()) + " ");
            vector<double> boxCoord = ((cell.boxes)[indi]).getBoxes();
            for (size_t indj = 0; indj < boxCoord.size() - 2; indj++) {
                /*outfile << boxCoord[indj++] << " " << boxCoord[indj] << " ";*/
                //this->strPoints.append(to_string(boxCoord[indj++]) + " " + to_string(boxCoord[indj + 1]) + " ");
                sprintf(point, "%1.4g %1.4g ", boxCoord[indj++] + xo, boxCoord[indj + 1] + yo);
                this->strPoints.append(point);
            }
            this->strPoints.append("\n");
            /*outfile << endl;*/
            this->numCdtIn++;
        }
        cout << "  List of " << numText << " text boxes:" << endl;
        cout << "  List of " << numSRef << " structure references:" << endl;
        for (size_t indi = 0; indi < numSRef; indi++) // Handle each structure reference
        {
            printall((cell.sreferences)[indi].getSRefName(), (((cell.sreferences)[indi]).getSRefs())[0] + xo, (((cell.sreferences)[indi]).getSRefs())[1] + yo);
        }

        // Close output file
        /*outfile.close();*/
    }

    // Print the ASCII database with the design geometry (fdtdMesh style)
    void print(vector<size_t> indCellPrint, fdtdMesh *sys)
    {
        // Delete existing file
        size_t indExtension = this->getFileName().find(".", 1);
        /*std::string polyFileName = this->getFileName().substr(0, indExtension) + "_polygon.txt";
        remove(polyFileName.c_str());*/

        // Analyze design and print to terminal
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
            std::string cellName = ((this->cells)[indCellPrint[indi]]).getCellName();
            cout << cellName << endl;
            this->printall(cellName, 0., 0., sys);   // the origin is the (0,0) point
            //(this->cells)[indCellPrint[indi]].printAlt();
        }
        cout << "------" << endl;
    }

    // Print the ASCII database with the design geometry (no fdtdMesh)
    void print(vector<size_t> indCellPrint)
    {
        // Delete existing file
        size_t indExtension = this->getFileName().find(".", 1);
        /*std::string polyFileName = this->getFileName().substr(0, indExtension) + "_polygon.txt";
        remove(polyFileName.c_str());*/

        // Analyze design and print to terminal
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
            std::string cellName = ((this->cells)[indCellPrint[indi]]).getCellName();
            cout << cellName << endl;
            this->printall(cellName, 0., 0.);   // the origin is the (0,0) point
            //(this->cells)[indCellPrint[indi]].printAlt();
        }
        cout << this->strPoints << endl;
        cout << "------" << endl;
    }

    // Dump ASCII database back into a GDSII file
    bool dump()
    {
        // Attempt to create GdsWriter variable
        GdsParser::GdsWriter gw((this->fileName).c_str());

        // Create library
        gw.create_lib(this->getLibName().c_str(), this->getdbUserUnits(), this->getdbUnits());

        // Create each cell
        size_t numCell = this->getNumCell();
        for (size_t indCell = 0; indCell < numCell; indCell++)
        {
            // Cell header
            gw.gds_write_bgnstr();
            gw.gds_write_strname(((this->cells)[indCell]).getCellName().c_str());

            // Create each boundary within the cell
            size_t numBound = ((this->cells)[indCell]).getNumBound();
            for (size_t indBound = 0; indBound < numBound; indBound++)
            {
                size_t numPt = (((this->cells)[indCell]).boundaries)[indBound].getNBoundPt();
                vector<double> bounds = (((this->cells)[indCell]).boundaries)[indBound].getBounds();
                bool has_last = true; // Boundaries always stored with first point repeated
                int xcoord[numPt];
                int ycoord[numPt];

                // Save boundary points to integer vectors (might not work)
                for (size_t indPt = 0; indPt < numPt; indPt++)
                {
                    xcoord[indPt] = (int)(bounds[2 * indPt + 0] / this->getdbUnits());
                    ycoord[indPt] = (int)(bounds[2 * indPt + 1] / this->getdbUnits());
                }

                // Write the boundary to file
                gw.gds_write_boundary();
                gw.gds_write_layer((((this->cells)[indCell]).boundaries)[indBound].getLayer());
                gw.gds_write_datatype(0);
                gw.gds_write_xy(xcoord, ycoord, numPt, has_last);
                gw.gds_write_endel();
            }

            // Create each path within the cell
            size_t numPath = ((this->cells)[indCell]).getNumPath();
            for (size_t indPath = 0; indPath < numPath; indPath++)
            {
                size_t numPt = (((this->cells)[indCell]).paths)[indPath].getNPathPt();
                vector<double> paths = (((this->cells)[indCell]).paths)[indPath].getPaths();
                int width = (int)((((this->cells)[indCell]).paths)[indPath].getWidth() / this->getdbUnits());
                int xcoord[numPt];
                int ycoord[numPt];

                // Save path points to integer vectors (might not work)
                for (size_t indPt = 0; indPt < numPt; indPt++)
                {
                    xcoord[indPt] = (int)(paths[2 * indPt + 0] / this->getdbUnits());
                    ycoord[indPt] = (int)(paths[2 * indPt + 1] / this->getdbUnits());
                }

                // Write the path to file
                gw.gds_write_path();
                gw.gds_write_layer((((this->cells)[indCell]).paths)[indPath].getLayer());
                gw.gds_write_datatype(0);
                gw.gds_write_pathtype((((this->cells)[indCell]).paths)[indPath].getType());
                gw.gds_write_width(width);
                gw.gds_write_xy(xcoord, ycoord, numPt);
                gw.gds_write_endel();
            }

            // Create each node within the cell (unfortunately unimplemented in GdsWriter.h)

            // Create each box outline within the cell
            size_t numBox = ((this->cells)[indCell]).getNumBox();
            for (size_t indBox = 0; indBox < numBox; indBox++)
            {
                size_t numPt = (((this->cells)[indCell]).boxes)[indBox].getNBoxPt();
                vector<double> boxes = (((this->cells)[indCell]).boxes)[indBox].getBoxes();
                bool has_last = true; // Box outlines always stored with first point repeated
                int xcoord[numPt];
                int ycoord[numPt];

                // Save box outline points to integer vectors (might not work)
                for (size_t indPt = 0; indPt < numPt; indPt++)
                {
                    xcoord[indPt] = (int)(boxes[2 * indPt + 0] / this->getdbUnits());
                    ycoord[indPt] = (int)(boxes[2 * indPt + 1] / this->getdbUnits());
                }

                // Write the box outline to file
                gw.gds_write_box();
                gw.gds_write_layer((((this->cells)[indCell]).boxes)[indBox].getLayer());
                gw.gds_write_boxtype((((this->cells)[indCell]).boxes)[indBox].getType());
                gw.gds_write_xy(xcoord, ycoord, numPt, has_last);
                gw.gds_write_endel();
            }

            // Create each textbox within the cell
            size_t numText = ((this->cells)[indCell]).getNumText();
            for (size_t indText = 0; indText < numText; indText++)
            {
                int xcoord[1] = { (int)(((((this->cells)[indCell]).textboxes)[indText].getTexts()[0]) / this->getdbUnits()) };
                int ycoord[1] = { (int)(((((this->cells)[indCell]).textboxes)[indText].getTexts()[1]) / this->getdbUnits()) };
                int width = (int)((((this->cells)[indCell]).textboxes)[indText].getWidth() / this->getdbUnits());

                // Write the textbox to file
                gw.gds_write_text();
                gw.gds_write_layer((((this->cells)[indCell]).textboxes)[indText].getLayer());
                gw.gds_write_texttype((((this->cells)[indCell]).textboxes)[indText].getType());
                gw.gds_write_presentation((((this->cells)[indCell]).textboxes)[indText].getFontID(), (((this->cells)[indCell]).textboxes)[indText].getJusts()[0], (((this->cells)[indCell]).textboxes)[indText].getJusts()[1]);
                gw.gds_write_width(width);
                gw.gds_write_strans(0, 0, 0); // Strans unimplemented until needed (reflect, mag, angle)
                gw.gds_write_xy(xcoord, ycoord, 1);
                gw.gds_write_string((((this->cells)[indCell]).textboxes)[indText].getTextStr().c_str());
                gw.gds_write_endel();
            }

            // Create each structure reference within the cell
            size_t numSRef = ((this->cells)[indCell]).getNumSRef();
            for (size_t indSRef = 0; indSRef < numSRef; indSRef++)
            {
                int xcoord[1] = { (int)(((((this->cells)[indCell]).sreferences)[indSRef].getSRefs()[0]) / this->getdbUnits()) };
                int ycoord[1] = { (int)(((((this->cells)[indCell]).sreferences)[indSRef].getSRefs()[1]) / this->getdbUnits()) };

                // Write the structure reference to file
                gw.gds_write_sref();
                gw.gds_write_sname((((this->cells)[indCell]).sreferences)[indSRef].getSRefName().c_str());
                //gw.gds_write_mag(1.0); // MAG unimplemented until needed
                //gw.gds_write_angle(0.0); // ANGLE unimplemented until needed
                gw.gds_write_xy(xcoord, ycoord, 1);
                gw.gds_write_endel();
            }

            // Cell closer
            gw.gds_write_endstr();
        }

        // Close library
        gw.gds_write_endlib();
        return true;
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
            for (size_t indj = 0; indj < data.size() + 1; indj++) // Only store printable characters
            {
                if (((int)data[indj] < 32) || ((int)data[indj] > 128))
                {
                    label[indj] = '\0';
                    break; // Break if unprintable character
                }
                label[indj] = (char)data[indj];
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
            for (size_t indj = 0; indj < data.size() + 1; indj++) // Only store printable characters
            {
                if (((int)data[indj] < 32) || ((int)data[indj] > 128))
                {
                    label[indj] = '\0';
                    break; // Break if unprintable character
                }
                label[indj] = (char)data[indj];
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
                boundary modBound = ((getCell(this->numCell)).boundaries).back(); // Get copy of boundary
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
            vector<double> coord;
            for (size_t indi = 0; indi < data.size(); indi++)
            {
                coord.push_back(unitFactor * data[indi]);
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
                ((this->cells)[this->numCell]).textboxes.emplace_back(textbox(coord, modText.getLayer(), modText.getProps(), modText.getType(), modText.getFontID(), modText.getJusts(), modText.getWidth(), modText.getTextStr()));
            }
            else if (this->getElement() == 's')
            {
                sref modSRef = getCell(this->numCell).sreferences.back();
                ((this->cells)[this->numCell]).sreferences.pop_back();
                ((this->cells)[this->numCell]).sreferences.emplace_back(sref(coord, modSRef.getProps(), modSRef.getSRefName()));
            }
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
            else if (this->getElement() == 't')
            {
                textbox modText = getCell(this->numCell).textboxes.back();
                ((this->cells)[this->numCell]).textboxes.pop_back();
                ((this->cells)[this->numCell]).textboxes.emplace_back(textbox(modText.getTexts(), modText.getLayer(), modText.getProps(), modText.getType(), modText.getFontID(), modText.getJusts(), unitFactor * data[0], modText.getTextStr()));
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
            for (size_t indj = 0; indj < data.size() + 1; indj++) // Only store printable characters
            {
                if (((int)data[indj] < 32) || ((int)data[indj] > 128))
                {
                    label[indj] = '\0';
                    break; // Break if unprintable character
                }
                label[indj] = (char)data[indj];
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
            for (size_t indj = 0; indj < data.size() + 1; indj++) // Only store printable characters
            {
                if (((int)data[indj] < 32) || ((int)data[indj] > 128))
                {
                    label[indj] = '\0';
                    break; // Break if unprintable character
                }
                label[indj] = (char)data[indj];
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
            for (size_t indj = 0; indj < data.size() + 1; indj++) // Only store printable characters
            {
                if (((int)data[indj] < 32) || ((int)data[indj] > 128))
                {
                    label[indj] = '\0';
                    break; // Break if unprintable character
                }
                label[indj] = (char)data[indj];
            }
            //cout << "Element property value: " << label << endl;

            // Store  property (might only handle a single one)
            if (this->getElement() == 'b')
            {
                boundary modBound = getCell(this->numCell).boundaries.back(); // Get copy of boundary
                vector<std::string> modProps = {}; // Use empty vector for new properties, except ...
                if (modBound.getProps().size() != 0)
                {
                    modProps = modBound.getProps(); // ... Replace with any existing properties of this boundary (to be used)
                }
                modProps.assign(this->numProp, ""); // Ensure enough elements available for property to go in correct location
                modProps[this->numProp - 1] = label; // Put new property value in place
                ((this->cells)[this->numCell]).boundaries.pop_back(); // Remove last boundary vector entry
                ((this->cells)[this->numCell]).boundaries.emplace_back(boundary(modBound.getBounds(), modBound.getLayer(), modProps)); // Put boundary back with properties update
            }
            else if (this->getElement() == 'p')
            {
                path modPath = getCell(this->numCell).paths.back();
                vector<std::string> modProps = {};
                if (modPath.getProps().size() != 0)
                {
                    modProps = modPath.getProps();
                }
                modProps.assign(this->numProp, "");
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).paths.pop_back();
                ((this->cells)[this->numCell]).paths.emplace_back(path(modPath.getPaths(), modPath.getLayer(), modProps, modPath.getType(), modPath.getWidth()));
            }
            else if (this->getElement() == 'n')
            {
                node modNode = getCell(this->numCell).nodes.back();
                vector<std::string> modProps = {};
                if (modNode.getProps().size() != 0)
                {
                    modProps = modNode.getProps();
                }
                modProps.assign(this->numProp, "");
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).nodes.pop_back();
                ((this->cells)[this->numCell]).nodes.emplace_back(node(modNode.getNodes(), modNode.getLayer(), modProps, modNode.getType()));
            }
            else if (this->getElement() == 'x')
            {
                box modBox = getCell(this->numCell).boxes.back();
                vector<std::string> modProps = {};
                if (modBox.getProps().size() != 0)
                {
                    modProps = modBox.getProps();
                }
                modProps.assign(this->numProp, "");
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).boxes.pop_back();
                ((this->cells)[this->numCell]).boxes.emplace_back(box(modBox.getBoxes(), modBox.getLayer(), modProps, modBox.getType()));
            }
            else if (this->getElement() == 't')
            {
                textbox modText = getCell(this->numCell).textboxes.back();
                vector<std::string> modProps = {};
                if (modText.getProps().size() != 0)
                {
                    modProps = modText.getProps();
                }
                modProps.assign(this->numProp, "");
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).textboxes.pop_back();
                ((this->cells)[this->numCell]).textboxes.emplace_back(textbox(modText.getTexts(), modText.getLayer(), modProps, modText.getType(), modText.getFontID(), modText.getJusts(), modText.getWidth(), modText.getTextStr()));
            }
            else if (this->getElement() == 's')
            {
                sref modSRef = getCell(this->numCell).sreferences.back();
                vector<std::string> modProps = {};
                if (modSRef.getProps().size() != 0)
                {
                    modProps = modSRef.getProps();
                }
                modProps.assign(this->numProp, "");
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).sreferences.pop_back();
                ((this->cells)[this->numCell]).sreferences.emplace_back(sref(modSRef.getSRefs(), modProps, modSRef.getSRefName()));
            }
        }
        else if (ascii_record_type == "ENDEL")
        {
            if (this->getElement() == 'b')
            {
                ((this->cells)[this->numCell]).boundaries.back().reorder(); // Put boundary points in preferred order
            }
            else if (this->getElement() == 'n')
            {
                ((this->cells)[this->numCell]).nodes.back().reorder(); // Put node points in preferred order
            }
            else if (this->getElement() == 'x')
            {
                ((this->cells)[this->numCell]).boxes.back().reorder(); // Put box points in preferred order
            }
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
