/**
* @file   limboint.h
* @author Michael R. Hayashi
* @date   18 October 2018
* @brief  Header for Limbo/Geometry Handling Interface for GDSII Files
*/

// Make sure this file is included only once
#ifndef limboint_h
#define limboint_h

#define _USE_MATH_DEFINES // Place before including <cmath> for e, log2(e), log10(e), ln(2), ln(10), pi, pi/2, pi/4, 1/pi, 2/pi, 2/sqrt(pi), sqrt(2), and 1/sqrt(2)
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
#include <stack>
#include <map>
#include <tr1/unordered_map>
#include <cmath>
#include <complex>
#include <limbo/parsers/gdsii/stream/GdsReader.h>
#include <limbo/parsers/gdsii/stream/GdsWriter.h>
#include "fdtd.h"

// PSLG macro
#define REL_INTERIOR_PRECISION (0.001) // Relative preceision to set interior points of a PSLG region

// Structure for convex hull comparison
struct compareObj {
    complex<double> p0; // Reference point

    // Binary compare function for convex hull (complex numbers by orientation then by distance)
    bool operator()(complex<double> p1, complex<double> p2)
    {
        // Find orientation angle with cross product of displacements
        double crossProd = (p1.imag() - this->p0.imag()) * (p2.real() - p1.real()) - (p1.real() - this->p0.real()) * (p2.imag() - p1.imag());
        int orientation = 0; // Assume collinear points
        if (crossProd > 0.)
        {
            orientation = 1; // Clockwise turn, so first point should NOT be ahead of second point
            return false;
        }
        else if (crossProd < 0.)
        {
            orientation = 2; // Counterclockwise turn, so first point should be ahead of second point
            return true;
        }
        else
        {
            // Must use square distance of displacements as tie-breaker for collinear points
            complex<double> disp1 = p1 - this->p0; // Displacement to first point
            complex<double> disp2 = p2 - this->p0; // Displacement to second point
            if (norm(disp2) >= norm(disp1))
            {
                return true; // Second point is further, so order is good
            }
            else
            {
                return false; // First point is further, so order needs swapping
            }
        }
    }
} convexCompare;

// Structure for planar straight-line graphs
struct pslg
{
public:
    vector<complex<double>> vertices;      // Vector of all (unique) vertex coordinates (m, x + iy format)
    vector<pair<size_t, size_t>> segments; // Vector of all segments as pairs connecting vertex indices
    vector<complex<double>> regions;       // Vector of all conductor region coordinates (m, x + iy format)
public:
    // Default constructor
    pslg()
    {
        this->vertices = {};
        this->segments = {};
        this->regions = {};
    }

    // Parametrized constructor
    pslg(vector<complex<double>> vertices, vector<pair<size_t, size_t>> segments, vector<complex<double>> regions)
    {
        this->vertices = vertices;
        this->segments = segments;
        this->regions = regions;
    }

    // Destructor
    ~pslg()
    {
        // Nothing
    }
};

// Helper class for linear transformations (reflection, rotation, and scaling) of GDSII file elements
class strans
{
private:
    bool mirrored;                 // Is the structure mirrored in the x-direction (before rotation)?
    bool absMagnify;               // Should the magnification be treated as absolute rather than as nested effects?
    bool absRotate;                // Should the rotation be treated as absolute rather than as nested effects?
    double magnify;                // Magnification of structure
    double rotation;               // Rotation angle of structure (rad)
public:
    // Default constructor
    strans()
    {
        // Default values for reflection, rotation, and scaling (note that squeezing, shearing, and projection are not included)
        this->mirrored = false;
        this->absMagnify = false;
        this->absRotate = false;
        this->magnify = 1.0;
        this->rotation = 0.0;
    }

    // Parametrized constructor
    strans(bool mirrored, bool absMagnify, bool absRotate, double magnify, double rotation)
    {
        this->mirrored = mirrored;
        this->absMagnify = absMagnify;
        this->absRotate = absRotate;
        this->magnify = magnify;
        this->rotation = rotation;
    }

    // Get if structure is mirrored in the x-direction before rotation
    bool getMirrored() const
    {
        return this->mirrored;
    }

    // Get if structure magnification should be treated as absolute rather than compounding
    bool getAbsMagnify() const
    {
        return this->absMagnify;
    }

    // Get if structure rotation should be treated as absolute rather than compounding
    bool getAbsRotate() const
    {
        return this->absRotate;
    }

    // Get magnification of structure
    double getMagnify() const
    {
        return this->magnify;
    }

    // Get rotation angle of structure after any mirroring (rad)
    double getRotation() const
    {
        return this->rotation;
    }

    // Set x-direction mirroring of structure before rotation
    void setMirrored(bool mirrored)
    {
        this->mirrored = mirrored;
    }

    // Set absolute interpretation of structure magnification
    void setAbsMagnify(bool absMagnify)
    {
        this->absMagnify = absMagnify;
    }

    // Set absolute interpretation of structure rotation
    void setAbsRotate(bool absRotate)
    {
        this->absRotate = absRotate;
    }

    // Set magnification of structure
    void setMagnify(double magnify)
    {
        this->magnify = magnify;
    }

    // Set rotation angle of structure to occur after any mirroring (rad)
    void setRotation(double rotation)
    {
        this->rotation = rotation;
    }

    // Compose this transformation with previous one to get net transformation
    strans composeTransform(strans oldTransform) const
    {
        /* Linear algebra combinatorial study:
        The rotation matrix is orthogonal, and the reflection matrix is orthogonal and involutory.
        Scaling has the effect of scalar multiplication and is simple to compose independent of other operations.
        The reflection matrix squared is the identity matrix because it is involutory. [2 factors: 1 scenario {1010}]
        The product of two different rotation matrices is a rotation matrix with the _sum of angles_ by trigonometric identities. [2 factors: 1 scenario {0101}]
        Reflection and then rotation is the rotation matrix with the first column negated (natural part of class). [2 factors: 3 scenarios {1100, 1001, 0011}]
        Rotation and then reflection is the rotation matrix with the first row negated (interpret as reflection then rotation by _negative_ angle). [2 factors: 1 scenario {0110}]
        Reflection, rotation, and then reflection again causes only rotation by the _negative_ angle. [3 factors: 1 scenario {1110}]
        Reflection, rotation, and then a different rotation causes reflection and then rotation by _sum of angles_. [3 factors: 1 scenario {1101}]
        Reflection, reflection again, and then rotation causes only rotation by the _positive_ angle because of involutory reflections. [3 factors: 1 scenario {1011}]
        Rotation, reflection, and then a different rotation causes reflection and then rotation by _second angle minus first_. [3 factors: 1 scenario {0111}]
        Reflection, rotation, reflection again, and then a different rotation yields only rotation by _second angle minus first_. [4 factors: 1 scenario {1111}]*/

        // Find net scaling/magnification
        double newMag = this->magnify;
        if (!this->absMagnify)
        {
            newMag *= oldTransform.getMagnify(); // Compounding only if the _present_ transformation allows it
        }

        // Find if net transformation requires x-direction mirroring
        bool newMirror = false;
        if ((!this->mirrored && oldTransform.getMirrored()) || (this->mirrored && !oldTransform.getMirrored()))
        {
            // The net transformation only requires mirroring if there are an _odd_ number of mirroring factors
            newMirror = true;
        }

        // Find the net rotation angle (rad)
        double newRot = this->rotation;
        if (!this->absRotate)
        {
            // Modification of the rotation angle if a rotation has occurred before (avoiding minimal logic for clarity)
            if (!oldTransform.getMirrored() && (oldTransform.getRotation() != 0.) && !this->mirrored && (this->rotation == 0.))
            {
                // Reflection and then rotation {1100}
                newRot = oldTransform.getRotation(); // Positive old angle (natural part of class, but need only this angle)
            }
            else if (!oldTransform.getMirrored() && (oldTransform.getRotation() != 0.) && !this->mirrored && (this->rotation != 0.))
            {
                // Two different rotations {0101}
                newRot += oldTransform.getRotation(); // Sum of angles
            }
            else if (!oldTransform.getMirrored() && (oldTransform.getRotation() != 0.) && this->mirrored && (this->rotation == 0.))
            {
                // Rotation and then reflection {0110}
                newRot = -1. * oldTransform.getRotation(); // Negative old angle
            }
            else if (oldTransform.getMirrored() && (oldTransform.getRotation() != 0.) && this->mirrored && (this->rotation == 0.))
            {
                // Reflection, rotation, and then reflection again {1110}
                newRot = -1. * oldTransform.getRotation(); // Negative old angle
            }
            else if (oldTransform.getMirrored() && (oldTransform.getRotation() != 0.) && !this->mirrored && (this->rotation != 0.))
            {
                //  Reflection, rotation, and then a different rotation {1101}
                newRot += oldTransform.getRotation(); // Sum of angles
            }
            else if (!oldTransform.getMirrored() && (oldTransform.getRotation() != 0.) && this->mirrored && (this->rotation != 0.))
            {
                // Rotation, reflection, and then a different rotation {0111}
                newRot -= oldTransform.getRotation(); // New angle minus old angle
            }
            else if (oldTransform.getMirrored() && (oldTransform.getRotation() != 0.) && this->mirrored && (this->rotation != 0.))
            {
                // Reflection, rotation, reflection again, and then a different rotation {1111}
                newRot -= oldTransform.getRotation(); // New angle minus old angle
            }
        } // Modification of the rotation angle is not possible without prior rotation
        newRot = fmod(newRot + M_PI, 2. * M_PI) - M_PI; // Wrap angle to interval [-pi, pi)

        // Return the composed linear transformation
        return strans(newMirror, this->absMagnify, this->absRotate, newMag, newRot);
    }

    // Apply this transformation to a point
    vector<double> applyTranform(vector<double> oldPt)
    {
        // Error checking on input point
        if (oldPt.size() != 2)
        {
            cerr << "Coordinates of point to transform must be a length-2 vector. Returning same point unaltered." << endl;
            return oldPt;
        }
        vector<double> newPt = oldPt;

        // Apply any x-axis mirroring
        if (this->mirrored)
        {
            newPt[1] = -newPt[1];
        }

        // Apply any rotation
        if (this->rotation != 0.)
        {
            vector<double> workPt = newPt; // Working point as temporary variable
            newPt[0] = +cos(this->rotation) * workPt[0] + -sin(this->rotation) * workPt[1]; // First row of rotation matrix-vector multiplication
            newPt[1] = +sin(this->rotation) * workPt[0] + +cos(this->rotation) * workPt[1]; // Second row of rotation matrix-vector multiplication
        }

        // Apply any magnification
        if (this->magnify != 1.)
        {
            newPt[0] *= this->magnify;
            newPt[1] *= this->magnify;
        }

        // Return the new point
        return newPt;
    }

    // Print function
    void print()
    {
        cout << "   Linear transformation: mirrored (" << (this->mirrored ? "Y" : "N") << "), magnification of " << this->magnify << "x (" << (this->absMagnify ? "absolute" : "relative") << "), rotation of " << this->rotation << " rad (" << (this->absRotate ? "absolute" : "relative") << ")" << endl;
    }

    // Destructor
    ~strans()
    {
        // Nothing
    }
};

// Custom classes for containing the 7 GDSII file elements
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

    // Set boundary points in cell coordinates
    void setBounds(vector<double> bounds)
    {
        this->bounds = bounds;
    }

    // Set boundary layer number
    void setLayer(int layer)
    {
        this->layer = layer;
    }

    // Set boundary properties
    void setProps(vector<std::string> props)
    {
        this->props = props;
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

    // Return perimeter of boundary
    double findPerimeter() const
    {
        double polyPeri = 0.;
        for (size_t indj = 0; indj < this->getNBoundPt() - 1; indj++)
        {
            polyPeri += hypot(this->bounds[2 * indj + 2] - this->bounds[2 * indj], this->bounds[2 * indj + 3] - this->bounds[2 * indj + 1]); // Running calculation of perimeter
        }
        return polyPeri; // Polygon perimeter (m)
    }

    // Return area of boundary
    double findArea() const
    {
        double polyArea = 0.;
        for (size_t indj = 0; indj < this->getNBoundPt() - 1; indj++)
        {
            polyArea += (this->bounds[2 * indj] * this->bounds[2 * indj + 3]) - (this->bounds[2 * indj + 2] * this->bounds[2 * indj + 1]); // Running calculation of closed polygon area
        }
        polyArea /= 2.0; // Closed polygon area (m^2)
        return fabs(polyArea);
    }

    // Return centroid of boundary
    vector<double> findCentroid(double xo, double yo) const
    {
        double centroidX = 0.;
        double centroidY = 0.;
        double polyArea = 0.;
        for (size_t indj = 0; indj < this->getNBoundPt() - 1; indj++) // Handle each coordinate comprising boundary except repeated point
        {
            centroidX += (this->bounds[2 * indj] + this->bounds[2 * indj + 2] + 2 * xo) * ((this->bounds[2 * indj] + xo) * (this->bounds[2 * indj + 3] + yo) - (this->bounds[2 * indj + 2] + xo) * (this->bounds[2 * indj + 1] + yo)); // Running calculation of polygon centroid x-coordinate
            centroidY += (this->bounds[2 * indj + 1] + this->bounds[2 * indj + 3] + 2 * yo) * ((this->bounds[2 * indj] + xo) * (this->bounds[2 * indj + 3] + yo) - (this->bounds[2 * indj + 2] + xo) * (this->bounds[2 * indj + 1] + yo)); // Running calculation of polygon centroid y-coordinate
            polyArea += (this->bounds[2 * indj] + xo) * (this->bounds[2 * indj + 3] + yo) - (this->bounds[2 * indj + 2] + xo) * (this->bounds[2 * indj + 1] + yo); // Running calculation of closed polygon area
        }
        polyArea /= 2.0; // Closed polygon area (m^2)
        centroidX /= (6.0 * polyArea); // Polygon centroid x-coordinate (m)
        centroidY /= (6.0 * polyArea); // Polygon centroid y-coordinate (m)
        return{ centroidX, centroidY };
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
            nBoundPt--; // Correct number of boundary points for this method
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
        if (((this->bounds)[2 * 0] > (this->bounds)[2 * 1]) && ((this->bounds)[2 * 0 + 1] >= (this->bounds)[2 * 1 + 1]))
        {
            vector<double> revBound;
            revBound.push_back((this->bounds)[0]); // Same first x-coordinate
            revBound.push_back((this->bounds)[1]); // Same first y-coordinate
            for (size_t indi = nBoundPt - 1; indi > 0; indi--)
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
        //this->layer = 0;
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
            cerr << "Path type must be 0, 1, or 2. Defaulting to type 2." << endl;
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

    // Set path vertices in cell coordinates
    void setPaths(vector<double> paths)
    {
        this->paths = paths;
    }

    // Set path layer number
    void setLayer(int layer)
    {
        this->layer = layer;
    }

    // Set path properties
    void setProps(vector<std::string> props)
    {
        this->props = props;
    }

    // Set path type
    // 0 = square ends at vertices, 1 = round ends, 2 = square ends overshoot vertices by half width
    void setType(int type)
    {
        if ((type == 0) || (type == 1) || (type == 2))
        {
            this->type = type;
        }
        else
        {
            cerr << "Path type must be 0, 1, or 2. Defaulting to type 2." << endl;
            this->type = 2; // Default to type 2 if invalid
        }
    }

    // Set path width
    // Negative values mean independent of structure scaling
    void setWidth(double width)
    {
        this->width = width;
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

    // Return path length
    double findLength() const
    {
        double pathLength = 0.;
        for (size_t indj = 0; indj < this->getNPathPt() - 1; indj++)
        {
            pathLength += hypot(this->paths[2 * indj + 2] - this->paths[2 * indj], this->paths[2 * indj + 3] - this->paths[2 * indj + 1]); // Running calculation of path length
        }
        return pathLength; // Path length (m)
    }

    // Destructor
    ~path()
    {
        //this->paths = {};
        //this->layer = 0;
        //this->type = 2;
        //this->width = 0.;
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

    // Set node vertices in cell coordinates
    void setNodes(vector<double> nodes)
    {
        this->nodes = nodes;
    }

    // Set node layer number
    void setLayer(int layer)
    {
        this->layer = layer;
    }

    // Set node properties
    void setProps(vector<std::string> props)
    {
        this->props = props;
    }

    // Set node type
    void setType(int type)
    {
        this->type = type;
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

    // Return perimeter of node
    double findPerimeter() const
    {
        double nodePeri = 0.;
        for (size_t indj = 0; indj < this->getNNodePt() - 1; indj++)
        {
            nodePeri += hypot(this->nodes[2 * indj + 2] - this->nodes[2 * indj], this->nodes[2 * indj + 3] - this->nodes[2 * indj + 1]); // Running calculation of perimeter
        }
        return nodePeri; // Node perimeter (m)
    }

    // Return area of node
    double findArea() const
    {
        double nodeArea = 0.;
        for (size_t indj = 0; indj < this->getNNodePt() - 1; indj++)
        {
            nodeArea += (this->nodes[2 * indj] * this->nodes[2 * indj + 3]) - (this->nodes[2 * indj + 2] * this->nodes[2 * indj + 1]); // Running calculation of closed node area
        }
        nodeArea /= 2.0; // Closed node area (m^2)
        return fabs(nodeArea);
    }

    // Return centroid of node
    vector<double> findCentroid(double xo, double yo) const
    {
        double centroidX = 0.;
        double centroidY = 0.;
        double nodeArea = 0.;
        for (size_t indj = 0; indj < this->getNNodePt() - 1; indj++) // Handle each coordinate comprising node except repeated point
        {
            centroidX += (this->nodes[2 * indj] + this->nodes[2 * indj + 2] + 2 * xo) * ((this->nodes[2 * indj] + xo) * (this->nodes[2 * indj + 3] + yo) - (this->nodes[2 * indj + 2] + xo) * (this->nodes[2 * indj + 1] + yo)); // Running calculation of node centroid x-coordinate
            centroidY += (this->nodes[2 * indj + 1] + this->nodes[2 * indj + 3] + 2 * yo) * ((this->nodes[2 * indj] + xo) * (this->nodes[2 * indj + 3] + yo) - (this->nodes[2 * indj + 2] + xo) * (this->nodes[2 * indj + 1] + yo)); // Running calculation of node centroid y-coordinate
            nodeArea += (this->nodes[2 * indj] + xo) * (this->nodes[2 * indj + 3] + yo) - (this->nodes[2 * indj + 2] + xo) * (this->nodes[2 * indj + 1] + yo); // Running calculation of closed node area
        }
        nodeArea /= 2.0; // Closed node area (m^2)
        centroidX /= (6.0 * nodeArea); // Node centroid x-coordinate (m)
        centroidY /= (6.0 * nodeArea); // Node centroid y-coordinate (m)
        return{ centroidX, centroidY };
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
            nNodePt--; // Correct number of node points for this method
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
        if (((this->nodes)[2 * 0] > (this->nodes)[2 * 1]) && ((this->nodes)[2 * 0 + 1] >= (this->nodes)[2 * 1 + 1]))
        {
            vector<double> revNode;
            revNode.push_back((this->nodes)[0]); // Same first x-coordinate
            revNode.push_back((this->nodes)[1]); // Same first y-coordinate
            for (size_t indi = nNodePt - 1; indi > 0; indi--)
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
        //this->layer = 0;
        //this->type = 0;
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

    // Set box outline vertices in cell coordinates
    void setBoxes(vector<double> boxes)
    {
        this->boxes = boxes;
    }

    // Set box outline layer number
    void setLayer(int layer)
    {
        this->layer = layer;
    }

    // Set box outline properties
    void setProps(vector<std::string> props)
    {
        this->props = props;
    }

    // Set box outline type
    void setType(int type)
    {
        this->type = type;
    }

    // Return number of box outline points (should be 5 for closed box with recommended repetition of first point or 4 without repetition)
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

    // Return perimeter of box outline
    double findPerimeter() const
    {
        double boxPeri = 0.;
        for (size_t indj = 0; indj < this->getNBoxPt() - 1; indj++)
        {
            boxPeri += hypot(this->boxes[2 * indj + 2] - this->boxes[2 * indj], this->boxes[2 * indj + 3] - this->boxes[2 * indj + 1]); // Running calculation of perimeter
        }
        return boxPeri; // Box outline perimeter (m)
    }

    // Return area of box
    double findArea() const
    {
        double boxArea = 0.;
        for (size_t indj = 0; indj < this->getNBoxPt() - 1; indj++)
        {
            boxArea += (this->boxes[2 * indj] * this->boxes[2 * indj + 3]) - (this->boxes[2 * indj + 2] * this->boxes[2 * indj + 1]); // Running calculation of closed box area
        }
        boxArea /= 2.0; // Closed box area (m^2)
        return fabs(boxArea);
    }

    // Return centroid of box outline
    vector<double> findCentroid(double xo, double yo) const
    {
        double centroidX = 0.;
        double centroidY = 0.;
        double boxArea = 0.;
        for (size_t indj = 0; indj < this->getNBoxPt() - 1; indj++) // Handle each coordinate comprising box outline except repeated point
        {
            centroidX += (this->boxes[2 * indj] + this->boxes[2 * indj + 2] + 2 * xo) * ((this->boxes[2 * indj] + xo) * (this->boxes[2 * indj + 3] + yo) - (this->boxes[2 * indj + 2] + xo) * (this->boxes[2 * indj + 1] + yo)); // Running calculation of box outline centroid x-coordinate
            centroidY += (this->boxes[2 * indj + 1] + this->boxes[2 * indj + 3] + 2 * yo) * ((this->boxes[2 * indj] + xo) * (this->boxes[2 * indj + 3] + yo) - (this->boxes[2 * indj + 2] + xo) * (this->boxes[2 * indj + 1] + yo)); // Running calculation of box outline centroid y-coordinate
            boxArea += (this->boxes[2 * indj] + xo) * (this->boxes[2 * indj + 3] + yo) - (this->boxes[2 * indj + 2] + xo) * (this->boxes[2 * indj + 1] + yo); // Running calculation of closed box area
        }
        boxArea /= 2.0; // Closed box area (m^2)
        centroidX /= (6.0 * boxArea); // Box outline centroid x-coordinate (m)
        centroidY /= (6.0 * boxArea); // Box outline centroid y-coordinate (m)
        return{ centroidX, centroidY };
    }

    // Reorder points with lower-right first, going counterclockwise
    void reorder()
    {
        // Remove last point if equal to first point
        size_t nBoxPt = this->getNBoxPt(); // Number of box points (coordinate pairs)
        if (((this->boxes)[2 * 0] == (this->boxes)[2 * (nBoxPt - 1)]) && ((this->boxes)[2 * 0 + 1] == (this->boxes)[2 * (nBoxPt - 1) + 1]))
        {
            (this->boxes).pop_back(); // Remove repeated y-coordinate
            (this->boxes).pop_back(); // Remove repeated x-coordinate
            nBoxPt--; // Correct number of box points for this method
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
        if (((this->boxes)[2 * 0] > (this->boxes)[2 * 1]) && ((this->boxes)[2 * 0 + 1] >= (this->boxes)[2 * 1 + 1]))
        {
            vector<double> revBox;
            revBox.push_back((this->boxes)[0]); // Same first x-coordinate
            revBox.push_back((this->boxes)[1]); // Same first y-coordinate
            for (size_t indi = nBoxPt - 1; indi > 0; indi--)
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
        //this->layer = 0;
        //this->type = 0;
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
    strans transform;              // Linear transformations (reflection, rotation, and scaling) of text
    std::string textStr;           // String in text box
public:
    // Default constructor
    textbox()
    {
        vector<double> texts = { 0., 0. };
        vector<int> justs = { 0, 0 };
        strans transform;
        this->texts = texts;
        this->layer = 0;
        this->props = {};
        this->type = 0;
        this->fontID = 0;
        this->justs = justs;
        this->width = 0.;
        this->transform = transform;
        this->textStr = "";
    }

    // Parametrized constructor
    textbox(vector<double> texts, int layer, vector<std::string> props, int type, int fontID, vector<int> justs, double width, strans transform, std::string textStr)
    {
        if (texts.size() != 2)
        {
            cerr << "Text box center coordinates must be stored in vector of length 2. Defaulting to origin (0, 0)." << endl;
            this->texts = { 0., 0. };
        }
        else
        {
            this->texts = texts;
        }
        this->layer = layer;
        this->props = props;
        this->type = type;
        this->fontID = fontID;
        if (justs.size() != 2)
        {
            cerr << "Justifications must be stored in vector of size 2. Defaulting to {1, 1} (middle, center)." << endl;
            this->justs = { 1, 1 };
        }
        else
        {
            if ((justs[0] < 0) || (justs[0] > 2) || (justs[1] < 0) || (justs[1] > 2))
            {
                cerr << "Justification values must each be 0, 1, or 2. Defaulting to {1, 1} (middle, center)." << endl;
                this->justs = { 1, 1 };
            }
            else
            {
                this->justs = justs;
            }
        }
        this->width = width;
        this->transform = transform;
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

    // Get the text linear transformations (reflection, rotation, and scaling)
    strans getTransform() const
    {
        return this->transform;
    }

    // Get string in text box
    std::string getTextStr() const
    {
        return this->textStr;
    }

    // Set text box center in cell coordinates
    void setTexts(vector<double> texts)
    {
        this->texts = texts;
    }

    // Set text box layer number
    void setLayer(int layer)
    {
        this->layer = layer;
    }

    // Set text box properties
    void setProps(vector<std::string> props)
    {
        this->props = props;
    }

    // Set text box type
    void setType(int type)
    {
        this->type = type;
    }

    // Set text box font ID (0-3)
    void setFontID(int fontID)
    {
        this->fontID = fontID;
    }

    // Set text box justifications
    // First (vertical): 0 = top, 1 = middle, 2 = bottom
    // Second (horizontal): 0 = left, 1 = center, 2 = right
    void setJusts(vector<int> justs)
    {
        if (justs.size() != 2)
        {
            cerr << "Justifications must be stored in vector of size 2. Defaulting to {1, 1} (middle, center)." << endl;
            this->justs = { 1, 1 };
        }
        else
        {
            if ((justs[0] < 0) || (justs[0] > 2) || (justs[1] < 0) || (justs[1] > 2))
            {
                cerr << "Justification values must each be 0, 1, or 2. Defaulting to {1, 1} (middle, center)." << endl;
                this->justs = { 1, 1 };
            }
            else
            {
                this->justs = justs;
            }
        }
    }

    // Set text width
    // Negative values mean independent of structure scaling
    void setWidth(double width)
    {
        this->width = width;
    }

    // Set text linear transformations (reflection, rotation, and scaling)
    void setTransform(strans transform)
    {
        this->transform = transform;
    }

    // Get string in text box
    void setTextStr(std::string textStr)
    {
        this->textStr = textStr;
    }

    // Destructor
    ~textbox()
    {
        /*vector<double> texts = { 0., 0. };
        vector<int> justs = { 0, 0 };
        strans transform;
        this->texts = texts;
        this->layer = 0;
        this->type = 2;
        this->fontID = 0;
        this->justs = justs;
        this->width = 0.;
        this->transform = transform;
        this->textStr = "";*/
    }
};

class sref
{
private:
    vector<double> srefs;          // Coordinate pair of placed structure reference
    vector<std::string> props;     // Properties of structure reference
    std::string srefName;          // Name of structure reference
    strans transform;              // Linear transformations (reflection, rotation, and scaling) of structure reference
public:
    // Default constructor
    sref()
    {
        vector<double> srefs = { 0., 0. };
        strans transform;
        this->srefs = srefs;
        this->props = {};
        this->srefName = "";
        this->transform = transform;
    }

    // Parametrized constructor
    sref(vector<double> srefs, vector<std::string> props, std::string srefName, strans transform)
    {
        if (srefs.size() != 2)
        {
            cerr << "Structure reference center coordinates must be stored in vector of length 2. Defaulting to origin (0, 0)." << endl;
            this->srefs = { 0., 0. };
        }
        else
        {
            this->srefs = srefs;
        }
        this->props = props;
        this->srefName = srefName;
        this->transform = transform;
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

    // Get name of referred geometric cell
    std::string getSRefName() const
    {
        return this->srefName;
    }

    // Get the linear transformations (reflection, rotation, and scaling) of the structure reference
    strans getTransform() const
    {
        return this->transform;
    }

    // Set placed structure reference center in cell coordinates
    void setSRefs(vector<double> srefs)
    {
        if (srefs.size() != 2)
        {
            cerr << "Structure reference center coordinates must be stored in vector of length 2. Defaulting to origin (0, 0)." << endl;
            this->srefs = { 0., 0. };
        }
        else
        {
            this->srefs = srefs;
        }
    }

    // Set structure reference properties
    void setProps(vector<std::string> props)
    {
        this->props = props;
    }

    // Set name of referred geometric cell
    void setSRefName(std::string srefName)
    {
        this->srefName = srefName;
    }

    // Set linear transformations (reflection, rotation, and scaling) of the structure reference
    void setTransform(strans transform)
    {
        this->transform = transform;
    }

    // Destructor
    ~sref()
    {
        // Nothing
    }
};

class aref
{
private:
    vector<double> arefs;          // Three coordinate pairs of arrayed structure references (corner instance center, last columnar instance, last row instance)
    vector<std::string> props;     // Properties of array reference
    std::string arefName;          // Name of arrayed structure reference
    vector<int> numColRow;         // Number of columns and rows in array
    strans transform;              // Linear transformations (reflection, rotation, and scaling) of arrayed reference instances
public:
    // Default constructor
    aref()
    {
        vector<double> arefs = { 0., 0., 0., 0., 0., 0. };
        vector<int> numColRow = { 0, 0 };
        strans transform;
        this->arefs = arefs;
        this->props = {};
        this->arefName = "";
        this->numColRow = numColRow;
        this->transform = transform;
    }

    // Parametrized constructor
    aref(vector<double> arefs, vector<std::string> props, std::string arefName, vector<int> numColRow, strans transform)
    {
        if (arefs.size() != 6)
        {
            cerr << "Array reference center coordinates must be stored in vector of length 6 for corner instance, last columnar instance, and last row instance. Defaulting to origin (0, 0) for all." << endl;
            this->arefs = { 0., 0., 0., 0., 0., 0. };
        }
        else
        {
            this->arefs = arefs;
        }
        this->props = props;
        this->arefName = arefName;
        if (numColRow.size() != 2)
        {
            cerr << "Array reference must have vector of 2 positive integers describing the number of columns and the number of rows in the array, respectively. Defaulting to 2 columns and 2 rows." << endl;
            this->numColRow = { 2, 2 };
        }
        else if ((numColRow[0] < 0) || (numColRow[1] < 0))
        {
            cerr << "Array reference must have a positive number of columns and a positive number of rows in the array. Defaulting to absolute value of given integers." << endl;
            this->numColRow = { abs(numColRow[0]), abs(numColRow[1]) };
        }
        else
        {
            this->numColRow = numColRow;
        }
        this->transform = transform;
    }

    // Get coordinates of the arrayed structure references (corner instance center, last columnar instance, last row instance)
    vector<double> getARefs() const
    {
        return this->arefs;
    }

    // Get array reference properties
    vector<std::string> getProps() const
    {
        return this->props;
    }

    // Get name of arrayed, referred geometric cell
    std::string getARefName() const
    {
        return this->arefName;
    }

    // Get number of columns and rows in array as size-2 vector of integers
    vector<int> getNumColRow() const
    {
        return this->numColRow;
    }

    // Get the linear transformations (reflection, rotation, and scaling) of each arrayed reference instance
    strans getTransform() const
    {
        return this->transform;
    }

    // Set coordinates of the arrayed structure references (corner instance center, last columnar instance, last row instance)
    void setARefs(vector<double> arefs)
    {
        if (arefs.size() != 6)
        {
            cerr << "Array reference center coordinates must be stored in vector of length 6 for corner instance, last columnar instance, and last row instance. Defaulting to origin (0, 0) for all." << endl;
            this->arefs = { 0., 0., 0., 0., 0., 0. };
        }
        else
        {
            this->arefs = arefs;
        }
    }

    // Set array reference properties
    void setProps(vector<std::string> props)
    {
        this->props = props;
    }

    // Set name of arrayed, referred geometric cell
    void setARefName(std::string arefName)
    {
        this->arefName = arefName;
    }

    // Set number of columns and rows in array as size-2 vector of integers
    void setNumColRow(vector<int> numColRow)
    {
        if (numColRow.size() != 2)
        {
            cerr << "Array reference must have vector of 2 positive integers describing the number of columns and the number of rows in the array, respectively. Defaulting to 2 columns and 2 rows." << endl;
            this->numColRow = { 2, 2 };
        }
        else if ((numColRow[0] < 0) || (numColRow[1] < 0))
        {
            cerr << "Array reference must have a positive number of columns and a positive number of rows in the array. Defaulting to absolute value of given integers." << endl;
            this->numColRow = { abs(numColRow[0]), abs(numColRow[1]) };
        }
        else
        {
            this->numColRow = numColRow;
        }
    }

    // Set linear transformations (reflection, rotation, and scaling) of the array reference instances
    void setTransform(strans transform)
    {
        this->transform = transform;
    }

    // Find columnar separation vector
    vector<double> findColSep() const
    {
        double numCol = (double)(this->numColRow[0]); // Number of columns in array
        if (numCol >= 2.) // Normal case
        {
            double xColSep = ((this->arefs)[2] - (this->arefs[0])) / (numCol - 1.); // x-component of separation between columnar instances
            double yColSep = ((this->arefs)[3] - (this->arefs[1])) / (numCol - 1.); // y-component of separation between columnar instances
            return{ xColSep, yColSep };
        }
        else // Separation undefined if singleton or nonexistant column
        {
            return{ 0., 0. };
        }
    }

    // Find row separation vector
    vector<double> findRowSep() const
    {
        double numRow = (double)(this->numColRow[1]); // Number of rows in array
        if (numRow >= 2.) // Normal case
        {
            double xRowSep = ((this->arefs)[4] - (this->arefs[0])) / (numRow - 1.); // x-component of separation between row instances
            double yRowSep = ((this->arefs)[5] - (this->arefs[1])) / (numRow - 1.); // y-component of separation between row instances
            return{ xRowSep, yRowSep };
        }
        else // Separation undefined if singleton or nonexistant row
        {
            return{ 0., 0. };
        }
    }

    // Find center coordinates of all referenced instances in array with translation
    vector<vector<double>> findInstances(vector<double> center) const
    {
        double xo, yo; // Coordinate offsets
        if (center.size() != 2)
        {
            cerr << "Coordinates of reference frame center must be a length-2 vector. Defaulting to (0, 0)." << endl;
            xo = 0.;
            yo = 0.;
        }
        else
        {
            xo = center[0];
            yo = center[1];
        }
        vector<vector<double>> instanceCoord; // Vector of all center (x,y)-coordinates of arrayed instances
        int numCol = (this->getNumColRow())[0]; // Number of columns
        int numRow = (this->getNumColRow())[1]; // Number of rows
        vector<double> colSep = this->findColSep(); // Column separation vector
        vector<double> rowSep = this->findRowSep(); // Row separation vector
        for (size_t indi = 0; indi < numCol - 1; indi++)
        {
            for (size_t indj = 0; indj < numRow - 1; indj++)
            {
                instanceCoord.push_back({ this->arefs[0] + indi * colSep[0] + indj * rowSep[0] + xo, this->arefs[1] + indi * colSep[0] + indj * rowSep[0] + yo }); // Calculate all points from corner
            }
        }
        return instanceCoord;
    }

    // Destructor
    ~aref()
    {
        // Nothing
    }
};

// Custom class for GDSII geometric cell information
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
    vector<aref> areferences;                    // All arrayed structure references in geometric cell
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
        vector<aref> areferences;
        this->cellName = cellName;
        this->dateCreate = dateCreate;
        this->dateMod = dateMod;
        this->boundaries = boundaries;
        this->paths = paths;
        this->nodes = nodes;
        this->boxes = boxes;
        this->textboxes = textboxes;
        this->sreferences = sreferences;
        this->areferences = areferences;
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

    // Get vector of array references
    vector<aref> getARefVec() const
    {
        return this->areferences;
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

    // Return number of array references
    int getNumARef() const
    {
        return (this->areferences).size();
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
        int numARef = GeoCell::getNumARef();

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
            (this->textboxes)[indi].getTransform().print();
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
            (this->sreferences)[indi].getTransform().print();
            //delete[] point; // Free array pointer
        }
        cout << "  List of " << numARef << " array references:" << endl;
        for (size_t indi = 0; indi < numARef; indi++) // Handle each array reference
        {
            vector<int> numColRow = ((this->areferences)[indi]).getNumColRow();
            if (((this->areferences)[indi]).getProps().size() != 0) // Check if property vector has nonzero length
            {
                cout << "   #" << indi + 1 << " using " << ((this->areferences)[indi]).getARefName() << " in " << numColRow[0] << " columns and " << numColRow[1] << " rows with " << (((this->areferences)[indi]).getProps())[0] << " as Property " << 1 << endl;
            }
            else
            {
                cout << "   #" << indi + 1 << " using " << ((this->areferences)[indi]).getARefName() << " in " << numColRow[0] << " columns and " << numColRow[1] << " rows" << endl;
            }
            char point[256];
            vector<double> arefCoord = ((this->areferences)[indi]).getARefs();
            sprintf(point, "corner instance: (%1.4g, %1.4g), last columnar: (%1.4g, %1.4g), last row: (%1.4g, %1.4g)", arefCoord[0], arefCoord[1], arefCoord[2], arefCoord[3], arefCoord[4], arefCoord[5]);
            cout << "    " << point << endl; // Print ordered pair of center
            (this->areferences)[indi].getTransform().print();
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
        int numARef = GeoCell::getNumARef();

        cout << " ------" << endl;
        cout << " Cell Details:" << endl;
        cout << "  Cell Name: " << this->cellName << endl;
        cout << "  List of " << numBound << " boundaries [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numBound; indi++) // Handle each boundary
        {
            cout << "   " << ((this->boundaries)[indi]).getLayer() << " ";
            vector<double> boundCoord = ((this->boundaries)[indi]).getBounds();
            for (size_t indj = 0; indj < boundCoord.size(); indj += 2)
            {
                cout << boundCoord[indj] << " " << boundCoord[indj + 1] << " ";
            }
            cout << ((this->boundaries)[indi]).getNBoundPt() << endl;
        }
        cout << "  List of " << numPath << " paths [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numPath; indi++) // Handle each path
        {
            cout << "   " << ((this->paths)[indi]).getLayer() << " ";
            vector<double> pathCoord = ((this->paths)[indi]).getPaths();
            for (size_t indj = 0; indj < pathCoord.size(); indj += 2)
            {
                cout << pathCoord[indj] << " " << pathCoord[indj + 1] << " ";
            }
            cout << ((this->paths)[indi]).getNPathPt() << endl;
        }
        cout << "  List of " << numNode << " nodes [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numNode; indi++) // Handle each node
        {
            cout << "   " << ((this->nodes)[indi]).getLayer() << " ";
            vector<double> nodeCoord = ((this->nodes)[indi]).getNodes();
            for (size_t indj = 0; indj < nodeCoord.size(); indj += 2) // C string for each ordered pair
            {
                cout << nodeCoord[indj] << " " << nodeCoord[indj + 1] << " ";
            }
            cout << ((this->nodes)[indi]).getNNodePt() << endl;
        }
        cout << "  List of " << numBox << " box outlines [layer <coord> numPt]:" << endl;
        for (size_t indi = 0; indi < numBox; indi++) // Handle each box outline
        {
            cout << "   " << ((this->boxes)[indi]).getLayer() << " ";
            vector<double> boxCoord = ((this->boxes)[indi]).getBoxes();
            for (size_t indj = 0; indj < boxCoord.size(); indj += 2) // C string for each ordered pair (should be 5, sometimes 4)
            {
                cout << boxCoord[indj] << " " << boxCoord[indj + 1] << " ";
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
        cout << "  List of " << numARef << " array references [name (x-corner, y-corner) (x-columnar, y-columnar) (x-row, y-row)]:" << endl;
        for (size_t indi = 0; indi < numARef; indi++) // Handle each structure reference
        {
            cout << "   " << ((this->areferences)[indi]).getARefName() << " ";
            char point[256];
            vector<double> arefCoord = ((this->areferences)[indi]).getARefs();
            sprintf(point, "(%1.4g, %1.4g) (%1.4g, %1.4g) (%1.4g, %1.4g)", arefCoord[0], arefCoord[1], arefCoord[2], arefCoord[3], arefCoord[4], arefCoord[5]);
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
        ~(this->sreferences);
        ~(this->areferences)*/
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
    vector<std::string> findNames() const
    {
        vector<std::string> names;
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
        // Initialize variables
        double xo, yo; // Coordinate offsets
        if (center.size() != 2)
        {
            cerr << "Coordinates of reference frame center must be a length-2 vector. Defaulting to (0, 0)." << endl;
            xo = 0.;
            yo = 0.;
        }
        else
        {
            xo = center[0];
            yo = center[1];
        }
        GeoCell thisCell = this->getCell(indCell);
        size_t indVia = this->locateCell(viaName);
        vector<vector<double>> viaList;

        // Handle each structure reference
        for (size_t indi = 0; indi < thisCell.getNumSRef(); indi++)
        {
            string srefName = (thisCell.sreferences)[indi].getSRefName();
            size_t indNextCell = this->locateCell(srefName);
            vector<double> thisSRef = (thisCell.sreferences)[indi].getSRefs();
            if (indNextCell == indVia) // Cell is a via ...
            {
                viaList.push_back({ thisSRef[0] + xo, thisSRef[1] + yo });
            }
            else if (indNextCell < this->getNumCell()) // Valid, non-via cell name was found ...
            {
                vector<vector<double>> newVias = this->findVias(indNextCell, { thisSRef[0] + xo, thisSRef[1] + yo }, viaName); // Recursion step
                viaList.insert(viaList.end(), newVias.begin(), newVias.end());
            }
        }

        // Handle each array reference
        for (size_t indi = 0; indi < thisCell.getNumARef(); indi++)
        {
            string arefName = (thisCell.areferences)[indi].getARefName();
            size_t indNextCell = this->locateCell(arefName);
            vector<vector<double>> instanceCoord = (thisCell.areferences)[indi].findInstances({ xo, yo });
            for (size_t indj = 0; indj < instanceCoord.size(); indj++)
            {
                if (indNextCell == indVia) // Cell is a via ...
                {
                    viaList.push_back({ instanceCoord[indj][0], instanceCoord[indj][1] });
                }
                else if (indNextCell < this->getNumCell()) // Valid, non-via cell name was found ...
                {
                    vector<vector<double>> newVias = this->findVias(indNextCell, { instanceCoord[indj][0], instanceCoord[indj][1] }, viaName); // Recursion step
                    viaList.insert(viaList.end(), newVias.begin(), newVias.end());
                }
            }
        }

        // Return the via list
        return viaList;
    }

    // Return all physical points for entire geometric cell
    vector<complex<double>> findPoints(std::string cellName, vector<double> center, strans transform)
    {
        // Error checking on input point
        double xo, yo; // Coordinate offsets
        if (center.size() != 2)
        {
            cerr << "Coordinates of reference frame center must be a length-2 vector. Defaulting to (0, 0)." << endl;
            xo = 0.;
            yo = 0.;
        }
        else
        {
            xo = center[0];
            yo = center[1];
        }

        // Create vector of all physical points in named cell for all layers
        const GeoCell cell = this->cells[this->locateCell(cellName)];
        vector<complex<double>> allPt; // Using complex numbers to represent ordered pairs
        for (size_t indi = 0; indi < cell.getNumBound(); indi++) // Handle each boundary
        {
            vector<double> boundCoord = (cell.boundaries[indi]).getBounds();
            for (size_t indj = 0; indj < boundCoord.size() / 2 - 1; indj++) // Handle each coordinate comprising boundary except repeated point
            {
                vector<double> boundPt = transform.applyTranform({ boundCoord[2 * indj], boundCoord[2 * indj + 1] }); // Apply the linear transformation of this cell reference
                complex<double> zBound(boundPt[0] + xo, boundPt[1] + yo);
                allPt.push_back(zBound); // Complex number with x- and y-coordinates (m)
            }
        }
        for (size_t indi = 0; indi < cell.getNumPath(); indi++) // Handle each path
        {
            vector<double> pathCoord = (cell.paths[indi]).getPaths();
            for (size_t indj = 0; indj < pathCoord.size() / 2; indj++) // Handle each coordinate comprising path
            {
                vector<double> pathPt = transform.applyTranform({ pathCoord[2 * indj], pathCoord[2 * indj + 1] }); // Apply the linear transformation of this cell reference
                complex<double> zPath(pathPt[0] + xo, pathPt[1] + yo);
                allPt.push_back(zPath); // Complex number with x- and y-coordinates (m)
            }
        }
        for (size_t indi = 0; indi < cell.getNumBox(); indi++) // Handle each box outline
        {
            vector<double> boxCoord = (cell.boxes[indi]).getBoxes();
            for (size_t indj = 0; indj < boxCoord.size() / 2 - 1; indj++) // Handle each coordinate comprising box outline except repeated point
            {
                vector<double> boxPt = transform.applyTranform({ boxCoord[2 * indj], boxCoord[2 * indj + 1] }); // Apply the linear transformation of this cell reference
                complex<double> zBox(boxPt[0] + xo, boxPt[1] + yo);
                allPt.push_back(zBox); // Complex number with x- and y-coordinates (m)
            }
        }
        for (size_t indi = 0; indi < cell.getNumSRef(); indi++) // Handle each structure reference recursively
        {
            vector<double> refPt = transform.applyTranform(((cell.sreferences)[indi]).getSRefs()); // Apply the linear transformation of this cell reference
            vector<complex<double>> newCoord = this->findPoints((cell.sreferences)[indi].getSRefName(), { refPt[0] + xo, refPt[1] + yo }, (cell.sreferences)[indi].getTransform().composeTransform(transform));
            allPt.insert(allPt.end(), newCoord.begin(), newCoord.end());
        }
        for (size_t indi = 0; indi < cell.getNumARef(); indi++) // Handle each array reference instance
        {
            vector<vector<double>> instanceCoord = (cell.areferences[indi]).findInstances({ 0., 0. });
            for (size_t indj = 0; indj < instanceCoord.size(); indj++) // Handle each instance in array reference recursively
            {
                vector<double> centPt = transform.applyTranform(instanceCoord[indj]); // Apply the linear transformation of this cell reference
                vector<complex<double>> newCoord = this->findPoints((cell.areferences)[indi].getARefName(), { centPt[0] + xo, centPt[1] + yo }, (cell.areferences)[indi].getTransform().composeTransform(transform));
                allPt.insert(allPt.end(), newCoord.begin(), newCoord.end());
            }
        }

        // Return the points
        allPt.erase(unique(allPt.begin(), allPt.end()), allPt.end()); // Removes only consecutive repeated points (those corresponding to zero-length changes)
        return allPt;
    }

    // Return all physical points, connecting segments, and regions (PSLG information) on given layer within a geometric cell
    pslg findPSLG(std::string name, int layer, vector<double> center, strans transform)
    {
        // Error checking on input point
        double xo, yo; // Coordinate offsets
        if (center.size() != 2)
        {
            cerr << "Coordinates of reference frame center must be a length-2 vector. Defaulting to (0, 0)." << endl;
            xo = 0.;
            yo = 0.;
        }
        else
        {
            xo = center[0];
            yo = center[1];
        }

        // Store all physical vertices for layer in vector of complex numbers
        const GeoCell cell = this->cells[this->locateCell(name)];
        vector<complex<double>> layerPt; // Using complex numbers to represent ordered pairs for vertices
        vector<pair<size_t, size_t>> layerSeg; // Using pairs of indices to represent the points connected by segments
        vector<complex<double>> layerReg; // Using complex numbers to represent an interior point of conductor regions
        //unordered_map<complex<double>, long long int, hash<complex<double>>> backwards; // Backwards map to ensure unique points entered (fails because complex numbers have no hash)
        //map<long long int, complex<double>> layerPt; // Forwards map has enumeration of unique points
        //long long int indPt = 0;
        for (size_t indi = 0; indi < cell.getNumBound(); indi++) // Handle each boundary (limitation: no GDSII element may overlap)
        {
            if ((cell.boundaries[indi]).getLayer() == layer)
            {
                vector<double> boundCoord = (cell.boundaries[indi]).getBounds();
                size_t numBoundPt = (cell.boundaries[indi]).getNBoundPt() - 1; // Number of unique points on boundary
                size_t indBoundStart = layerPt.size(); // Index of first point on boundary
                for (size_t indj = 0; indj < numBoundPt; indj++) // Handle each coordinate comprising boundary except repeated point
                {
                    vector<double> boundPt = transform.applyTranform({ boundCoord[2 * indj], boundCoord[2 * indj + 1] }); // Apply the linear transformation of this cell reference
                    complex<double> zBound(boundPt[0] + xo, boundPt[1] + yo); // Complex number with x- and y-coordinates (m)
                    layerPt.push_back(zBound);
                    layerSeg.push_back(pair<size_t, size_t>(indBoundStart + indj, indBoundStart + ((indj + 1) % numBoundPt))); // Segment connects this point to next cyclically
                    //bool couldInsert = backwards.insert(pair<complex<double>, long long int>(zBound, indPt)).second; // Check if point already exists
                    //if (couldInsert)
                    //{
                    //    layerPt.insert(pair<long long int, complex<double>>(indPt, zBound)); // Enumerate the unique point
                    //    indPt++;
                    //}
                }
                //vector<double> centroidCoord = (cell.boundaries[indi]).findCentroid(xo, yo); // Boundary centroid as {xcent, ycent} (m)
                complex<double> boundInterior = { boundCoord[0] - REL_INTERIOR_PRECISION * fabs(boundCoord[0] - boundCoord[2 * numBoundPt - 2]) + xo, boundCoord[1] + REL_INTERIOR_PRECISION * fabs(boundCoord[1] - boundCoord[3]) + yo }; // Point just left and above the bottom-right point (small offset based on nearest points)
                layerReg.push_back(boundInterior); // New region for this polygon
            }
        }
        for (size_t indi = 0; indi < cell.getNumPath(); indi++) // Handle each path (does not work because implementation does not work around overlapping points)
        {
            if ((cell.paths[indi]).getLayer() == layer)
            {
                vector<double> pathCoord = (cell.paths[indi]).getPaths();
                double width = (cell.paths[indi]).getWidth();
                int type = (cell.paths[indi]).getType();
                size_t numPathPt = (cell.paths[indi]).getNPathPt(); // Number of points along path (these are the path points, not the points outlining a finite-width path)
                size_t indPathStart = layerPt.size();
                //vector<complex<double>> returnPt; // Vector of path points needed to close the figure

                // Handle the points along a path of finite width
                for (size_t indj = 0; indj < numPathPt - 1; indj++) // Handle each coordinate comprising path
                {
                    double overshoot1 = 1.0; // Assume path has round ends (type = 1) approximated as square or square ends with overshoot (type = 2)
                    double overshoot2 = 1.0; // Assume path has round ends (type = 1) approximated as square or square ends with overshoot (type = 2)
                    if ((type == 0) && (indj == 0))
                    {
                        overshoot1 = 0.0; // Path actually has square ends at terminal vertices
                    }
                    if ((type == 0) && (indj == numPathPt - 1))
                    {
                        overshoot2 = 0.0; // Path actually has square ends at terminal vertices
                    }
                    if (pathCoord[2 * indj] == pathCoord[2 * indj + 2]) // Path segment along y-axis
                    {
                        if (pathCoord[2 * indj + 1] > pathCoord[2 * indj + 3]) // First point is above the second point
                        {
                            vector<double> pathPt = transform.applyTranform({ pathCoord[2 * indj] - width / 2., pathCoord[2 * indj + 1] + overshoot1 * width / 2. }); // Apply the linear transformation of this cell reference
                            complex<double> zPath(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-left point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj] + width / 2., pathCoord[2 * indj + 1] + overshoot1 * width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-right point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] + width / 2., pathCoord[2 * indj + 3] - overshoot2 * width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-right point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] - width / 2., pathCoord[2 * indj + 3] - overshoot2 * width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-left point of path segment
                            layerPt.push_back(zPath);
                            vector<double> centroidPt = transform.applyTranform({ pathCoord[2 * indj], (pathCoord[2 * indj + 1] + pathCoord[2 * indj + 3]) / 2. + (+overshoot1 - overshoot2) * width / 4. }); // Apply the linear transformation
                            complex<double> centroidCoord(centroidPt[0] + xo, centroidPt[1] + yo);
                            layerReg.push_back(centroidCoord); // New region for this path segment
                        }
                        else // Second point is above the first point
                        {
                            vector<double> pathPt = transform.applyTranform({ pathCoord[2 * indj] - width / 2., pathCoord[2 * indj + 1] - overshoot1 * width / 2. }); // Apply the linear transformation of this cell reference
                            complex<double> zPath(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-left point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj] + width / 2., pathCoord[2 * indj + 1] - overshoot1 * width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-right point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] + width / 2., pathCoord[2 * indj + 3] + overshoot2 * width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-right point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] - width / 2., pathCoord[2 * indj + 3] + overshoot2 * width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-left point of path segment
                            layerPt.push_back(zPath);
                            vector<double> centroidPt = transform.applyTranform({ pathCoord[2 * indj], (pathCoord[2 * indj + 1] + pathCoord[2 * indj + 3]) / 2. + (-overshoot1 + overshoot2) * width / 4. }); // Apply the linear transformation
                            complex<double> centroidCoord(centroidPt[0] + xo, centroidPt[1] + yo);
                            layerReg.push_back(centroidCoord); // New region for this path segment
                        }
                    }
                    else // Path segment along x-axis
                    {
                        if (pathCoord[indj] > pathCoord[indj + 2]) // First point is on the right of the second point
                        {
                            vector<double> pathPt = transform.applyTranform({ pathCoord[2 * indj] + overshoot1 * width / 2., pathCoord[2 * indj + 1] + width / 2. }); // Apply the linear transformation of this cell reference
                            complex<double> zPath(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-right point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj] + overshoot1 * width / 2., pathCoord[2 * indj + 1] - width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-right point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] - overshoot2 * width / 2., pathCoord[2 * indj + 3] - width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-left point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] - overshoot2 * width / 2., pathCoord[2 * indj + 3] + width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-left point of path segment
                            layerPt.push_back(zPath);
                            vector<double> centroidPt = transform.applyTranform({ (pathCoord[2 * indj] + pathCoord[2 * indj + 2]) / 2. + (+overshoot1 - overshoot2) * width / 4., pathCoord[2 * indj + 1] }); // Apply the linear transformation
                            complex<double> centroidCoord(centroidPt[0] + xo, centroidPt[1] + yo);
                            layerReg.push_back(centroidCoord); // New region for this path segment
                        }
                        else // Second point is on the right of the first point
                        {
                            vector<double> pathPt = transform.applyTranform({ pathCoord[2 * indj] - overshoot1 * width / 2., pathCoord[2 * indj + 1] + width / 2. }); // Apply the linear transformation of this cell reference
                            complex<double> zPath(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-left point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj] - overshoot1 * width / 2., pathCoord[2 * indj + 1] - width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-left point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] + overshoot2 * width / 2., pathCoord[2 * indj + 3] - width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of lower-right point of path segment
                            layerPt.push_back(zPath);
                            pathPt = transform.applyTranform({ pathCoord[2 * indj + 2] + overshoot2 * width / 2., pathCoord[2 * indj + 3] + width / 2. });
                            zPath = complex<double>(pathPt[0] + xo, pathPt[1] + yo); // Complex number of upper-right point of path segment
                            layerPt.push_back(zPath);
                            vector<double> centroidPt = transform.applyTranform({ (pathCoord[2 * indj] + pathCoord[2 * indj + 2]) / 2. + (-overshoot1 + overshoot2) * width / 4., pathCoord[2 * indj + 1] }); // Apply the linear transformation
                            complex<double> centroidCoord(centroidPt[0] + xo, centroidPt[1] + yo);
                            layerReg.push_back(centroidCoord); // New region for this path segment
                        }
                    }

                    // Handle the segments connecting the points on each path segment
                    layerSeg.push_back(pair<size_t, size_t>(indPathStart + 4 * indj + 0, indPathStart + 4 * indj + 1)); // Connect first and second points in path segment
                    layerSeg.push_back(pair<size_t, size_t>(indPathStart + 4 * indj + 1, indPathStart + 4 * indj + 2)); // Connect second and third points in path segment
                    layerSeg.push_back(pair<size_t, size_t>(indPathStart + 4 * indj + 2, indPathStart + 4 * indj + 3)); // Connect third and fourth points in path segment
                    layerSeg.push_back(pair<size_t, size_t>(indPathStart + 4 * indj + 3, indPathStart + 4 * indj + 0)); // Connect fourth and first points in path segment

                    //bool couldInsert = backwards.insert(pair<complex<double>, long long int>(zPath, indPt)).second; // Check if point already exists
                    //if (couldInsert)
                    //{
                    //    layerPt.insert(pair<long long int, complex<double>>(indPt, zPath)); // Enumerate the unique point
                    //    indPt++;
                    //}
                }
            }
        }
        for (size_t indi = 0; indi < cell.getNumBox(); indi++) // Handle each box outline
        {
            if ((cell.boxes[indi]).getLayer() == layer)
            {
                vector<double> boxCoord = (cell.boxes[indi]).getBoxes();
                size_t numBoxPt = (cell.boxes[indi]).getNBoxPt() - 1; // Number of unique points on box outline
                size_t indBoxStart = layerPt.size(); // Index of first point on box outline
                for (size_t indj = 0; indj < numBoxPt; indj++) // Handle each coordinate comprising box outline except repeated point
                {
                    vector<double> boxPt = transform.applyTranform({ boxCoord[2 * indj], boxCoord[2 * indj + 1] }); // Apply the linear transformation of this cell reference
                    complex<double> zBox(boxPt[0] + xo, boxPt[1] + yo); // Complex number with x- and y-coordinates (m)
                    layerPt.push_back(zBox);
                    layerSeg.push_back(pair<size_t, size_t>(indBoxStart + indj, indBoxStart + ((indj + 1) % numBoxPt))); // Segment connects this point to next cyclically
                    //bool couldInsert = backwards.insert(pair<complex<double>, long long int>(zBox, indPt)).second; // Check if point already exists
                    //if (couldInsert)
                    //{
                    //    layerPt.insert(pair<long long int, complex<double>>(indPt, zBox)); // Enumerate the unique point
                    //    indPt++;
                    //}
                }
                vector<double> centroidCoord = (cell.boxes[indi]).findCentroid(xo, yo);
                complex<double> centroid(centroidCoord[0], centroidCoord[1]); // Box outline centroid (m, xcent + i*ycent)
                layerReg.push_back(centroid); // New region for this box
            }
        }
        for (size_t indi = 0; indi < cell.getNumSRef(); indi++) // Handle each structure reference recursively
        {
            vector<double> refPt = transform.applyTranform(((cell.sreferences)[indi]).getSRefs()); // Apply the linear transformation of this cell reference
            pslg newPSLG = this->findPSLG((cell.sreferences)[indi].getSRefName(), layer, { refPt[0] + xo, refPt[1] + yo }, (cell.sreferences)[indi].getTransform().composeTransform(transform));
            vector<complex<double>> newVertices = newPSLG.vertices;
            vector<pair<size_t, size_t>> newSegments = newPSLG.segments; // Segments have zero-referenced indices for referenced geometric cell
            vector<complex<double>> newRegions = newPSLG.regions;
            size_t refCellIndOffset = layerPt.size();
            for (size_t indj = 0; indj < newSegments.size(); indj++)
            {
                newSegments[indj] = pair<size_t, size_t>(newSegments[indj].first + refCellIndOffset, newSegments[indj].second + refCellIndOffset); // Add offset to all values in segment pairs
            }
            layerPt.insert(layerPt.end(), newVertices.begin(), newVertices.end());
            layerSeg.insert(layerSeg.end(), newSegments.begin(), newSegments.end());
            layerReg.insert(layerReg.end(), newRegions.begin(), newRegions.end());
            //map<long long int, complex<double>> newPt = this->findLayerPt((cell.sreferences)[indi].getSRefName(), layer, (((cell.sreferences)[indi]).getSRefs())[0] + xo, (((cell.sreferences)[indi]).getSRefs())[1] + yo);
            //for (auto it = newPt.begin(); it != newPt.end(); ++it)
            //{
            //    layerPt.insert(pair<long long int, complex<double>>(indPt, it->second)); // Include the unique points from the structure reference
            //    indPt++;
            //}
        }
        for (size_t indi = 0; indi < cell.getNumARef(); indi++) // Handle each array reference instance
        {
            vector<vector<double>> instanceCoord = (cell.areferences[indi]).findInstances({ 0., 0. });
            for (size_t indj = 0; indj < instanceCoord.size(); indj++) // Handle each instance in array reference recursively
            {
                vector<double> centPt = transform.applyTranform(instanceCoord[indj]); // Apply the linear transformation of this cell reference
                pslg newPSLG = this->findPSLG((cell.areferences)[indi].getARefName(), layer, { centPt[0] + xo, centPt[1] + yo }, (cell.areferences)[indi].getTransform().composeTransform(transform));
                vector<complex<double>> newVertices = newPSLG.vertices;
                vector<pair<size_t, size_t>> newSegments = newPSLG.segments; // Segments have zero-referenced indices for referenced geometric cell
                vector<complex<double>> newRegions = newPSLG.regions;
                size_t refCellIndOffset = layerPt.size();
                for (size_t indk = 0; indk < newSegments.size(); indk++)
                {
                    newSegments[indk] = pair<size_t, size_t>(newSegments[indk].first + refCellIndOffset, newSegments[indk].second + refCellIndOffset); // Add offset to all values in segment pairs
                }
                layerPt.insert(layerPt.end(), newVertices.begin(), newVertices.end());
                layerSeg.insert(layerSeg.end(), newSegments.begin(), newSegments.end());
                layerReg.insert(layerReg.end(), newRegions.begin(), newRegions.end());
                //map<long long int, complex<double>> newPt = this->findLayerPt((cell.sreferences)[indi].getSRefName(), layer, (((cell.sreferences)[indi]).getSRefs())[0] + xo, (((cell.sreferences)[indi]).getSRefs())[1] + yo);
                //for (auto it = newPt.begin(); it != newPt.end(); ++it)
                //{
                //    layerPt.insert(pair<long long int, complex<double>>(indPt, it->second)); // Include the unique points from the structure reference
                //    indPt++;
                //}
            }
        }

        // Return structure
        return pslg(layerPt, layerSeg, layerReg);
    }

    // Compute the convex hull of points in geometric cell with Graham scan
    vector<complex<double>> convexHull(std::string cellName)
    {
        // Retrieve vector of all physical points in named cell for all layers
        vector<complex<double>> allPt = this->findPoints(cellName, { 0., 0. }, strans());
        cout << "Identified " << allPt.size() << " points along PCB outline" << endl;

        // Find bottommost point and put at front
        size_t indPtMin = 0;
        double yMin = allPt[0].imag(); // Initialize to y-coordinate of first point
        for (size_t indi = 1; indi < allPt.size(); indi++)
        {
            // Tie breaker is leftmost point if compared points have same y-coordinate
            if ((allPt[indi].imag() < yMin) || ((allPt[indi].imag() == yMin) && (allPt[indi].real() < allPt[indPtMin].real())))
            {
                indPtMin = indi;
                yMin = allPt[indi].imag();
            }
        }
        iter_swap(allPt.begin(), allPt.begin() + indPtMin); // Bottommost point moved to front

        // Sort points relative to first point using custom structure
        convexCompare.p0 = allPt[0];
        sort(allPt.begin() + 1, allPt.end(), convexCompare);

        // Prepare vector of convex hull point candidates while removing those at same angle to first point
        vector<complex<double>> candPt;
        candPt.push_back(allPt[0]); // Include bottommost point
        for (size_t indi = 1; indi < allPt.size(); indi++)
        {
            int orientation = 0; // Assume collinear points

            // Points with same orientation angle are not added as candidates
            while ((indi < allPt.size() - 1) && (orientation == 0))
            {
                // Find orientation angle with cross product of displacements
                double crossProd = (allPt[indi].imag() - allPt[0].imag()) * (allPt[indi + 1].real() - allPt[indi].real()) - (allPt[indi].real() - allPt[0].real()) * (allPt[indi + 1].imag() - allPt[indi].imag());
                if (crossProd > 0.)
                {
                    orientation = 1; // Clockwise turn
                    break;
                }
                else if (crossProd < 0.)
                {
                    orientation = 2; // Counterclockwise turn
                    break;
                }
                else
                {
                    orientation = 0; // Still collinear
                }
                indi++;
            }

            // All remaining points added as candidates
            candPt.push_back(allPt[indi]);
        }

        // Convex hull not possible if fewer than 3 points are candidates
        if (candPt.size() < 3)
        {
            cerr << "Not enough points remain to construct a convex hull of PCB outline" << endl;
            return{};
        }

        // Convex hull will start with first three candidate points
        vector<complex<double>> hullPt = { candPt[0], candPt[1], candPt[2] };

        // Convex hull points change based on orientation angle of remaining candidates
        for (size_t indi = 3; indi < candPt.size(); indi++)
        {
            int orientation = 0; // Assume collinear

            // Points with that make a non-left turn need to be removed
            while (orientation != 2)
            {
                // Find orientation angle with cross product of displacements of next-to-last, last, and next candidate points
                double crossProd = (hullPt[hullPt.size() - 1].imag() - hullPt[hullPt.size() - 2].imag()) * (candPt[indi].real() - hullPt[hullPt.size() - 1].real()) - (hullPt[hullPt.size() - 1].real() - hullPt[hullPt.size() - 2].real()) * (candPt[indi].imag() - hullPt[hullPt.size() - 1].imag());
                if (crossProd > 0.)
                {
                    orientation = 1; // Clockwise turn
                    hullPt.pop_back();
                }
                else if (crossProd < 0.)
                {
                    orientation = 2; // Counterclockwise turn
                    break;
                }
                else
                {
                    orientation = 0; // Collinear points
                    hullPt.pop_back();
                }
            }

            // Add the current point as the last convex hull point for now
            hullPt.push_back(candPt[indi]);
        }

        // Return the valid list of convex hull points as a complex vector
        cout << "Convex hull of PCB outline is ready" << endl;
        return hullPt;
    }

    // Produce planar straight-line graph (PSLG) file for layer within geometric cell
    bool convertPSLG(std::string cellName, int layer, vector<complex<double>> outlinePt)
    {
        // Attempt to open planar straight-line graph (PSLG) file for writing
        std::string polyName = this->fileName.substr(0, this->fileName.find_last_of(".")) + "_" + to_string(layer) + ".poly";
        ofstream polyFile(polyName.c_str());
        if (polyFile.is_open())
        {
            // Retrieve list of points, segments, and regions on layer
            pslg layerPSLG = this->findPSLG(cellName, layer, { 0., 0. }, strans());
            vector<complex<double>> layerPt = layerPSLG.vertices;
            vector<pair<size_t, size_t>> layerSeg = layerPSLG.segments;
            vector<complex<double>> layerReg = layerPSLG.regions;
            size_t numOutlinePt = outlinePt.size();
            size_t numLayerPt = layerPt.size();
            size_t numLayerSeg = layerSeg.size();
            size_t numLayerReg = layerReg.size();

            // Create fictitious node of board outline to guess interior point (not robust)
            vector<double> outlineCoord;
            for (size_t indi = 0; indi < outlinePt.size(); indi++)
            {
                outlineCoord.push_back(outlinePt[indi].real()); // x-coordinate of an outline point
                outlineCoord.push_back(outlinePt[indi].imag()); // y-coordinate of an outline point
            }
            node outlineNode = node(outlineCoord, 1, {}, 0);
            outlineNode.reorder();
            vector<double> outlineCentroid = outlineNode.findCentroid(0., 0.);
            size_t indUnusedOutline = 0;

            // Write comment to file
            polyFile << "# File generated by limboint.h" << endl;
            polyFile << "# Layer " << layer << " has area of " << outlineNode.findArea() << " m^2" << endl;

            // Write vertices to file
            polyFile << "# vertices (vert_num, x, y, bound)" << endl;
            polyFile << numOutlinePt + numLayerPt << " 2 0 1" << endl; // Number of vertices (do not zero pad), 2 dimensions, 0 attributes, with boundary markers
            for (size_t indi = 0; indi < numOutlinePt; indi++)
            {
                polyFile << std::right << std::setw(6) << indi + 1 << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << outlinePt[indi].real() << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << outlinePt[indi].imag() << " 1" << endl; // Vertex # (do not zero pad), x-coordinate, y-coordinate, boundary marker
            }
            for (size_t indi = 0; indi < numLayerPt; indi++)
            {
                polyFile << std::right << std::setw(6) << numOutlinePt + indi + 1 << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << layerPt[indi].real() << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << layerPt[indi].imag() << " 0" << endl; // Vertex #, x-coordinate, y-coordinate, boundary marker
                if ((fabs((layerPt[indi].real() - outlineCoord[2 * indUnusedOutline]) / outlineCoord[2 * indUnusedOutline]) <= REL_INTERIOR_PRECISION) && (fabs((layerPt[indi].imag() - outlineCoord[2 * indUnusedOutline + 1]) / outlineCoord[2 * indUnusedOutline + 1]) <= REL_INTERIOR_PRECISION))
                {
                    indUnusedOutline++; // Future layer points must be checked against the next outline point (not robust if next outline point was a previous layer point)
                }
            }

            // Write segments to file
            polyFile << "# segments (seg_num, vert1, vert2, bound)" << endl;
            polyFile << numOutlinePt + numLayerSeg << " 1" << endl; // Number of segments, with boundary markers
            for (size_t indi = 0; indi < numOutlinePt; indi++)
            {
                polyFile << std::right << std::setw(6) << indi + 1 << " " << std::right << std::setw(6) << (indi % numOutlinePt) + 1 << " " << std::right << std::setw(6) << ((indi + 1) % numOutlinePt) + 1 << " 1" << endl; // Segment #, first vertex #, second vertex #, boundary marker
            }
            for (size_t indi = 0; indi < numLayerSeg; indi++)
            {
                polyFile << std::right << std::setw(6) << (numOutlinePt + indi + 1) << " " << std::right << std::setw(6) << (numOutlinePt + layerSeg[indi].first + 1) << " " << std::right << std::setw(6) << (numOutlinePt + layerSeg[indi].second + 1) << " 0" << endl; // Segment #, first vertex #, second vertex #, boundary marker
            }

            // Write holes to file
            polyFile << "# holes (hole_num, x, y)" << endl;
            polyFile << "0" << endl; // No holes in this design (future consideration)

            // Write regions to file
            polyFile << "# regional attributes (1 = conductor, 2 = dielectric)" << endl;
            polyFile << numLayerReg + 1 << endl; // Number of regions
            for (size_t indi = 0; indi < numLayerReg; indi++)
            {
                polyFile << std::right << std::setw(6) << indi + 1 << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << layerReg[indi].real() << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << layerReg[indi].imag() << " 1" << endl; // Region #, x-coordinate, y-coordinate, attribute number
            }
            vector<double> outlineInterior = { (1. - REL_INTERIOR_PRECISION) * outlineNode.getNodes()[2 * indUnusedOutline] + REL_INTERIOR_PRECISION * outlineCentroid[2 * indUnusedOutline], (1. - REL_INTERIOR_PRECISION) * outlineNode.getNodes()[2 * indUnusedOutline + 1] + REL_INTERIOR_PRECISION * outlineCentroid[2 * indUnusedOutline + 1] }; // Point just left and above the unused outine point (weighted average, default from bottom-right point)
            polyFile << std::right << std::setw(6) << numLayerReg + 1 << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << outlineInterior[0] << " " << std::setfill(' ') << std::showpos << std::setw(13) << std::setprecision(7) << outlineInterior[1] << " 2" << endl; // Region #, x-coordinate, y-coordinate, attribute number
            polyFile << endl; // Blank line for good measure

            // Close PSLG file
            polyFile.close();
            return true;
        }
        else
        {
            // File could not be opened
            return false;
        }
    }

    // Save to fdtdMesh conductor information
    void saveToMesh(std::string name, vector<double> center, strans transform, fdtdMesh *sys)
    {
        // Error checking on input point
        double xo, yo; // Coordinate offsets
        if (center.size() != 2)
        {
            cerr << "Coordinates of reference frame center must be a length-2 vector. Defaulting to (0, 0)." << endl;
            xo = 0.;
            yo = 0.;
        }
        else
        {
            xo = center[0];
            yo = center[1];
        }

        // Instantiate fdtdMesh variables
        int si, condj;    // Size of the vector sys->conductorIn, index of conductor in loop

        // Get information about this cell in ASCII database
        const GeoCell cell = this->cells[this->locateCell(name)];
        int numBound = cell.getNumBound();
        int numPath = cell.getNumPath();
        int numBox = cell.getNumBox();
        int numSRef = cell.getNumSRef();
        int numARef = cell.getNumARef();

        // Set conductor information from all elements in this geometric cell
        for (size_t indi = 0; indi < numBound; indi++) // Handle each boundary
        {
            vector<double> boundCoord = ((cell.boundaries)[indi]).getBounds();
            this->numCdtIn++;
            sys->conductorIn.push_back(fdtdOneCondct());
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
            for (size_t indj = 0; indj < boundCoord.size() - 2; indj += 2) // Iterate over each ordered pair, -2 because the last point is the starting point
            {
                vector<double> boundPt = transform.applyTranform({ boundCoord[indj], boundCoord[indj + 1] }); // Apply the linear transformation of this cell reference
                sys->conductorIn[si].x[condj] = boundPt[0] + xo;
                sys->conductorIn[si].y[condj] = boundPt[1] + yo;
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
        for (size_t indi = 0; indi < numPath; indi++) // Handle each path
        {
            vector<double> pathCoord = ((cell.paths)[indi]).getPaths();
            double width = ((cell.paths)[indi]).getWidth();
            int type = (cell.paths[indi]).getType();
            for (size_t indj = 0; indj < pathCoord.size() - 2; indj += 2) // Iterate over each ordered pair
            {
                // Account for the path termination type
                double overshoot1 = 1.0; // Assume path has round ends (type = 1) approximated as square or square ends with overshoot (type = 2)
                double overshoot2 = 1.0; // Assume path has round ends (type = 1) approximated as square or square ends with overshoot (type = 2)
                if ((type == 0) && (indj == 0))
                {
                    overshoot1 = 0.0; // Path actually has square ends at terminal vertices
                }
                if ((type == 0) && (indj == pathCoord.size() - 4))
                {
                    overshoot2 = 0.0; // Path actually has square ends at terminal vertices
                }

                // Initialize the conductorIn information
                sys->conductorIn.push_back(fdtdOneCondct());
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
                if (pathCoord[indj] == pathCoord[indj + 2]) // Path segment along y-axis
                {
                    if (pathCoord[indj + 1] > pathCoord[indj + 3]) // First point is above the second point
                    {
                        vector<double> pathPt = transform.applyTranform({ pathCoord[indj] - width / 2., pathCoord[indj + 1] + overshoot1 * width / 2. }); // Apply the linear transformation of this cell reference
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmin = sys->conductorIn[si].x[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj] + width / 2., pathCoord[indj + 1] + overshoot1 * width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmax = sys->conductorIn[si].x[condj];
                        sys->conductorIn[si].ymax = sys->conductorIn[si].y[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] + width / 2., pathCoord[indj + 3] - overshoot2 * width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] - width / 2., pathCoord[indj + 3] - overshoot2 * width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].ymin = sys->conductorIn[si].y[condj];
                        condj++;
                        this->numCdtIn++;
                    }
                    else // Second point is above the first point
                    {
                        vector<double> pathPt = transform.applyTranform({ pathCoord[indj] - width / 2., pathCoord[indj + 1] - overshoot1 * width / 2. }); // Apply the linear transformation of this cell reference
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmin = sys->conductorIn[si].x[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj] + width / 2., pathCoord[indj + 1] - overshoot1 * width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmax = sys->conductorIn[si].x[condj];
                        sys->conductorIn[si].ymin = sys->conductorIn[si].y[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] + width / 2., pathCoord[indj + 3] + overshoot2 * width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] - width / 2., pathCoord[indj + 3] + overshoot2 * width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].ymax = sys->conductorIn[si].y[condj];
                        condj++;
                        this->numCdtIn++;
                    }
                }
                else // Path segment along x-axis
                {
                    if (pathCoord[indj] > pathCoord[indj + 2]) // First point is on the right of the second point
                    {
                        vector<double> pathPt = transform.applyTranform({ pathCoord[indj] + overshoot1 * width / 2., pathCoord[indj + 1] + width / 2. }); // Apply the linear transformation of this cell reference
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmax = sys->conductorIn[si].x[condj];
                        sys->conductorIn[si].ymax = sys->conductorIn[si].y[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj] + overshoot1 * width / 2., pathCoord[indj + 1] - width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].ymin = sys->conductorIn[si].y[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] - overshoot2 * width / 2., pathCoord[indj + 3] - width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] - overshoot2 * width / 2., pathCoord[indj + 3] + width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmin = sys->conductorIn[si].x[condj];
                        condj++;
                        this->numCdtIn++;
                    }
                    else // Second point is on the right of the first point
                    {
                        vector<double> pathPt = transform.applyTranform({ pathCoord[indj] - overshoot1 * width / 2., pathCoord[indj + 1] + width / 2. }); // Apply the linear transformation of this cell reference
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmin = sys->conductorIn[si].x[condj];
                        sys->conductorIn[si].ymax = sys->conductorIn[si].y[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj] - overshoot1 * width / 2., pathCoord[indj + 1] - width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-left point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].ymin = sys->conductorIn[si].y[condj];
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] + overshoot2 * width / 2., pathCoord[indj + 3] - width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The lower-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        condj++;
                        pathPt = transform.applyTranform({ pathCoord[indj + 2] + overshoot2 * width / 2., pathCoord[indj + 3] + width / 2. });
                        sys->conductorIn[si].x[condj] = pathPt[0] + xo; // The upper-right point of path segment
                        sys->conductorIn[si].y[condj] = pathPt[1] + yo;
                        sys->conductorIn[si].xmax = sys->conductorIn[si].x[condj];
                        condj++;
                        this->numCdtIn++;
                    }
                }
            }
        }
        for (size_t indi = 0; indi < numBox; indi++) // Handle each box outline
        {
            vector<double> boxCoord = ((cell.boxes)[indi]).getBoxes();
            this->numCdtIn++;
            sys->conductorIn.push_back(fdtdOneCondct());
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
            for (size_t indj = 0; indj < boxCoord.size() - 2; indj += 2) // Iterate over each ordered pair, -2 because the last point is the starting point
            {
                vector<double> boxPt = transform.applyTranform({ boxCoord[indj], boxCoord[indj + 1] }); // Apply the linear transformation of this cell reference
                sys->conductorIn[si].x[condj] = boxPt[0] + xo;
                sys->conductorIn[si].y[condj] = boxPt[1] + yo;
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
        for (size_t indi = 0; indi < numSRef; indi++) // Handle each structure reference recursively
        {
            vector<double> refPt = transform.applyTranform(((cell.sreferences)[indi]).getSRefs()); // Apply the linear transformation of this cell reference
            saveToMesh((cell.sreferences)[indi].getSRefName(), { refPt[0] + xo, refPt[1] + yo }, (cell.sreferences)[indi].getTransform().composeTransform(transform), sys);
        }
        for (size_t indi = 0; indi < cell.getNumARef(); indi++) // Handle each array reference instance
        {
            vector<vector<double>> instanceCoord = (cell.areferences[indi]).findInstances({ 0., 0. });
            for (size_t indj = 0; indj < instanceCoord.size(); indj++) // Handle each instance in array reference recursively
            {
                vector<double> centPt = transform.applyTranform(instanceCoord[indj]); // Apply the linear transformation of this cell reference
                saveToMesh((cell.areferences)[indi].getARefName(), { centPt[0] + xo, centPt[1] + yo }, (cell.areferences)[indi].getTransform().composeTransform(transform), sys);
            }
        }
    }

    // Print the ASCII database with the design geometry
    void print(vector<size_t> indCellPrint)
    {
        // Analyze design and print to terminal
        int numCell = getNumCell();

        cout << "------" << endl;
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
            cout << "   Counts: " << (this->cells)[indi].getNumBound() << " boundaries, " << (this->cells)[indi].getNumPath() << " paths, " << (this->cells)[indi].getNumNode() << " nodes, " << (this->cells)[indi].getNumBox() << " boxes," << endl << "     " << (this->cells)[indi].getNumText() << " text boxes, " << (this->cells)[indi].getNumSRef() << " structure references, and " << (this->cells)[indi].getNumARef() << " array references" << endl;
        }
        for (size_t indi = 0; indi < indCellPrint.size(); indi++)
        {
            std::string cellName = ((this->cells)[indCellPrint[indi]]).getCellName();
            cout << cellName << endl;
            (this->cells)[indCellPrint[indi]].printAlt();
        }
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
                strans linTrans = (((this->cells)[indCell]).textboxes)[indText].getTransform();

                // Write the textbox to file
                gw.gds_write_text();
                gw.gds_write_layer((((this->cells)[indCell]).textboxes)[indText].getLayer());
                gw.gds_write_texttype((((this->cells)[indCell]).textboxes)[indText].getType());
                gw.gds_write_presentation((((this->cells)[indCell]).textboxes)[indText].getFontID(), (((this->cells)[indCell]).textboxes)[indText].getJusts()[0], (((this->cells)[indCell]).textboxes)[indText].getJusts()[1]);
                gw.gds_write_width(width);
                gw.gds_write_strans(linTrans.getMirrored(), linTrans.getAbsMagnify(), linTrans.getAbsRotate());
                gw.gds_write_mag(linTrans.getMagnify());
                gw.gds_write_angle(180. * M_1_PI * linTrans.getRotation());
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
                strans linTrans = (((this->cells)[indCell]).sreferences)[indSRef].getTransform();

                // Write the structure reference to file
                gw.gds_write_sref();
                gw.gds_write_sname((((this->cells)[indCell]).sreferences)[indSRef].getSRefName().c_str());
                gw.gds_write_strans(linTrans.getMirrored(), linTrans.getAbsMagnify(), linTrans.getAbsRotate());
                gw.gds_write_mag(linTrans.getMagnify());
                gw.gds_write_angle(180. * M_1_PI * linTrans.getRotation());
                gw.gds_write_xy(xcoord, ycoord, 1);
                gw.gds_write_endel();
            }

            // Create each array reference within the cell
            size_t numARef = ((this->cells)[indCell]).getNumARef();
            for (size_t indARef = 0; indARef < numARef; indARef++)
            {
                vector<double> arefs = (((this->cells)[indCell]).areferences)[indARef].getARefs();
                int xcoord[3] = { (int)(arefs[0] / this->getdbUnits()), (int)(arefs[2] / this->getdbUnits()), (int)(arefs[4] / this->getdbUnits()) };
                int ycoord[3] = { (int)(arefs[1] / this->getdbUnits()), (int)(arefs[3] / this->getdbUnits()), (int)(arefs[5] / this->getdbUnits()) };
                vector<int> numColRow = (((this->cells)[indCell]).areferences)[indARef].getNumColRow();
                strans linTrans = (((this->cells)[indCell]).areferences)[indARef].getTransform();

                // Write the array reference to file
                gw.gds_write_aref();
                gw.gds_write_sname((((this->cells)[indCell]).areferences)[indARef].getARefName().c_str());
                gw.gds_write_colrow(numColRow[0], numColRow[1]);
                gw.gds_write_strans(linTrans.getMirrored(), linTrans.getAbsMagnify(), linTrans.getAbsRotate());
                gw.gds_write_mag(linTrans.getMagnify());
                gw.gds_write_angle(180. * M_1_PI * linTrans.getRotation());
                gw.gds_write_xy(xcoord, ycoord, 3);
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
        /*if (true) //(this->getElement() == 'a')
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
        else if (ascii_record_type == "AREF")
        {
            (this->element) = 'a';
            ((this->cells)[this->numCell]).areferences.emplace_back(aref());
        }
        else if (ascii_record_type == "LAYER")
        {
            if (this->getElement() == 'b')
            {
                ((this->cells)[this->numCell]).boundaries.back().setLayer(data[0]); // Use setter to update layer of last boundary
            }
            else if (this->getElement() == 'p')
            {
                ((this->cells)[this->numCell]).paths.back().setLayer(data[0]);
            }
            else if (this->getElement() == 'n')
            {
                ((this->cells)[this->numCell]).nodes.back().setLayer(data[0]);
            }
            else if (this->getElement() == 'x')
            {
                ((this->cells)[this->numCell]).boxes.back().setLayer(data[0]);
            }
            else if (this->getElement() == 't')
            {
                ((this->cells)[this->numCell]).textboxes.back().setLayer(data[0]);
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
                ((this->cells)[this->numCell]).boundaries.back().setBounds(coord); // Use setter to update boundary points in cell coordinates of last boundary
            }
            else if (this->getElement() == 'p')
            {
                ((this->cells)[this->numCell]).paths.back().setPaths(coord);
            }
            else if (this->getElement() == 'n')
            {
                ((this->cells)[this->numCell]).nodes.back().setNodes(coord);
            }
            else if (this->getElement() == 'x')
            {
                ((this->cells)[this->numCell]).boxes.back().setBoxes(coord);
            }
            else if (this->getElement() == 't')
            {
                ((this->cells)[this->numCell]).textboxes.back().setTexts(coord);
            }
            else if (this->getElement() == 's')
            {
                ((this->cells)[this->numCell]).sreferences.back().setSRefs(coord);
            }
            else if (this->getElement() == 'a')
            {
                ((this->cells)[this->numCell]).areferences.back().setARefs(coord);
            }
        }
        else if (ascii_record_type == "DATATYPE") // Unimplemented in the GDSII standard and custom boundary and path classes
        {
            //cout << "Datatype for " << (this->element) << ": " << data[0] << endl;
        }
        else if (ascii_record_type == "PATHTYPE") // Necessary record in the GDSII standard to determine how path segments terminate
        {
            if (this->getElement() == 'p')
            {
                ((this->cells)[this->numCell]).paths.back().setType(data[0]);
            } // Unimplemented in custom textbox class
        }
        else if (ascii_record_type == "NODETYPE") // Should be identically zero in the GDSII standard
        {
            if (this->getElement() == 'n')
            {
                ((this->cells)[this->numCell]).nodes.back().setType(data[0]);
            }
        }
        else if (ascii_record_type == "BOXTYPE") // Should be identically zero in the GDSII standard
        {
            if (this->getElement() == 'x')
            {
                ((this->cells)[this->numCell]).boxes.back().setType(data[0]);
            }
        }
        else if (ascii_record_type == "TEXTTYPE") // Should be identically zero in the GDSII standard
        {
            if (this->getElement() == 't')
            {
                ((this->cells)[this->numCell]).textboxes.back().setType(data[0]);
            }
        }
        else if (ascii_record_type == "WIDTH")
        {
            // Scaling for widths
            double unitFactor = this->getdbUnits();

            // Store widths (negative values mean independent of structure scaling)
            if (this->getElement() == 'p')
            {
                ((this->cells)[this->numCell]).paths.back().setWidth(unitFactor * data[0]);
            }
            else if (this->getElement() == 't')
            {
                ((this->cells)[this->numCell]).textboxes.back().setWidth(unitFactor * data[0]);
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
                // Find font identifier based on bit pattern
                int fontID = ((int)data[0] & fontIDMask) >> 4;

                // Find vertical and horizontal justification based on bit pattern
                int vertJust = ((int)data[0] & vertJustMask) >> 2;
                int horizJust = ((int)data[0] & horizJustMask);
                if (vertJust == 3) // Only illegal value possible with bit mask
                {
                    cerr << "Illegal text box vertical justification, resetting to middle" << endl;
                    vertJust = 1;
                }
                if (horizJust == 3)
                {
                    cerr << "Illegal text box horizontal justification, resetting to center" << endl;
                    vertJust = 1;
                }

                // Update presentation in last text box
                ((this->cells)[this->numCell]).textboxes.back().setFontID(fontID);
                ((this->cells)[this->numCell]).textboxes.back().setJusts({ vertJust, horizJust });
            }
        }
        else if (ascii_record_type == "STRANS")
        {
            // Declare binary masks
            int mirrorMask = 0b1000000000000000;
            int absMagMask = 0b0000000000000100;
            int absRotMask = 0b0000000000000010;

            // Find if mirroring, absolute magnification, or absolute rotation used
            int mirrorInX = ((int)data[0] & mirrorMask) >> 15;
            int absMagnify = ((int)data[0] & absMagMask) >> 2;
            int absRotate = ((int)data[0] & absRotMask) >> 1;

            // Store mirroring, absolute magnification, and absolute rotation information
            if (this->getElement() == 't')
            {
                ((this->cells)[this->numCell]).textboxes.back().setTransform(strans((bool)mirrorInX, (bool)absMagnify, (bool)absRotate, 1.0, 0.0));
            }
            else if (this->getElement() == 's')
            {
                ((this->cells)[this->numCell]).sreferences.back().setTransform(strans((bool)mirrorInX, (bool)absMagnify, (bool)absRotate, 1.0, 0.0));
            }
            else if (this->getElement() == 'a')
            {
                ((this->cells)[this->numCell]).areferences.back().setTransform(strans((bool)mirrorInX, (bool)absMagnify, (bool)absRotate, 1.0, 0.0));
            }
        }
        else if (ascii_record_type == "MAG")
        {
            if (this->getElement() == 't')
            {
                // Get copy of existing linear transformation of last textbox
                strans modTrans = ((this->cells)[this->numCell]).textboxes.back().getTransform();

                // Modify linear transform and store
                modTrans.setMagnify(data[0]);
                ((this->cells)[this->numCell]).textboxes.back().setTransform(modTrans);
            }
            else if (this->getElement() == 's')
            {
                // Get copy of existing linear transformation of last structure reference
                strans modTrans = ((this->cells)[this->numCell]).sreferences.back().getTransform();

                // Modify linear transform and store
                modTrans.setMagnify(data[0]);
                ((this->cells)[this->numCell]).sreferences.back().setTransform(modTrans);
            }
            else if (this->getElement() == 'a')
            {
                // Get copy of existing linear transformation of last array reference
                strans modTrans = ((this->cells)[this->numCell]).areferences.back().getTransform();

                // Modify linear transform and store
                modTrans.setMagnify(data[0]);
                ((this->cells)[this->numCell]).areferences.back().setTransform(modTrans);
            }
        }
        else if (ascii_record_type == "ANGLE")
        {
            if (this->getElement() == 't')
            {
                // Get copy of existing linear transformation of last textbox
                strans modTrans = ((this->cells)[this->numCell]).textboxes.back().getTransform();

                // Modify linear transform and store
                modTrans.setRotation(M_PI / 180. * data[0]);
                ((this->cells)[this->numCell]).textboxes.back().setTransform(modTrans);
            }
            else if (this->getElement() == 's')
            {
                // Get copy of existing linear transformation of last structure reference
                strans modTrans = ((this->cells)[this->numCell]).sreferences.back().getTransform();

                // Modify linear transform and store
                modTrans.setRotation(M_PI / 180. * data[0]);
                ((this->cells)[this->numCell]).sreferences.back().setTransform(modTrans);
            }
            else if (this->getElement() == 'a')
            {
                // Get copy of existing linear transformation of last array reference
                strans modTrans = ((this->cells)[this->numCell]).areferences.back().getTransform();

                // Modify linear transform and store
                modTrans.setRotation(M_PI / 180. * data[0]);
                ((this->cells)[this->numCell]).areferences.back().setTransform(modTrans);
            }
        }
        else if (ascii_record_type == "STRING")
        {
            // Print and store text box string
            char label[512];
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
                ((this->cells)[this->numCell]).textboxes.back().setTextStr(label);
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
            //cout << "Reference cell name: " << label << endl;
            if (this->getElement() == 's')
            {
                ((this->cells)[this->numCell]).sreferences.back().setSRefName(label);

            }
            else if (this->getElement() == 'a')
            {
                ((this->cells)[this->numCell]).areferences.back().setARefName(label);
            }
        }
        else if (ascii_record_type == "COLROW")
        {
            if (this->getElement() == 'a')
            {
                ((this->cells)[this->numCell]).areferences.back().setNumColRow({ (int)data[0], (int)data[1] });
            }
        }
        else if (ascii_record_type == "PROPATTR")
        {
            this->numProp = data[0]; // Store position for attribute value to be saved
        }
        else if (ascii_record_type == "PROPVALUE")
        {
            // Fix and print property
            char label[128];
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

            // Store property
            if (this->getElement() == 'b')
            {
                vector<std::string> modProps = ((this->cells)[this->numCell]).boundaries.back().getProps(); // Get copy of existing properties of last boundary
                //cout << " Retreived copy of properties of size " << modProps.size() << " to try and place \"" << label << "\" at position " << this->getNumProp() << endl;
                if (this->numProp > modProps.size())
                {
                    for (size_t indProp = modProps.size(); indProp < this->numProp; indProp++)
                    {
                        // Keep pushing back empty strings until the correct property index is reached
                        modProps.push_back("");
                    }
                }
                modProps[this->numProp - 1] = label; // Add new property at specified location to copy
                ((this->cells)[this->numCell]).boundaries.back().setProps(modProps); // Use setter to update properties of last boundary
            }
            else if (this->getElement() == 'p')
            {
                vector<std::string> modProps = ((this->cells)[this->numCell]).paths.back().getProps();
                if (this->numProp > modProps.size())
                {
                    for (size_t indProp = modProps.size(); indProp < this->numProp; indProp++)
                    {
                        // Keep pushing back empty strings until the correct property index is reached
                        modProps.push_back("");
                    }
                }
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).paths.back().setProps(modProps);
            }
            else if (this->getElement() == 'n')
            {
                vector<std::string> modProps = ((this->cells)[this->numCell]).nodes.back().getProps();
                if (this->numProp > modProps.size())
                {
                    for (size_t indProp = modProps.size(); indProp < this->numProp; indProp++)
                    {
                        // Keep pushing back empty strings until the correct property index is reached
                        modProps.push_back("");
                    }
                }
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).nodes.back().setProps(modProps);
            }
            else if (this->getElement() == 'x')
            {
                vector<std::string> modProps = ((this->cells)[this->numCell]).boxes.back().getProps();
                if (this->numProp > modProps.size())
                {
                    for (size_t indProp = modProps.size(); indProp < this->numProp; indProp++)
                    {
                        // Keep pushing back empty strings until the correct property index is reached
                        modProps.push_back("");
                    }
                }
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).boxes.back().setProps(modProps);
            }
            else if (this->getElement() == 't')
            {
                vector<std::string> modProps = ((this->cells)[this->numCell]).textboxes.back().getProps();
                if (this->numProp > modProps.size())
                {
                    for (size_t indProp = modProps.size(); indProp < this->numProp; indProp++)
                    {
                        // Keep pushing back empty strings until the correct property index is reached
                        modProps.push_back("");
                    }
                }
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).textboxes.back().setProps(modProps);
            }
            else if (this->getElement() == 's')
            {
                vector<std::string> modProps = ((this->cells)[this->numCell]).sreferences.back().getProps();
                if (this->numProp > modProps.size())
                {
                    for (size_t indProp = modProps.size(); indProp < this->numProp; indProp++)
                    {
                        // Keep pushing back empty strings until the correct property index is reached
                        modProps.push_back("");
                    }
                }
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).sreferences.back().setProps(modProps);
            }
            else if (this->getElement() == 'a')
            {
                vector<std::string> modProps = ((this->cells)[this->numCell]).areferences.back().getProps();
                if (this->numProp > modProps.size())
                {
                    for (size_t indProp = modProps.size(); indProp < this->numProp; indProp++)
                    {
                        // Keep pushing back empty strings until the correct property index is reached
                        modProps.push_back("");
                    }
                }
                modProps[this->numProp - 1] = label;
                ((this->cells)[this->numCell]).areferences.back().setProps(modProps);
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
