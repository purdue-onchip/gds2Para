#pragma once
#ifndef AUTOPORTFROMDEFLEF_H_
#define AUTOPORTFROMDEFLEF_H_

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>

#include <limbo/parsers/def/adapt/DefDriver.h>  // Limbo DEF parser
#include <limbo/parsers/lef/adapt/LefDriver.h>  // Limbo LEF parser
//#include "solnoutclass.hpp"                     // class Port
#include "limboint.hpp"                         // class strans

using namespace std;

enum OrientType {
    N = 0,
    S = 1,
    W = 2,
    E = 3,
    FN = 4,
    FS = 5,
    FW = 6,
    FE = 7
};

struct ComponentInfo {
    double xInUm = 0;   // x, y of origin, converted to default DEF unit in um.
    double yInUm = 0;
    string orient = "N";
    string cellName = "";
};

struct PinInfo {
    string direct = "INPUT";    // direction, {INPUT, OUTPUT, INOUT, FEEDTHRU}
    vector<string> vLayer;      // layers  

    PinInfo() {}
    PinInfo(string dir, vector<string> vLay)
        : direct{ dir }, vLayer{ vLay } {}
};

struct DefPinInfo : public PinInfo {
    double xInUm = 0;   // x, y placement of this Pin, converted to default DEF unit in um.
    double yInUm = 0;

    DefPinInfo() {}
    DefPinInfo(string dir, vector<string> vLay, double x, double y)
        : PinInfo{ dir, vLay }, xInUm{ x }, yInUm{ y } {}
};

struct NetInfo {        // same as DefParser::Net, redefine here for easier reference.
    string netName = "";
    int netWeight = 1;  // net weight, used in automatic layout tools, not used in gds2Para.
    vector<std::pair<string, string> > vNodenamePin;    // array of (node, pin) pair. vNetPin[i].first = componentName.
};

/// @brief Custom class that inheritates @ref DefParser::DefDataBase 
/// with all the required callbacks defined. 
class DefDataBase : public DefParser::DefDataBase
{
public:
    vector<NetInfo> allNets;    // need to track the order of nets in DEF file, so use vector.
    unordered_map<string, ComponentInfo> allComponents; // map[componentName] = struct{cellName, x, y, orient}
    unordered_map<string, DefPinInfo> allDefPins;       // map[pinName] = struct{x, y, direction, layer}
    string defVersion = "";
    string defDesign = "";
    int defUnit = 1;        // int coordinate in DEF divided by this->defUnit will be true physical coordinate in unit um. 
    double dieAreaInUm[4];  // {xmin, ymin, xmax, ymax} in um of this design 

    DefDataBase() {}

    //////////////////// Custom member functions defined here ///////////////////

    void print_allNets();
    void print_allDefPins();
    void print_allComponents();

    // check if all nodeNames defined in nets are valid component or valid external pin
    bool areAllNetNodesValidComponentOrValidPin();

    //////////////////// required callbacks from abstract DefParser::DefDataBase ///////////////////
    /// @param token divider characters 
    virtual void set_def_dividerchar(string const& token);
    /// @param token BUS bit characters 
    virtual void set_def_busbitchars(string const& token);
    /// @param token DEF version 
    virtual void set_def_version(string const& token);
    /// @param token design name 
    virtual void set_def_design(string const& token);
    /// @param token DEF unit 
    virtual void set_def_unit(int token);
    /// @param t1, t2, t3, t4 die area (xl, yl, xh, yh)
    virtual void set_def_diearea(int t1, int t2, int t3, int t4);
    /// @brief add row 
    virtual void add_def_row(DefParser::Row const&);
    /// @brief add component 
    /// @param c component 
    virtual void add_def_component(DefParser::Component const& c);
    /// @param token number of components 
    virtual void resize_def_component(int token);
    /// @brief add pin 
    /// @param p pin 
    virtual void add_def_pin(DefParser::Pin const& p);
    /// @brief set number of pins 
    /// @param token number of pins 
    virtual void resize_def_pin(int token);
    /// @brief add net 
    /// @param n net 
    virtual void add_def_net(DefParser::Net const& n);
    /// @brief set number of nets 
    /// @param token number of nets 
    virtual void resize_def_net(int token);
    /// @brief set number of blockages 
    /// @param n number of blockages 
    virtual void resize_def_blockage(int n);
    /// @brief add placement blockages 
    /// @param vBbox array of boxes with xl, yl, xh, yh
    virtual void add_def_placement_blockage(std::vector<std::vector<int> > const& vBbox);
    /// @brief end of design 
    virtual void end_def_design();
};

struct LefPinInfo : public PinInfo {
    vector<vector<double>> vRectsInUm;      // rectangles. vRects[i] = {x0,y0,x1,y1} are two opposite corners of the rect

    LefPinInfo() {}
    LefPinInfo(string dir, vector<string> vLay, vector<vector<double>> vRect)
        : PinInfo{ dir, vLay }, vRectsInUm{ vRect } {}
};

struct LefCellInfo {
    double originXInUm = 0; // origin of this cell to align with a DEF COMPONENT placement point, in unit um.
    double originYInUm = 0;
    double sizeXInUm = 0;   // width BY height of the placement bounding rectangle, in um.
    double sizeYInUm = 0;
    unordered_map<string, LefPinInfo> allPinsInCell;    // map[pinName] = struct LefPinInfo.
};

/// @brief test LefParser
class LefDataBase : public LefParser::LefDataBase
{
public:
    unordered_map<string, LefCellInfo> allCells;        // map[cellName] = struct LefCellInfo.
    string lefVersion = "";

protected:
    string tempCellName = "";       // temp data to store info of current cell
    unordered_map<string, LefPinInfo> tempPinsInCell;


public:
    //////////////////// Custom member functions defined here ///////////////////

    void appendCellMap(const unordered_map<string, LefCellInfo>& newCells);
    void print_allCells();

    /// base type 
    typedef LefParser::LefDataBase base_type;
    /// @brief constructor 
    LefDataBase() : base_type() {}
    //////////////////// required callbacks from abstract LefParser::LefDataBase ///////////////////
    /// @brief set LEF version 
    /// @param v string of LEF version 
    virtual void lef_version_cbk(string const& v);
    /// @brief set LEF version 
    /// @param v floating point number of LEF version 
    virtual void lef_version_cbk(double v);
    /// @brief set divider characters 
    /// @param v divider characters
    virtual void lef_dividerchar_cbk(string const& v);
    /// @brief set unit 
    /// @param v an object for unit 
    virtual void lef_units_cbk(LefParser::lefiUnits const& v);
    /// @brief set manufacturing entry 
    /// @param v manufacturing entry 
    virtual void lef_manufacturing_cbk(double v);
    /// @brief set bus bit characters 
    /// @param v but bit characters 
    virtual void lef_busbitchars_cbk(string const& v);
    /// @brief add layer 
    /// @param v an object for layer 
    virtual void lef_layer_cbk(LefParser::lefiLayer const& v);
    /// @brief add via 
    /// @param v an object for via 
    virtual void lef_via_cbk(LefParser::lefiVia const& v);
    /// @brief add via rule 
    /// @param v an object for via rule 
    virtual void lef_viarule_cbk(LefParser::lefiViaRule const& v);
    /// @brief spacing callback 
    /// @param v an object for spacing 
    virtual void lef_spacing_cbk(LefParser::lefiSpacing const& v);
    /// @brief site callback 
    /// @param v an object for site 
    virtual void lef_site_cbk(LefParser::lefiSite const& v);
    /// @brief macro begin callback, describe standard cell type 
    /// @param v name of macro 
    virtual void lef_macrobegin_cbk(std::string const& v);
    /// @brief macro callback, describe standard cell type 
    /// @param v an object for macro 
    virtual void lef_macro_cbk(LefParser::lefiMacro const& v);
    /// @brief property callback 
    /// @param v an object for property 
    virtual void lef_prop_cbk(LefParser::lefiProp const& v);
    /// @brief noise margin callback 
    /// @param v an object for noise margin 
    virtual void lef_maxstackvia_cbk(LefParser::lefiMaxStackVia const& v);
    /// @brief obstruction callback 
    /// @param v an object for obstruction 
    virtual void lef_obstruction_cbk(LefParser::lefiObstruction const& v);
    /// @brief pin callback, describe pins in a standard cell 
    /// @param v an object for pin 
    virtual void lef_pin_cbk(lefiPin const& v);
};

// check if the cellName and LefPinName used as net node are correctly defined in LEF files
bool areAllComponentsInNetsValidCell(
    const unordered_map<string, ComponentInfo>& allComponentsDEF,
    const vector<NetInfo>& allNetsDEF,
    const unordered_map<string, LefCellInfo>& allCellsLEF);


struct LayerMapInfo {
    int layerNameInNum = -1;    // used in GDSII
    double zminInUm = 0.0;
    double zmaxInUm = 0.0;
};
unordered_map<string, LayerMapInfo> readLayerMap(string inFileName);
void print_layerMap(const unordered_map<string, LayerMapInfo>& layerMap);

/* Directly reusing class "Port" in "solnoutclass.hpp" needs inporting 
the entire big "solnoutclass.hpp" and extra packages. So, define a 
simple struct NetPortCoor here, which can be easily converted to class Port when needed.*/
struct NetPortCoor : public PinInfo {
    string portName;    // netName + NodeName + pinName
    double xInUm = 0;   // x, y placement of this Pin, converted to default DEF unit in um.
    double yInUm = 0;
    double zInUm = 0;
    int gdsiiNum = -1;  // Layer number in GDSII file on which port exists
};

class AutoPorts {
public:
    unordered_map<string, LayerMapInfo> layerMap;
    unordered_map<string, vector<NetPortCoor>> netName_to_vPortCoor;


    AutoPorts();

    void readLayerMap_cbk(string inFileName);
    void print_layerMap_cbk();
    void getPortCoordinate(
        const unordered_map<string, ComponentInfo>& allComponentsDEF,
        const unordered_map<string, DefPinInfo>& allDefPinsDEF,
        const vector<NetInfo>& allNetsDEF,
        const unordered_map<string, LefCellInfo>& allCellsLEF);
    void print_netName_to_vPortCoor();
};
#endif
