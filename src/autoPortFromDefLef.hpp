#pragma once
#ifndef AUTOPORTFROMDEFLEF_H_
#define AUTOPORTFROMDEFLEF_H_

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

#include <limbo/parsers/def/adapt/DefDriver.h>  // Limbo DEF parser
#include <limbo/parsers/lef/adapt/LefDriver.h>  // Limbo LEF parser
//#include "solnoutclass.hpp"                     // class Port
//#include "limboint.hpp"                         // class strans

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

struct PinInfo {
    string direct = "INPUT";    // direction, {INPUT, OUTPUT, INOUT, FEEDTHRU}
    vector<string> vLayer;      // layers  

    PinInfo() {}
    PinInfo(string dir, vector<string> vLay)
        : direct{ dir }, vLayer{ vLay } {}
};

struct CompPinInfo : public PinInfo {
    // rectangles of each pin in the component, converted to global coordinate in um.
    vector<vector<double>> vRectsInUm_globCoor; // vRects[i] = {x0,y0,x1,y1} are two opposite corners of the rect
};

struct ComponentInfo {
    double xInUm = 0;   // x, y of origin, converted to default DEF unit in um.
    double yInUm = 0;
    string orient = "N";
    string cellName = "";
    unordered_map<string, CompPinInfo> allPinsInComp;   // map[pinName] = struct CompPinInfo.
};

struct DefPinInfo : public PinInfo {
    double xInUm = 0;   // x, y placement of this Pin, converted to default DEF unit in um.
    double yInUm = 0;

    DefPinInfo() {}
    DefPinInfo(string dir, vector<string> vLay, double x, double y)
        : PinInfo{ dir, vLay }, xInUm{ x }, yInUm{ y } {}
};

struct ViaInfo {
    string viaName = "";
    double xInUm;
    double yInUm;
};

struct NetInfo {        // same as DefParser::Net, redefine here for easier reference.
    string netName = "";
    int netWeight = 1;  // net weight, used in automatic layout tools, not used in gds2Para.
    vector<std::pair<string, string> > vNodenamePin;    // array of (node, pin) pair. vNetPin[i].first = componentName.
};

class DefDataBase : public DefParser::DefDataBase
{
public:
    vector<NetInfo> allNets;    // need to track the order of nets in DEF file, so use vector.
    unordered_map<string, ComponentInfo> allComponents; // map[componentName] = struct{cellName, x, y, orient}
    unordered_map<string, DefPinInfo> allDefPins;       // map[pinName] = struct{x, y, direction, layer}
    unordered_map<string, vector<ViaInfo>> netName_to_vVias;
    string defVersion = "";
    string defDesign = "";
    int defUnit = 1;        // int coordinate in DEF divided by this->defUnit will be true physical coordinate in unit um. 
    double dieAreaInUm[4];  // {xmin, ymin, xmax, ymax} in um of this design 

    DefDataBase() {}

    void print_allNets();
    void print_allDefPins();
    void print_allComponents();

    // check if all nodeNames defined in nets are valid component or valid external pin
    bool areAllNetNodesValidComponentOrValidPin();

    // required callbacks from abstract DefParser::DefDataBase
    virtual void set_def_dividerchar(string const& token);
    virtual void set_def_busbitchars(string const& token);
    virtual void set_def_version(string const& token);
    virtual void set_def_design(string const& token);
    virtual void set_def_unit(int token);
    virtual void set_def_diearea(int t1, int t2, int t3, int t4);
    virtual void add_def_row(DefParser::Row const&);
    virtual void add_def_component(DefParser::Component const& c);
    virtual void resize_def_component(int token);
    virtual void add_def_pin(DefParser::Pin const& p);
    virtual void resize_def_pin(int token);
    virtual void add_def_net(DefParser::Net const& n);
    virtual void resize_def_net(int token);
    virtual void resize_def_blockage(int n);
    virtual void add_def_placement_blockage(std::vector<std::vector<int> > const& vBbox);
    virtual void end_def_design();
};

namespace CustomDefParser {
    class CustomDefDriver : public DefParser::Driver {
    public:
        CustomDefDriver(DefDataBase& dbDef);
        bool custom_parse_file(const string& filename);
    };

    bool customDefRead(DefDataBase& dbDef, const string& defFile);
}   // end of namespace CustomDefParser

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

class LefDataBase : public LefParser::LefDataBase
{
public:
    unordered_map<string, LefCellInfo> allCells;        // map[cellName] = struct LefCellInfo.
    string lefVersion = "";

protected:
    string tempCellName = "";       // temp data to store info of current cell
    unordered_map<string, LefPinInfo> tempPinsInCell;

public:
    void appendCellMap(const unordered_map<string, LefCellInfo>& newCells);
    void print_allCells();

    typedef LefParser::LefDataBase base_type;
    LefDataBase() : base_type() {}

    // required callbacks from abstract LefParser::LefDataBase
    virtual void lef_version_cbk(string const& v);
    virtual void lef_version_cbk(double v);
    virtual void lef_dividerchar_cbk(string const& v);
    virtual void lef_units_cbk(LefParser::lefiUnits const& v);
    virtual void lef_manufacturing_cbk(double v);
    virtual void lef_busbitchars_cbk(string const& v);
    virtual void lef_layer_cbk(LefParser::lefiLayer const& v);
    virtual void lef_via_cbk(LefParser::lefiVia const& v);
    virtual void lef_viarule_cbk(LefParser::lefiViaRule const& v);
    virtual void lef_spacing_cbk(LefParser::lefiSpacing const& v);
    virtual void lef_site_cbk(LefParser::lefiSite const& v);
    virtual void lef_macrobegin_cbk(std::string const& v);
    virtual void lef_macro_cbk(LefParser::lefiMacro const& v);
    virtual void lef_prop_cbk(LefParser::lefiProp const& v);
    virtual void lef_maxstackvia_cbk(LefParser::lefiMaxStackVia const& v);
    virtual void lef_obstruction_cbk(LefParser::lefiObstruction const& v);
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
        const unordered_map<string, vector<ViaInfo>> netName_to_vVias,
        const unordered_map<string, ComponentInfo>& allComponentsDEF,
        const unordered_map<string, DefPinInfo>& allDefPinsDEF,
        const vector<NetInfo>& allNetsDEF);
    void print_netName_to_vPortCoor(const vector<NetInfo>& allNetsDEF);
};

vector<double> localLefCoorToGlobalDefCoor(double localLefCoor[2], 
    double sizeBB[2], double origin[2], double placementDef[2], string orient);
void localCellPinRect_to_globalCompPinRect(
    const unordered_map<string, LefCellInfo>& allCellsLEF,
    unordered_map<string, ComponentInfo>& allComponentsDEF);
string allDigitsInString(const string& str);

#endif
