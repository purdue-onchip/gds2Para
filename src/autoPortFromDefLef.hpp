#pragma once
#ifndef AUTOPORTFROMDEFLEF_H_
#define AUTOPORTFROMDEFLEF_H_

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>

#include <limbo/parsers/def/adapt/DefDriver.h>

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
    double xInUm = 0;   // x, y of origin, converted to unit um, as default DEF unit.
    double yInUm = 0;
    string orient = "";
    string cellName = "";
};

struct PinInfo {    // external pin defined in DEF, copied from DefParser::Pin
    double xInUm = 0;   // x, y of this Pin, converted to unit um, as default DEF unit.
    double yInUm = 0;
    //string pin_name; ///< pin name 
    //string net_name; ///< net name 
    string direct; //< direction, {INPUT, OUTPUT, INOUT, FEEDTHRU}
    //string status; ///< placement status , {COVER, FIXED, PLACED}
    //int32_t origin[2]; ///< offset to node origin 
    //string orient; ///< orientation 
    vector<string> vLayer; ///< layers  
    //vector<vector<int32_t> > vBbox; ///< bounding box on each layer, pin geometry
    //string use; ///< "use" token in DEF file, {SIGNAL, POWER, GROUND, CLOCK, etc.}
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
    unordered_map<string, PinInfo> allPins;             // map[pinName] = struct{x, y, direction, layer}
    string defVersion = "";
    string defDesign = "";
    int defUnit = 1;        // int coordinate in DEF divided by this->defUnit will be true physical coordinate in unit um. 
    double dieAreaInUm[4];  // {xmin, ymin, xmax, ymax} in um of this design 

    DefDataBase()
    {
        //cout << "DefDataBase::" << __func__ << endl;
    }

    //////////////////// Custom member functions defined here ///////////////////

    void print_allNets();
    void print_allPins();
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


#endif
