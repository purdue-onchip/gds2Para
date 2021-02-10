#include "autoPortFromDefLef.hpp"

namespace CustomDefParser {
    CustomDefDriver::CustomDefDriver(DefDataBase& dbDef)
        : Driver(dbDef) {
    }

    // define as global to access retrieved data inside customized callback function
    DefDataBase* defDB = NULL;

    int custom_netf(defrCallbackType_e c, defiNet* net, defiUserData ud) {
#define IGNORE_DEBUG_PRINT_defiNet_
#ifndef IGNORE_DEBUG_PRINT_defiNet_
        cout << endl << endl;
        net->print(stdout);
        cout << "numWires = " << net->numWires() << endl;
#endif
        
        string netName = net->name();
        vector<ViaInfo> vVias;

        // Via information from Regular Wiring Statement in DEF-NETS
        /*      refer to Limbo\limbo\thirdparty\lefdef\5.8\def\def\defiNet.hpp  */
        int x, y;       // current point
        ViaInfo curVia;
        for (int i = 0; i < net->numWires(); i++) {
            int newLayer = 0;
            const defiWire* w = net->wire(i);
#ifndef IGNORE_DEBUG_PRINT_defiNet_
            fprintf(stdout, "+ %s ", w->wireType());
#endif
            for (int j = 0; j < w->numPaths(); j++) {
                const defiPath* p = w->path(j);
                p->initTraverse();
                int path;
                while ((path = (int)(p->next())) != DEFIPATH_DONE) {
                    switch (path) {
                    case DEFIPATH_LAYER:
                        if (newLayer == 0) {
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                            fprintf(stdout, "%s ", p->getLayer());
#endif
                            newLayer = 1;
                        }
                        else {
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                            fprintf(stdout, "\n  NEW %s ", p->getLayer());
#endif
                        }
                        break;
                    case DEFIPATH_VIA:
                        curVia.viaName = p->getVia();
                        curVia.xInUm = x * 1.0 / defDB->defUnit;
                        curVia.yInUm = y * 1.0 / defDB->defUnit;
                        vVias.push_back(curVia);
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                        fprintf(stdout, "%s ", p->getVia());
#endif
                        break;
                    case DEFIPATH_VIAROTATION:
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                        fprintf(stdout, "%d\n", p->getViaRotation());
#endif
                        break;
                    case DEFIPATH_VIADATA:
                        int numX, numY, stepX, stepY;
                        p->getViaData(&numX, &numY, &stepX, &stepY);
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                        fprintf(stdout, "%d %d %d %d\n", numX, numY, stepX, stepY);
#endif
                        break;
                    case DEFIPATH_WIDTH:
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                        fprintf(stdout, "%d\n", p->getWidth());
#endif
                        break;
                    case DEFIPATH_POINT:
                        p->getPoint(&x, &y);
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                        fprintf(stdout, "( %d %d ) ", x, y);
#endif
                        break;
                    case DEFIPATH_TAPER:
#ifndef IGNORE_DEBUG_PRINT_defiNet_
                        fprintf(stdout, "TAPER\n");
#endif
                        break;
                    }
                }
            }
        }

        defDB->netName_to_vVias[netName] = vVias;
        
        return 0;
    }

    bool CustomDefDriver::custom_parse_file(const string& filename) {
        // Setup DEF reader API. Doc in Limbo/limbo/thirdparty/lefdef/5.8/def/doc/defapi.pdf
        defrInitSession(1);
        defrSetUserData((void*)3);  // Set user data
        defrSetAddPathToNet();      // Add path data to the appropriate net data
        
        // Set customized callback function to retrieve info from thirdparty DEF reader
        /* Available function cbk and function input type are defined in 
                Limbo\limbo\thirdparty\lefdef\5.8\def\def\defrReader.hpp
           The usage example in Limbo is in
                Limbo\limbo\parsers\def\adapt\DefDriver.cc
        */
        defrSetNetCbk(custom_netf);

        FILE* f = fopen(filename.c_str(), "r");
        if (!f) {
            std::cerr << "Could not open input file " << filename << "\n";
            return false;
        }

        // Set case sensitive to 0 to start with, in History & PropertyDefinition
        // reset it to 1.
        void* userData = NULL;
        int res = defrRead(f, filename.c_str(), userData, 1);

        (void)defrPrintUnusedCallbacks(stdout);
        (void)defrReleaseNResetMemory();
        (void)defrUnsetNonDefaultCbk();
        (void)defrUnsetNonDefaultStartCbk();
        (void)defrUnsetNonDefaultEndCbk();

        // Unset all the callbacks
        defrUnsetNetCbk();

        // Release allocated singleton data.
        defrClear();
        fclose(f);

        return true;
    }

    bool customDefRead(DefDataBase& dbDef, const string& defFile) {
        defDB = &dbDef;      // set dbDef as global defDB 

        CustomDefDriver driver(dbDef);
        bool res = driver.custom_parse_file(defFile);

        defDB = NULL;       // reset to NULL
        return res;
    }
}   // end of namespace CustomDefParser


void DefDataBase::print_allNets() {
    cout << "Total " << this->allNets.size() << " nets. netName: (nodeName, pin), viaName, viaCoor in um \n";
    for (const auto& net : this->allNets) {
        string netName = net.netName;
        cout << "Net \"" << netName << "\": \n";
        for (const auto& NodenamePin : net.vNodenamePin)
            cout << "  (" << NodenamePin.first << ", " << NodenamePin.second << ") ";
        cout << endl;
        if (this->netName_to_vVias.find(netName) == this->netName_to_vVias.end()) {
            cout << "  No via found in net " << netName << endl;
        }
        else {
            for (const ViaInfo& via : this->netName_to_vVias[netName]) {
                cout << "  " << via.viaName << "  ( " << via.xInUm << ", " << via.yInUm << " ) \n";
            }
        }
    }
}

void DefDataBase::print_allDefPins() {
    cout << "Total " << this->allDefPins.size() << " external pins. pinName (unordered): x(um), y(um), direction, layers\n";
    for (const auto& pin : this->allDefPins) {
        cout << pin.first << ": " << pin.second.xInUm << ", " << pin.second.yInUm << ", "
            << pin.second.direct << ", ";
        for (const auto& layer : pin.second.vLayer) cout << layer << "  ";
        cout << endl;
    }
}

void DefDataBase::print_allComponents() {
    cout << "Total " << this->allComponents.size() << " components. compName (unordered): cellName, x(um), y(um), orient \n";
    for (const auto& comp : this->allComponents) {
        cout << comp.first << ": " << comp.second.cellName << ", " << comp.second.xInUm << ", "
            << comp.second.yInUm << ", " << comp.second.orient << endl;
    }
}

// check if all nodeNames defined in nets are valid component or valid external pin
bool DefDataBase::areAllNetNodesValidComponentOrValidPin() {
    bool isValidCom = true;
    for (const auto& net : this->allNets) {
        for (const auto& NodenamePin : net.vNodenamePin) {
            const string& nodeName = NodenamePin.first;
            const string& pinName = NodenamePin.second;

            // Node type 1: defined as component
            if ((nodeName != "PIN") && (this->allComponents.find(nodeName) == this->allComponents.end())) {
                isValidCom = false;
                cout << "Net \"" << net.netName << "\" has invalid node name \"" << nodeName << "\"\n";
            }

            // Node type 2: defined as external pin
            if ((nodeName == "PIN") && (this->allDefPins.find(pinName) == this->allDefPins.end())) {
                isValidCom = false;
                cout << "Net \"" << net.netName << "\" has invalid node name at pin \"" << nodeName << "\"\n";
            }
        }
    }

    return isValidCom;
}

// check if the cellName and LefPinName used as net node are correctly defined in LEF files
bool areAllComponentsInNetsValidCell(
    const unordered_map<string, ComponentInfo>& allComponentsDEF,
    const vector<NetInfo>& allNetsDEF,
    const unordered_map<string, LefCellInfo>& allCellsLEF) 
{
    bool isValid = true;

    for (const auto& net : allNetsDEF) {
        for (const auto& NodenamePin : net.vNodenamePin) {
            const string& nodeName = NodenamePin.first;
            const string& pinName = NodenamePin.second;

            // Node type 1: defined as component
            if (nodeName != "PIN") {
                // check cellName
                string cellName = allComponentsDEF.at(nodeName).cellName;
                if (allCellsLEF.find(cellName) == allCellsLEF.end()) {
                    isValid = false;
                    cout << "Node \"" << nodeName << "\" of net \"" << net.netName << "\" has invalid cell \"" << cellName << "\"\n";
                }

                // check pinName
                unordered_map<string, LefPinInfo> allPinsInCell = allCellsLEF.at(cellName).allPinsInCell;
                if (allPinsInCell.find(pinName) == allPinsInCell.end()) {
                    isValid = false;
                    cout << "Node \"" << nodeName << "\" of net \"" << net.netName << "\" has invalid pin \"" << pinName << "\"\n";
                }
            }
        }
    }

    return isValid;
}

//////////////////// required callbacks from abstract DefParser::DefDataBase ///////////////////
/// @param token divider characters 
void DefDataBase::set_def_dividerchar(string const& token)
{
    //cout << __func__ << " => " << token << endl;
}
/// @param token BUS bit characters 
void DefDataBase::set_def_busbitchars(string const& token)
{
    //cout << __func__ << " => " << token << endl;
}
/// @param token DEF version 
void DefDataBase::set_def_version(string const& token)
{
    //cout << __func__ << " => " << token << endl;
    this->defVersion = token;
}
/// @param token design name 
void DefDataBase::set_def_design(string const& token)
{
    //cout << __func__ << " => " << token << endl;
    this->defDesign = token;
}
/// @param token DEF unit 
void DefDataBase::set_def_unit(int token)
{
    //cout << __func__ << " => " << token << endl;
    this->defUnit = token;
}
/// @param t1, t2, t3, t4 die area (xl, yl, xh, yh)
void DefDataBase::set_def_diearea(int t1, int t2, int t3, int t4)
{
    //cout << __func__ << " => " << t1 << "," << t2 << "," << t3 << "," << t4 << endl;
    this->dieAreaInUm[0] = t1 * 1.0 / this->defUnit;    // xmin in um
    this->dieAreaInUm[1] = t2 * 1.0 / this->defUnit;    // ymin in um
    this->dieAreaInUm[2] = t3 * 1.0 / this->defUnit;    // xmax in um
    this->dieAreaInUm[3] = t4 * 1.0 / this->defUnit;    // ymax in um
}
/// @brief add row 
void DefDataBase::add_def_row(DefParser::Row const&)
{
    //cout << __func__ << endl;
}
/// @brief add component 
/// @param c component 
void DefDataBase::add_def_component(DefParser::Component const& c)
{
    //cout << __func__ << ": " << c.comp_name << ": status = " << c.status << endl;
    double xInUm = c.origin[0] * 1.0 / this->defUnit;
    double yInUm = c.origin[1] * 1.0 / this->defUnit;
    this->allComponents[c.comp_name] = { xInUm, yInUm, c.orient, c.macro_name };
}
/// @param token number of components 
void DefDataBase::resize_def_component(int token)
{
    //cout << __func__ << " => " << token << endl;
}
/// @brief add pin 
/// @param p pin 
void DefDataBase::add_def_pin(DefParser::Pin const& p)
{
    //cout << __func__ << ": " << p.pin_name << endl;
    double xInUm = p.origin[0] * 1.0 / this->defUnit;
    double yInUm = p.origin[1] * 1.0 / this->defUnit;
    DefPinInfo tempDefPinInfo(p.direct, p.vLayer, xInUm, yInUm);
    this->allDefPins[p.pin_name] = tempDefPinInfo;
}
/// @brief set number of pins 
/// @param token number of pins 
void DefDataBase::resize_def_pin(int token)
{
    //cout << __func__ << " => " << token << endl;
}
/// @brief add net 
/// @param n net 
void DefDataBase::add_def_net(DefParser::Net const& n)
{
    //cout << __func__ << ": " << n.net_name << ": weight " << n.net_weight << endl;
    NetInfo curNet = { n.net_name, n.net_weight, n.vNetPin };
    this->allNets.push_back(curNet);
}
/// @brief set number of nets 
/// @param token number of nets 
void DefDataBase::resize_def_net(int token)
{
    //cout << __func__ << " => " << token << endl;
}
/// @brief set number of blockages 
/// @param n number of blockages 
void DefDataBase::resize_def_blockage(int n)
{
    //cout << __func__ << " => " << n << endl;
}
/// @brief add placement blockages 
/// @param vBbox array of boxes with xl, yl, xh, yh
void DefDataBase::add_def_placement_blockage(std::vector<std::vector<int> > const& vBbox)
{
    /*
    cout << __func__ << " => ";
    for (std::vector<std::vector<int> >::const_iterator it = vBbox.begin(); it != vBbox.end(); ++it)
    cout << "(" << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << ", " << (*it)[3] << ") ";
    cout << endl;
    */
}
/// @brief end of design 
void DefDataBase::end_def_design()
{
    //cout << __func__ << endl;
}


void LefDataBase::appendCellMap(const unordered_map<string, LefCellInfo>& newCells) {
    this->allCells.insert(newCells.begin(), newCells.end());
};


void LefDataBase::print_allCells() {
    cout << "Total " << this->allCells.size() << " cells in all LEF files. (Cells & pins printed below are unordered)\n";
    for (const auto& [cellName, cellInfo] : this->allCells) {
        cout << "\nCell \"" << cellName << "\": origin " << cellInfo.originXInUm << "  " << cellInfo.originYInUm << " um, size "
            << cellInfo.sizeXInUm << "  " << cellInfo.sizeYInUm << " um \n";
        for (const auto& [pinName, pinInfo] : cellInfo.allPinsInCell) {
            cout << "PIN " << pinName << ", " << pinInfo.direct << ", Layer ";
            for (const auto& layer : pinInfo.vLayer) {
                cout << layer << "  ";
            }
            cout << endl;
            for (const auto& rect : pinInfo.vRectsInUm) {
                cout << "    RECT " << rect[0] << "  " << rect[1] << "  " << rect[2] << "  " << rect[3] << endl;
            }
        }
    }
}

//////////////////// required callbacks from abstract LefParser::LefDataBase ///////////////////
/// @brief set LEF version 
/// @param v string of LEF version 
void LefDataBase::lef_version_cbk(string const& v)
{
    //cout << "lef version = " << v << endl;
    this->lefVersion = v;
}
/// @brief set LEF version 
/// @param v floating point number of LEF version 
void LefDataBase::lef_version_cbk(double v)
{
    //cout << "lef version = " << v << endl;
    this->lefVersion = to_string(v);
}
/// @brief set divider characters 
/// @param v divider characters
void LefDataBase::lef_dividerchar_cbk(string const& v)
{
    //cout << "lef dividechar = " << v << endl;
}
/// @brief set unit 
/// @param v an object for unit 
void LefDataBase::lef_units_cbk(LefParser::lefiUnits const& v)
{
    //v.print(stdout);      // typically defined in tech LEF, not in cell LEF
}
/// @brief set manufacturing entry 
/// @param v manufacturing entry 
void LefDataBase::lef_manufacturing_cbk(double v)
{
    //cout << "lef manufacturing = " << v << endl;
}
/// @brief set bus bit characters 
/// @param v but bit characters 
void LefDataBase::lef_busbitchars_cbk(string const& v)
{
    //cout << "lef busbitchars = " << v << endl;
}
/// @brief add layer 
/// @param v an object for layer 
void LefDataBase::lef_layer_cbk(LefParser::lefiLayer const& v)
{
    //v.print(stdout);
}
/// @brief add via 
/// @param v an object for via 
void LefDataBase::lef_via_cbk(LefParser::lefiVia const& v)
{
    //v.print(stdout);
}
/// @brief add via rule 
/// @param v an object for via rule 
void LefDataBase::lef_viarule_cbk(LefParser::lefiViaRule const& v)
{
    //v.print(stdout);
}
/// @brief spacing callback 
/// @param v an object for spacing 
void LefDataBase::lef_spacing_cbk(LefParser::lefiSpacing const& v)
{
    //v.print(stdout);
}
/// @brief site callback 
/// @param v an object for site 
void LefDataBase::lef_site_cbk(LefParser::lefiSite const& v)
{
    //v.print(stdout);
}
/// @brief macro begin callback, describe standard cell type 
/// @param v name of macro 
void LefDataBase::lef_macrobegin_cbk(std::string const& v)
{
    //cout << __func__ << " => " << v << endl;

    // start reading a new cell
    this->tempCellName = v;
    this->tempPinsInCell.clear();   // clear stored pin info in prev cell
}
/// @brief macro callback, describe standard cell type 
/// @param v an object for macro 
void LefDataBase::lef_macro_cbk(LefParser::lefiMacro const& v)
{
    //v.print(stdout);

    // refer to "/limbo/thirdparty/lefdef/5.8/lef/lef/lefiMacro.hpp" for detailed lefiMacro
    double originX = v.originX();
    double originY = v.originY();
    double sizeX = v.sizeX();
    double sizeY = v.sizeY();
    this->allCells[this->tempCellName] = { originX, originY, sizeX, sizeY, this->tempPinsInCell };

    // here, end reading current cell 
}
/// @brief property callback 
/// @param v an object for property 
void LefDataBase::lef_prop_cbk(LefParser::lefiProp const& v)
{
    //v.print(stdout);
}
/// @brief noise margin callback 
/// @param v an object for noise margin 
void LefDataBase::lef_maxstackvia_cbk(LefParser::lefiMaxStackVia const& v)
{
    //v.print(stdout);
}
/// @brief obstruction callback 
/// @param v an object for obstruction 
void LefDataBase::lef_obstruction_cbk(LefParser::lefiObstruction const& v)
{
    //v.print(stdout);
}
/// @brief pin callback, describe pins in a standard cell 
/// @param v an object for pin 
void LefDataBase::lef_pin_cbk(lefiPin const& v)
{
    //v.print(stdout);

    // refer to "/limbo/thirdparty/lefdef/5.8/lef/lef/lefiMacro.hpp" for detailed lefiPin
    string pinName = v.name();
    string pinDirect(v.direction());

#define IGNORE_DEBUG_PRINT_lefiPin_
#ifndef IGNORE_DEBUG_PRINT_lefiPin_
    cout << "MACRO " << this->tempCellName << ": Pin " << pinName << "  ";
#endif

    // find geometries "Layer" and "Rect" of this pin, details in "5.8/lef/lef/lefiMisc.hpp"
    vector<string> vLayer = {};
    vector<vector<double>> vRect;

    lefiGeometries  *geo = v.port(0);       // only consider the first port of each pin, ignore other ports of this pin
    int numItems = geo->numItems();         // num of items like Layer, Rect, etc.
    lefiGeomRect *rect;
    for (int i = 0; i < numItems; i++) {
        switch (geo->itemType(i)) {

        case lefiGeomLayerE:
            vLayer.push_back((string)geo->getLayer(i));
#ifndef IGNORE_DEBUG_PRINT_lefiPin_
            fprintf(stdout, "Layer %s\n", geo->getLayer(i));
#endif
            break;
        case lefiGeomRectE:
            rect = geo->getRect(i);
            vRect.push_back({ rect->xl, rect->yl, rect->xh, rect->yh });
#ifndef IGNORE_DEBUG_PRINT_lefiPin_
            if (rect->colorMask) {
                fprintf(stdout, "Rect MASK %d, %g,%g  %g,%g\n", rect->colorMask,
                    rect->xl, rect->yl,
                    rect->xh, rect->yh);
            }
            else {
                fprintf(stdout, "Rect %g,%g  %g,%g\n", rect->xl, rect->yl,
                    rect->xh, rect->yh);
            }
#endif
            break;
        default:
#ifndef IGNORE_DEBUG_PRINT_lefiPin_
            lefiError(0, 1375, "ERROR (LEFPARS-1375): unknown geometry type");
            fprintf(stdout, "Unknown geometry type %d\n",
                (int)(geo->itemType(i)));
#endif
            break;
        }
    }

    LefPinInfo tempLefPinInfo(pinDirect, vLayer, vRect);
    this->tempPinsInCell[pinName] = tempLefPinInfo;
}

unordered_map<string, LayerMapInfo> readLayerMap(string inFileName) {
    /*  return map[layerName] = struct LayerMapInfo */
    ifstream inFile(inFileName);
    if (!inFile.is_open()) {
        cerr << "ERROR in opening file " << inFileName << endl;
        exit(-1);
    }

    unordered_map<string, LayerMapInfo> layerMap;
    
    string curLine;
    string firstWord;
    int layerNameInNum = -1;
    double zmin = 0, zmax = 0;
    double lengthUnit = 1.0;

    bool notFoundLayerMapYet = true;
    while (notFoundLayerMapYet && getline(inFile, curLine)) {
        stringstream curLineStream(curLine);
        curLineStream >> firstWord;

        if (firstWord == "lengthUnit") {
            string temp;
            curLineStream >> temp;
            curLineStream >> lengthUnit;
        }
        if (firstWord == "BEGIN_LAYER_MAP") {
            notFoundLayerMapYet = false;
            while (getline(inFile, curLine)) {  // read the layer map section
                stringstream ls(curLine);
                ls >> firstWord;
                if (firstWord == "END_LAYER_MAP") break;

                ls >> layerNameInNum >> zmin >> zmax;
                
                double zminInUm = zmin * lengthUnit * 1.0e6;
                double zmaxInUm = zmax * lengthUnit * 1.0e6;
                layerMap[firstWord] = { layerNameInNum, zminInUm, zmaxInUm };
            }
        }
    }
    return layerMap;
}

void print_layerMap(const unordered_map<string, LayerMapInfo>& layerMap) {
    cout << "Layer Map: (layerName layerNameInNumUsedInGDS zminInUm zmaxInUm) \n";
    for (auto&[layerName, layerInfo] : layerMap) {
        cout << layerName << "  " << layerInfo.layerNameInNum << "  "
            << layerInfo.zminInUm << "  " << layerInfo.zmaxInUm << endl;
    }
}


AutoPorts::AutoPorts() {
}

void AutoPorts::readLayerMap_cbk(string inFileName) {
    this->layerMap = readLayerMap(inFileName);
}

void AutoPorts::print_layerMap_cbk() {
    print_layerMap(this->layerMap);
}

// Transform local coordinate within a LEF cell to global coordinate at entire DEF design structure
vector<double> localLefCoorToGlobalDefCoor(
    double localLefCoor[2],     // local LEF coordinate
    double sizeBB[2],           // size of bounding box of the cell
    double origin[2],           // origin of the cell to align with a DEF COMPONENT placement point
    double placementDef[2],     // the DEF COMPONENT placement point
    string orient)              // orientation of the cell defined in DEF COMPONENTS
{
    // step 1: find local coordinate {xc, yc} relative to cell center
    double xc = localLefCoor[0] - sizeBB[0] / 2;
    double yc = localLefCoor[1] - sizeBB[1] / 2;

    // step 2: mirror & rotation with respect to cell center. => {xmr, ymr}
    double xmr, ymr;
    double xll = -sizeBB[0] / 2;// new lower left corner w.r.t. cell center
    double yll = -sizeBB[1] / 2;
    if (orient == "N") {        // R0
        xmr = xc; 
        ymr = yc;
    }
    else if (orient == "S") {   // R180
        xmr = -xc;
        ymr = -yc;
    }
    else if (orient == "FN") {  // MY
        xmr = -xc;
        ymr = yc;
    }
    else if (orient == "FS") {  // MX
        xmr = xc;
        ymr = -yc;
    }
    else {
        // new lower left corner w.r.t. cell center when orient = W, E, FW, FE
        xll = -sizeBB[1] / 2;
        yll = -sizeBB[0] / 2;

        if (orient == "W") {        // R90
            xmr = -yc;
            ymr = xc;
        }
        else if (orient == "E") {   // R270
            xmr = yc;
            ymr = -xc;
        }
        else if (orient == "FW") {  // MX90
            xmr = yc;
            ymr = xc;
        }
        else if (orient == "FE") {  // MY90
            xmr = -yc;
            ymr = -xc;
        }
        else {
            cerr << "ERROR! Invalid orientation type: " << orient
                << ". Reset to default orient \"N\"." << endl;
            xmr = xc;
            ymr = yc;
            double xll = -sizeBB[0] / 2;
            double yll = -sizeBB[1] / 2;
        }
    }

    // step 3: find global DEF coordinate {xg, yg}
    double xg = (xmr - xll) - origin[0] + placementDef[0];
    double yg = (ymr - yll) - origin[1] + placementDef[1];
    return {xg, yg};
}

void AutoPorts::getPortCoordinate(
    const unordered_map<string, ComponentInfo>& allComponentsDEF,
    const unordered_map<string, DefPinInfo>& allDefPinsDEF,
    const vector<NetInfo>& allNetsDEF,
    const unordered_map<string, LefCellInfo>& allCellsLEF)
{
    for (const auto& net : allNetsDEF) {                        // loop over all nets
        const string& netName = net.netName;
        vector<NetPortCoor> vPortCoor;
        for (const auto& NodenamePin : net.vNodenamePin) {      // loop over all nodes in cur net
            const string& nodeName = NodenamePin.first;
            const string& pinName = NodenamePin.second;

            // Node type 1: defined as component
            if (nodeName != "PIN") {
                // corresponding cell of cur component
                const ComponentInfo& defComp = allComponentsDEF.at(nodeName);
                string cellName = defComp.cellName;
                const LefCellInfo& lefCell = allCellsLEF.at(cellName);

                // cur cell info: size of bounding box, origin, placement, and orient
                double sizeBBInUm[2] = { lefCell.sizeXInUm, lefCell.sizeYInUm };
                double originInUm[2] = { lefCell.originXInUm, lefCell.originYInUm };
                double placementDefInUm[2] = { defComp.xInUm, defComp.yInUm };
                string cellOrient = defComp.orient;

                // pin of cur cell and its local coordinate within the cell
                const LefPinInfo& lefPin = lefCell.allPinsInCell.at(pinName);
                const vector<double>& rect = lefPin.vRectsInUm[0];      // use the first Rect of cur pin
                double localCoor_xInUm = (rect[0] + rect[2]) / 2.0;     // use the center point of 1st rect as cur pin location
                double localCoor_yInUm = (rect[1] + rect[3]) / 2.0;

                // apply transformation from local cell coordinate to global port coordinate
                double localLefCoorInUm[2] = { localCoor_xInUm, localCoor_yInUm };
                vector<double> globalCoorInUm = localLefCoorToGlobalDefCoor(
                    localLefCoorInUm,   // local LEF coordinate
                    sizeBBInUm,         // size of bounding box of the cell
                    originInUm,         // origin of the cell to align with a DEF COMPONENT placement point
                    placementDefInUm,   // the DEF COMPONENT placement point
                    cellOrient);        // orientation of the cell defined in DEF COMPONENTS

                // organize port coordinate together in struct NetPortCoor
                NetPortCoor portCoor;
                portCoor.direct = lefPin.direct;
                portCoor.vLayer = lefPin.vLayer;
                portCoor.portName = "port_" + netName + "_" + nodeName + "_" + pinName;
                portCoor.xInUm = globalCoorInUm[0];
                portCoor.yInUm = globalCoorInUm[1];
                string layerName = lefPin.vLayer[0];                    // use the first found layerName as the port layer
                portCoor.zInUm = this->layerMap[layerName].zminInUm;    // put port at layer bottom
                portCoor.gdsiiNum = this->layerMap[layerName].layerNameInNum;

                vPortCoor.push_back(portCoor);
            }

            // Node type 2: defined as external pin
            if (nodeName == "PIN") {
                const DefPinInfo& defPin = allDefPinsDEF.at(pinName);
                NetPortCoor portCoor;
                portCoor.direct = defPin.direct;
                portCoor.vLayer = defPin.vLayer;
                portCoor.portName = "port_" + netName + "_PIN_" + pinName;
                portCoor.xInUm = defPin.xInUm;
                portCoor.yInUm = defPin.yInUm;
                string layerName = defPin.vLayer[0];                    // use the first found layerName as the port layer
                portCoor.zInUm = this->layerMap[layerName].zminInUm;    // put port at layer bottom
                portCoor.gdsiiNum = this->layerMap[layerName].layerNameInNum;

                vPortCoor.push_back(portCoor);
            }
        }
        this->netName_to_vPortCoor[netName] = vPortCoor;
    }

}

void AutoPorts::print_netName_to_vPortCoor() {
    cout << "Port Coordinate in um: (portName layerName, layerNameInNum, x, y, z) \n";
    for (auto&[netName, vPortCoor] : this->netName_to_vPortCoor) {
        cout << "\nNet \"" << netName << "\": \n"; 
        for (auto & portCoor : vPortCoor) {
            cout << portCoor.portName << "  " << portCoor.vLayer[0] << "  " << portCoor.gdsiiNum << "  "
                << portCoor.xInUm << "  " << portCoor.yInUm << "  " << portCoor.zInUm << endl;
        }
    }
}

