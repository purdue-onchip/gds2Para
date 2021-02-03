#include "autoPortFromDefLef.hpp"



void DefDataBase::print_allNets() {
    cout << "Total " << this->allNets.size() << " nets. netName: (nodeName, pin) \n";
    for (const auto& net : this->allNets) {
        cout << net.netName << ": ";
        for (const auto& NodenamePin : net.vNodenamePin)
            cout << "(" << NodenamePin.first << ", " << NodenamePin.second << ") ";
        cout << endl;
    }
}

void DefDataBase::print_allPins() {
    cout << "Total " << this->allPins.size() << " external pins. pinName (unordered): x(um), y(um), direction, layers\n";
    for (const auto& pin : this->allPins) {
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
            if ((nodeName == "PIN") && (this->allPins.find(pinName) == this->allPins.end())) {
                isValidCom = false;
                cout << "Net \"" << net.netName << "\" has invalid node name at pin \"" << nodeName << "\"\n";
            }
        }
    }

    return isValidCom;
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
    this->allPins[p.pin_name] = { xInUm, yInUm, p.direct, p.vLayer };
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