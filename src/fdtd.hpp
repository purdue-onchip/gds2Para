#ifndef _FDTD_H
#define _FDTD_H
//#include "stdafx.h"
#define _USE_MATH_DEFINES // Place before including <cmath> for e, log2(e), log10(e), ln(2), ln(10), pi, pi/2, pi/4, 1/pi, 2/pi, 2/sqrt(pi), sqrt(2), and 1/sqrt(2)
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <queue>
#include <stack>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <utility>

// HYPRE and MKL data type control
#define LARGE_SYSTEM (1) // Must leave defined with the current makefile
#ifdef LARGE_SYSTEM
#define MKL_ILP64 (1) // Must define before including mkl.h if using long long int 
typedef long long int myint;
#else
typedef int myint;
#endif
#include <mkl.h>
#include <mkl_spblas.h>

// Manipulate namespace
using namespace std;

// Fundamental physical constant macros
#define MU (4*M_PI*1.e-7)
#define CSPED (299792458.)
#define EPSILON0 (1./(CSPED*CSPED*MU))

// Solver discretization control macros
#define SIGMA (5.8e+7) // Default conductivity for conductors is copper (S/m)
#define DOUBLEMAX (1.e+30)
#define DOUBLEMIN (-1.e+30)
#define MINDISFRACX (3e-3) // Fraction setting minimum discretization retained in x-directions after node merging in terms of smaller of x-extent
#define MINDISFRACY (2.5e-3) // Fraction setting minimum discretization retained in y-directions after node merging in terms of smaller of y-extent
#define MINDISFRACZ (0.05) // Fraction setting minimum discretization retained in z-direction after node merging in terms of distance between closest layers
#define MAXDISFRACX (0.1) // Fraction setting largest discretization in x-direction in terms of x-extent
#define MAXDISFRACY (0.1) // Fraction setting largest discretization in y-direction in terms of y-extent
#define MAXDISLAYERZ (2.) // Largest discretization in z-direction represented as fewest nodes placed between closest layers (1. = distance between closest layers, 2. = half distance between closest layers)
#define DT (1.e-15) // Time step for finding high-frequency modes (s)

// Debug testing macros (comment out if not necessary)
//#define UPPER_BOUNDARY_PEC
//#define LOWER_BOUNDARY_PEC
#define PRINT_NODE_COORD
#define PRINT_DIS_COUNT (1)
#define SKIP_MARK_CELL
#define PRINT_VERBOSE_TIMING // Terminal output has extra runtime clock information
//#define PRINT_PORT_SET
//#define PRINT_V0D_BLOCKS
//#define PRINT_V0_Z_PARAM
//#define PRINT_V0_Vh_Z_PARAM
#define SKIP_PARDISO // Remove PARDISO solver code
#define GENERATE_V0_SOLUTION
//#define SKIP_VH
//#define SKIP_GENERATE_STIFF
#define SKIP_STIFF_REFERENCE 

// Disable layered FDTD code (comment out if you want to test layered FDTD)
#define SKIP_WRITE_SYS_TO_FILE        // Skip writing sys obj to txt files
#define SKIP_LAYERED_FDTD             // Skip the main function to call layeredFdtd code

// Function-like macros
#define NELEMENT(x) ((sizeof x) / (sizeof x[0]))


class fdtdOneCondct {
public:
    int cdtName;
    int vcc;
    int plaStackId;
    int portId;
    myint numNode;
    int numVert;
    double *x;
    double *y;
    double *z;
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double sigma;
    myint *cdtInNode;
    int *markCdtInNode;
    int layer;
};

class fdtdCdt {
public:
    myint *node;
    int markPort;
    int *portNode;
    int portind;
    myint cdtNodeind;
};

class fdtdBound {
public:
    int numBound;
    int numPeri; /*number of open surrounding boundaries */
    int topNum; /*top open boundary index*/
    int botNum; /*bot open boundary index*/
    int *p1;
    int *p2;
    int *bound;
    int *ukwnBegEdge; /*starting edge of each bound*/
    int *ukwnBegNode; /*starting node of each bound*/
    int *ukwnNumEdge; /*number of edges of each bound*/
    int *ukwnNumNode;  /*number of nodes of each bound*/
    int top;
    int bot;
    double **norm;/*bound norm direction*/
    int *numPort;
    int **port;
};

class fdtdPort{
    /* Port information */
public:
    double *x, *y, *z;
    int multiplicity;                        // Number of port sides to excite
    vector<double> x1, x2, y1, y2, z1, z2;   // Supply/return coordinates for each port side
    vector<myint> portCnd;                   // Conductor index of supply point for each port side
    vector<vector<myint>> portEdge;          // Edge indices from return to supply with nonzero current density for each port side
    vector<double> portArea;                 // Area of each port side (m^2)
    vector<int> portDirection;               // Sign of current density component for each port side
    int *node;

    /* Default Constructor */
    fdtdPort() {
        this->x = NULL;
        this->y = NULL;
        this->z = NULL;
        this->multiplicity = 0;
        this->x1 = {};
        this->y1 = {};
        this->z1 = {};
        this->x2 = {};
        this->y2 = {};
        this->z2 = {};
        this->portCnd = {};
        this->portEdge = {};
        this->portArea = {};
        this->portDirection = {};
        this->node = NULL;
    }

    /* Print Function */
    void print();

    /* Destructor */
    ~fdtdPort() {
        free(this->x);
        free(this->y);
        free(this->z);
        free(this->node);
    }
};

class fdtdPatch {
public:
    double patchArea;
    int *node;
};

class fdtdMesh {
    /* Mesh information */
public:
    double lengthUnit;
    myint outedge;
    myint inedge;

    /* Frequency parameters */
    double freqUnit;
    double freqStart;
    double freqEnd;
    int nfreq;
    int freqScale;

    /* Discretization information*/
    myint nx, ny, nz;    // number of nodes along x, y, z

    double *xn, *yn, *zn;    // coordinates of the nodes along x, y, z
    double *xnu, *ynu, *znu;

    myint N_cell_x;
    myint N_cell_y;
    myint N_cell_z;

    double xlim1, xlim2;
    double ylim1, ylim2;
    double zlim1, zlim2;

    myint N_edge;
    myint N_edge_s;
    myint N_edge_v;

    myint N_node;
    myint N_node_s;

    myint N_patch;
    myint N_patch_s;
    myint N_patch_v;

    //double *nodepos;
    double *Epoints;
    //myint *edgelink;
    double *Hpoints;
    vector<vector<pair<myint, double>>> nodeEdge;    // for each node which edge is connected with this node
    vector<vector<pair<myint, double>>> nodeEdgea;    // for each node which edge is connected with this node

    /* Upper and lower PEC */
    int *bd_node1;   //lower PEC
    int *bd_node2;   //upper PEC
    int *bd_edge;
    set<myint> ubde, lbde;    // upper boundary edge and lower boundary edge
    set<myint> ubdn, lbdn;    // upper boundary node and lower boundary node
    myint* mapEdge;   // map the original edges to the new edge # with upper and lower PEC boundaries
    myint* mapEdgeR;    // map the new edge # to the original edges
    int bden;    // boundary edge number

    /* Layer stack up parameters */
    int numStack;
    vector<double> stackEps;
    vector<double> stackSig;
    vector<double> stackBegCoor;
    vector<double> stackEndCoor;
    vector<string> stackName;
    vector<double> eps;
    vector<double> stackEpsn;
    vector<double> stackSign;
    double *stackCdtMark;

    /* Conductor parameters */
    vector<fdtdOneCondct> conductorIn;
    myint numCdtRow;                      // how many input rows
    myint numCdt;                         // number of isolated conductors in design
    myint *markEdge;                      // mark if this edge is inside a conductor
    myint *markCell;
    myint *cdtNumNode;                    // number of nodes along each isolated conductor
    double *sig;
    fdtdCdt *conductor;                   // information about isolated conductors
    myint *markNode;                      // mark this node if it is inside the conductor
    vector<vector<int>> edgeCell;         // for each cell which edge is around it
    vector<vector<double>> edgeCellArea;  // for each cell the area of the perpendicular rectangle
    vector<int> acu_cnno;                 // accumulated conductor number of nodes
    vector<int> cindex;
    int *exciteCdtLayer;
    unordered_set<int> cond2condIn;       // put the active conductors' corresponding conductorIn #
    vector<bool> markProSide;             // mark nodes near excited conductors for less aggressive node merging in rapid field change area

    /* Patch information */
    fdtdPatch *patch;

    /* Boundary information */
    fdtdBound *bound;

    /* V0c row, column, value */
    myint *v0cRowId;
    myint *v0cColId;
    myint *v0cColIdo;
    double *v0cval;
    double *v0cvalo;

    myint *v0caRowId;
    myint *v0caColId;
    myint *v0caColIdo;
    double *v0caval;
    double *v0cavalo;

    /* V0c2 and yc row, column, value */
    double *v0c2y0c2;
    double *v0c2y0c2o;
    double *yc;
    double *v0cy0c;

    /* V0c'*D_sig*V0c row, column, value */
    myint *AcRowId;
    myint *AcRowId1;
    myint *AcColId;
    double *Acval;
    myint *AdRowId;
    myint *AdRowId1;
    myint *AdColId;
    double *Adval;

    double *crhs;

    /* V0d row, column, value */
    myint *v0d1RowId;
    myint *v0d1ColId;
    myint *v0d1ColIdo;
    double *v0d1val;
    double *v0d1valo;

    myint *v0d1aRowId;
    myint *v0d1aColId;
    myint *v0d1aColIdo;
    double *v0d1aval;
    double *v0d1avalo;

    myint *v0d2RowId;
    myint *v0d2ColId;
    myint *v0d2ColIdo;
    double *v0d2val;
    double *v0d2valo;

    myint *v0d2aRowId;
    myint *v0d2aColId;
    myint *v0d2aColIdo;
    double *v0d2aval;
    double *v0d2avalo;

    double *yd;

    /* Vh */
    lapack_complex_double *Vh;
    myint leng_Vh;

    /* Se and Sh */
    myint *SRowId;
    myint *SColId;
    double *Sval;
    myint leng_S;

	/* Laplacian decomposed v0 xh left matrix */
	myint* LdrowId;
	myint* LdcolId;
	double* Ldval;
	myint leng_Ld;

    /* Solution storage */
    complex<double> *y;
    vector<complex<double>> x;    // the solution involving all the sourcePorts

    /* Port information */
    int numPorts;
    vector<fdtdPort> portCoor;

    /* Current sources */
    double *J;

    /* Current V0c,s^T*I matrix */
    complex<double> *v0csJ;
    complex<double> *Y;

    /* Default Constructor */
    fdtdMesh(){
        // Set some important numbers to zero
        this->outedge = (myint)0;
        this->inedge = (myint)0;
        this->numCdtRow = (myint)0;
        this->leng_Vh = (myint)0;
        this->leng_S = (myint)0;

        // Set all pointers to NULL for safety
        this->xn = NULL;
        this->yn = NULL;
        this->zn = NULL;
        this->xnu = NULL;
        this->ynu = NULL;
        this->znu = NULL;
        //this->nodepos = NULL;
        this->Epoints = NULL;
        //this->edgelink = NULL;
        this->Hpoints = NULL;
        this->bd_node1 = NULL;
        this->bd_node2 = NULL;
        this->bd_edge = NULL;
        this->stackCdtMark = NULL;
        this->markEdge = NULL;
        this->markCell = NULL;
        this->cdtNumNode = NULL;
        this->sig = NULL;
        this->conductor = NULL;
        this->markNode = NULL;
        this->exciteCdtLayer = NULL;
        this->patch = NULL;
        this->bound = NULL;
        this->v0cRowId = NULL;
        this->v0cColId = NULL;
        this->v0cColIdo = NULL;
        this->v0cval = NULL;
        this->v0cvalo = NULL;
        this->v0caRowId = NULL;
        this->v0caColId = NULL;
        this->v0caColIdo = NULL;
        this->v0caval = NULL;
        this->v0cavalo = NULL;
        this->v0c2y0c2 = NULL;
        this->v0c2y0c2o = NULL;
        this->yc = NULL;
        this->v0cy0c = NULL;
        this->AcRowId = NULL;
        this->AcRowId1 = NULL;
        this->AcColId = NULL;
        this->Acval = NULL;
        this->AdRowId = NULL;
        this->AdRowId1 = NULL;
        this->AdColId = NULL;
        this->Adval = NULL;
        this->crhs = NULL;
        this->v0d1RowId = NULL;
        this->v0d1ColId = NULL;
        this->v0d1ColIdo = NULL;
        this->v0d1val = NULL;
        this->v0d1valo = NULL;
        this->v0d1aRowId = NULL;
        this->v0d1aColId = NULL;
        this->v0d1aColIdo = NULL;
        this->v0d1aval = NULL;
        this->v0d1avalo = NULL;
        this->v0d2RowId = NULL;
        this->v0d2ColId = NULL;
        this->v0d2ColIdo = NULL;
        this->v0d2val = NULL;
        this->v0d2valo = NULL;
        this->v0d2aRowId = NULL;
        this->v0d2aColId = NULL;
        this->v0d2aColIdo = NULL;
        this->v0d2aval = NULL;
        this->v0d2avalo = NULL;
        this->yd = NULL;
        this->Vh = NULL;
        this->SRowId = NULL;
        this->SColId = NULL;
        this->Sval = NULL;
        this->y = NULL;
        this->J = NULL;
        this->v0csJ = NULL;
        this->Y = NULL;

        // Set all vectors to empty vectors
        this->nodeEdge = {};
        this->nodeEdgea = {};
        this->numStack = 0;
        this->stackEps = {};
        this->stackSig = {};
        this->stackBegCoor = {};
        this->stackEndCoor = {};
        this->stackName = {};
        this->eps = {};
        this->stackEpsn = {};
        this->stackSign = {};
        this->conductorIn = {};
        this->edgeCell = {};
        this->edgeCellArea = {};
        this->acu_cnno = {};
        this->cindex = {};
        this->markProSide = {};
        this->x = {};
        this->numPorts = 0;
        this->portCoor = {};

        // Set all other collection data types as empty
        this->cond2condIn = unordered_set<int>();
    }

    /* Print Function */
    void print();

    /* Print conductorIn */
    void printConductorIn(){
        int i, j;

        cout << "Print conductorIn information: " << endl;
        for (i = 0; i < this->numCdtRow; i++){
            for (j = 0; j < this->conductorIn[i].numVert - 1; j++) {
                if (this->conductorIn[i].layer == 5){
                    cout << this->conductorIn[i].x[j] << " " << this->conductorIn[i].y[j] << " " << this->conductorIn[i].zmin << " " << this->conductorIn[i].zmax << endl;
                }
            }
        }
    }

    /* Check whether the point's markNode */
    void checkPoint(double x, double y, double z, unordered_map<double, int> & xi, unordered_map<double, int> & yi, unordered_map<double, int> & zi){
        int inx, iny, inz;

        inx = xi[x];
        iny = yi[y];
        inz = zi[z];
        cout << "inx " << inx << " iny " << iny << " inz " << inz << endl;
        cout << "x " << this->xn[inx] << " y " << this->yn[iny] << " z " << this->zn[inz] << endl;
        cout << "This node's mark is " << this->markNode[inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny] << endl;
    }

    /* Is point (x,y) within the polygon? */
    bool polyIn(double x, double y, int inPoly) {
        int npol;
        myint indi = 0, indj = 0;
        bool isCond = false;
        double disMin = 1.e-10;

        npol = this->conductorIn[inPoly].numVert;

        for (indi = 0, indj = npol - 1; indi < npol; indj = indi++) {
            if ((abs(y - this->conductorIn[inPoly].y[indj]) < disMin && abs(y - this->conductorIn[inPoly].y[indi]) < disMin &&
                ((x >= this->conductorIn[inPoly].x[indj] && x <= this->conductorIn[inPoly].x[indi]) ||
                (x >= this->conductorIn[inPoly].x[indi] && x <= this->conductorIn[inPoly].x[indj])))) {    // on x direction edge
                return true;
            }
            else if (abs(x - this->conductorIn[inPoly].x[indj]) < disMin && abs(x - this->conductorIn[inPoly].x[indi]) < disMin &&
                ((y >= this->conductorIn[inPoly].y[indj] && y <= this->conductorIn[inPoly].y[indi]) ||
                (y >= this->conductorIn[inPoly].y[indi] && y <= this->conductorIn[inPoly].y[indj]))) {    // on y direction edge
                return true;
            }

            else if ((abs(this->conductorIn[inPoly].y[indi] - this->conductorIn[inPoly].y[indj]) > disMin &&
                (((this->conductorIn[inPoly].y[indi] <= y) && (y < this->conductorIn[inPoly].y[indj])) ||
                ((this->conductorIn[inPoly].y[indj] <= y) && (y < this->conductorIn[inPoly].y[indi])))) &&
                (x < (this->conductorIn[inPoly].x[indj] - this->conductorIn[inPoly].x[indi]) * (y - this->conductorIn[inPoly].y[indi]) /
                (this->conductorIn[inPoly].y[indj] - this->conductorIn[inPoly].y[indi]) + this->conductorIn[inPoly].x[indi])) {
                isCond = !isCond;
            }
        }
        return isCond;
    }


    /* Find nodes inside conductors by judging whether the point in each small window is inside the polygon or not */
    void findInsideCond(unordered_map<double, int> & xi, unordered_map<double, int> & yi, unordered_map<double, int> & zi){
        int indi, indj, indl, indk;
        myint xrange_max;
        unordered_map<myint, myint> xrange;
        vector<myint> xcoorv;
        set<myint> xcoor;
        unordered_map<int, vector<double>> xcoory;    // store the start coordinate of the edge mapping the y min and y max
        myint ss, ee;
        myint x1, x2;
        double y1, y2;
        int y1in, y2in;
        int ymax, ymin;
        double y;
        int mark1, mini_k, mark;
        double mini;
        myint node1, node2;



        for (indi = 0; indi < this->numCdtRow; indi++) {
            if (this->conductorIn[indi].zmax == this->conductorIn[indi].zmin)
                continue;
            xrange.clear();
            xcoor.clear();
            xcoorv.clear();
            xcoory.clear();
            mark1 = 0;
            ymax = 0;    // this polygon's max y index
            ymin = this->ny - 1;    // this polygon's min y index
            //cout << endl;
            //cout << endl;
            for (indj = 0; indj < this->conductorIn[indi].numVert - 1; indj++) {
                if (max(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[indj + 1]]) > ymax) {
                    ymax = max(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[indj + 1]]);
                }
                if (min(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[indj + 1]]) < ymin) {
                    ymin = min(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[indj + 1]]);
                }
                if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[indj + 1]) {
                    x1 = xi[this->conductorIn[indi].x[indj]];   // smaller x
                    x2 = xi[this->conductorIn[indi].x[indj + 1]];    // larger x
                }
                else {
                    x1 = xi[this->conductorIn[indi].x[indj + 1]];
                    x2 = xi[this->conductorIn[indi].x[indj]];
                }

                for (indl = x1; indl <= x2; indl++) {
                    xcoor.insert(indl);
                }

            }
            if (max(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[0]]) > ymax) {
                ymax = max(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[0]]);
            }
            if (min(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[0]]) < ymin) {
                ymin = min(yi[this->conductorIn[indi].y[indj]], yi[this->conductorIn[indi].y[0]]);
            }
            if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[0]) {
                x1 = xi[this->conductorIn[indi].x[indj]];   // smaller x
                x2 = xi[this->conductorIn[indi].x[0]];    // larger x
            }
            else {
                x1 = xi[this->conductorIn[indi].x[0]];
                x2 = xi[this->conductorIn[indi].x[indj]];
            }
            for (indl = x1; indl <= x2; indl++) {
                xcoor.insert(indl);
            }

            for (auto xcoori : xcoor) {
                xcoorv.push_back(xcoori);
            }
            xrange_max = xcoorv.back();   // the maximal x coordinate

            for (indj = 0; indj < xcoorv.size() - 1; indj++) {
                mark1 = 1;    // the x coordinates are more than 1
                xrange[xcoorv[indj]] = xcoorv[indj + 1];
            }

            if (xcoorv.size() == 1) {    // If it has only one value
                xrange[xcoorv[0]] = xcoorv[0];
            }
            for (indj = 0; indj < this->conductorIn[indi].numVert - 1; indj++) {
                if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[indj + 1]) {
                    ss = xi[this->conductorIn[indi].x[indj]];
                    ee = xi[this->conductorIn[indi].x[indj + 1]];
                    if (ss == ee && mark1 == 1) {
                        continue;
                    }

                    while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                        if (xcoory.find(ss) == xcoory.end()) {
                            xcoory[ss].push_back(DOUBLEMAX);
                            xcoory[ss].push_back(DOUBLEMIN);
                        }

                        y1 = (this->conductorIn[indi].x[indj + 1] - this->xn[ss]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[ss] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj + 1];
                        y2 = (this->conductorIn[indi].x[indj + 1] - this->xn[xrange[ss]]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[xrange[ss]] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj + 1];
                        if (xcoory[ss][0] > min(y1, y2)){
                            xcoory[ss][0] = min(y1, y2);
                        }
                        if (xcoory[ss][1] < max(y1, y2)){
                            xcoory[ss][1] = max(y1, y2);
                        }
                        if (mark1 == 0) {   // xrange only has one value, break
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
                else if (this->conductorIn[indi].x[indj] > this->conductorIn[indi].x[indj + 1]) {
                    ss = xi[this->conductorIn[indi].x[indj + 1]];
                    ee = xi[this->conductorIn[indi].x[indj]];
                    if (ss == ee && mark1 == 1) {
                        continue;
                    }
                    while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                        if (xcoory.find(ss) == xcoory.end()) {
                            xcoory[ss].push_back(DOUBLEMAX);
                            xcoory[ss].push_back(DOUBLEMIN);
                        }

                        y1 = (this->conductorIn[indi].x[indj + 1] - this->xn[ss]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[ss] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj + 1];
                        y2 = (this->conductorIn[indi].x[indj + 1] - this->xn[xrange[ss]]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[xrange[ss]] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj + 1];
                        if (xcoory[ss][0] > min(y1, y2)){
                            xcoory[ss][0] = min(y1, y2);
                        }
                        if (xcoory[ss][1] < max(y1, y2)){
                            xcoory[ss][1] = max(y1, y2);
                        }
                        if (mark1 == 0) {    // xrange only has one value, break
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
            }
            if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[0]) {
                ss = xi[this->conductorIn[indi].x[indj]];
                ee = xi[this->conductorIn[indi].x[0]];
                if (!(ss == ee && mark1 == 1)) {
                    

                    while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                        if (xcoory.find(ss) == xcoory.end()) {
                            xcoory[ss].push_back(DOUBLEMAX);
                            xcoory[ss].push_back(DOUBLEMIN);
                        }

                        y1 = (this->conductorIn[indi].x[0] - this->xn[ss]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[ss] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[0];
                        y2 = (this->conductorIn[indi].x[0] - this->xn[xrange[ss]]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[xrange[ss]] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[0];
                        if (xcoory[ss][0] > min(y1, y2)){
                            xcoory[ss][0] = min(y1, y2);
                        }
                        if (xcoory[ss][1] < max(y1, y2)){
                            xcoory[ss][1] = max(y1, y2);
                        }
                        if (mark1 == 0) {   // xrange only has one value, break
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
            }
            else if (this->conductorIn[indi].x[indj] > this->conductorIn[indi].x[0]) {
                ss = xi[this->conductorIn[indi].x[0]];
                ee = xi[this->conductorIn[indi].x[indj]];
                if (!(ss == ee && mark1 == 1)) {
                    
                    while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                        if (xcoory.find(ss) == xcoory.end()) {
                            xcoory[ss].push_back(DOUBLEMAX);
                            xcoory[ss].push_back(DOUBLEMIN);
                        }

                        y1 = (this->conductorIn[indi].x[0] - this->xn[ss]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[ss] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[0];
                        y2 = (this->conductorIn[indi].x[0] - this->xn[xrange[ss]]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[xrange[ss]] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[0];
                        if (xcoory[ss][0] > min(y1, y2)){
                            xcoory[ss][0] = min(y1, y2);
                        }
                        if (xcoory[ss][1] < max(y1, y2)){
                            xcoory[ss][1] = max(y1, y2);
                        }
                        if (mark1 == 0) {    // xrange only has one value, break
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
            }

            if (mark1 == 0) {    // only has one x
                for (auto xcooryi : xcoory){
                    x1 = xcooryi.first;    // x index for this window
                    y1 = xcooryi.second[0];    // this window's smallest y
                    y2 = xcooryi.second[1];    // this window's largest y
                    y1in = ymin;    // this window's smallest y's index
                    y2in = ymax;    // this window's largest y's index
                    for (indj = ymin; indj <= ymax; indj++) {
                        if (this->yn[indj] >= y2) {
                            y2in = indj;
                            break;
                        }
                        if (this->yn[indj] < y1) {
                            y1in = indj;
                        }
                    }
                    for (indj = y1in; indj <= y2in; indj++) {
                        if (polyIn(this->xn[x1], this->yn[indj], indi)) {   // if this point is inside the polygon, set the z direction edge in markEdge with 1 and markNode with 1
                            for (indl = zi[this->conductorIn[indi].zmin]; indl < zi[this->conductorIn[indi].zmax]; indl++) {
                                this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + x1 * (this->N_cell_y + 1) + indj] = indi + 1;
                                compute_edgelink(indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + x1 * (this->N_cell_y + 1) + indj, node1, node2);
                                this->markNode[node1] = 1;   // if the edge is inside the conductor, then both its two end nodes are inside the conductor
                                this->markNode[node2] = 1;
                            }
                        }
                        if (indj < y2in) {
                            if (polyIn(this->xn[x1], (this->yn[indj] + this->yn[indj + 1]) / 2, indi)) {   // if the edge along y axis is inside the ploygon, set the y direction edge with 1 in markEdge
                                for (indl = zi[this->conductorIn[indi].zmin]; indl <= zi[this->conductorIn[indi].zmax]; indl++) {
                                    this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + x1 * (this->N_cell_y) + indj] = indi + 1;
                                    compute_edgelink(indl * (this->N_edge_s + this->N_edge_v) + x1 * (this->N_cell_y) + indj, node1, node2);
                                    this->markNode[node1] = 1;   // if the edge is inside the conductor, then both its two end nodes are inside the conductor
                                    this->markNode[node2] = 1;
                                }
                            }
                        }
                    }
                }
                continue;
            }
            for (auto xcooryi : xcoory) {
                x1 = xcooryi.first;    // this window's small x index
                x2 = xrange[x1];    // this window's large x index
                y1 = xcooryi.second[0];    // this window's smallest y
                y2 = xcooryi.second[1];    // this window's largest y
                y1in = ymin;    // this window's smallest y index
                y2in = ymax;    // this window's largest y index
                for (indj = ymin; indj <= ymax; indj++) {
                    if (this->yn[indj] >= y2) {
                        y2in = indj;
                        break;
                    }
                    if (this->yn[indj] < y1) {
                        y1in = indj;
                    }
                }
                for (indj = y1in; indj <= y2in; indj++) {
                    for (indk = x1; indk <= x2; indk++) {
                        if (polyIn(this->xn[indk], this->yn[indj], indi)) {   // if this point is inside the polygon, set the z direction edge in markEdge with 1 and markNode with 1
                            for (indl = zi[this->conductorIn[indi].zmin]; indl < zi[this->conductorIn[indi].zmax]; indl++) {
                                this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indk * (this->N_cell_y + 1) + indj] = indi + 1;
                                compute_edgelink(indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indk * (this->N_cell_y + 1) + indj, node1, node2);
                                this->markNode[node1] = 1;
                                this->markNode[node2] = 1;
                            }
                        }
                        if (indk < x2) {
                            if (polyIn((this->xn[indk] + this->xn[indk + 1]) / 2, this->yn[indj], indi)) {   // if the edge along x axis is inside the polygon, set the x direction edge with 1 in markEdge
                                for (indl = zi[this->conductorIn[indi].zmin]; indl <= zi[this->conductorIn[indi].zmax]; indl++) {
                                    this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indk * (this->N_cell_y + 1) + indj] = indi + 1;
                                    compute_edgelink(indl * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indk * (this->N_cell_y + 1) + indj, node1, node2);
                                    this->markNode[node1] = 1;
                                    this->markNode[node2] = 1;
                                }
                            }
                        }
                        if (indj < y2in) {
                            if (polyIn(this->xn[indk], (this->yn[indj] + this->yn[indj + 1]) / 2, indi)) {    // if the edge along y axis is inside the polygon, set the y directio edge with 1 in markEdge
                                for (indl = zi[this->conductorIn[indi].zmin]; indl <= zi[this->conductorIn[indi].zmax]; indl++) {
                                    this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + indk * this->N_cell_y + indj] = indi + 1;
                                    compute_edgelink(indl * (this->N_edge_s + this->N_edge_v) + indk * this->N_cell_y + indj, node1, node2);
                                    this->markNode[node1] = 1;
                                    this->markNode[node2] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }


    }



    /* Find nodes inside conductors with linear complexity by recording the y inside each x range */
    void findInsideCond_xrangey(unordered_map<double, int> & xi, unordered_map<double, int> & yi, unordered_map<double, int> & zi){
        /* Find all the x ranges and put all the shown up y in the x ranges */
        int indi, indj, indl, indk;
        myint xrange_max;
        unordered_map<myint, myint> xrange;
        vector<myint> xcoorv;
        set<myint> xcoor;
        unordered_map<myint, vector<myint>> xcoory;    // store the start coordinate of the range, the end can be checked from xrange
        myint ss, ee;
        myint x1, x2;
        myint y1, y2;
        double y3;
        int mark1, mini_k, mark;
        double mini;



        for (indi = 0; indi < this->numCdtRow; indi++) {
            if (this->conductorIn[indi].zmax == this->conductorIn[indi].zmin)
                continue;
            xrange.clear();
            xcoor.clear();
            xcoorv.clear();
            xcoory.clear();
            mark1 = 0;
            //cout << endl;
            //cout << endl;
            for (indj = 0; indj < this->conductorIn[indi].numVert - 1; indj++) {
                //cout << this->conductorIn[indi].x[indj] << " " << this->conductorIn[indi].y[indj] << " ";
                if (this->conductorIn[indi].x[indj] != this->conductorIn[indi].x[indj + 1]) {   // line spanning across x

                    if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[indj + 1]) {
                        x1 = xi[this->conductorIn[indi].x[indj]];   // smaller x
                        x2 = xi[this->conductorIn[indi].x[indj + 1]];    // larger x
                    }
                    else {
                        x1 = xi[this->conductorIn[indi].x[indj + 1]];
                        x2 = xi[this->conductorIn[indi].x[indj]];
                    }

                    for (indl = x1; indl <= x2; indl++) {
                        xcoor.insert(indl);
                    }
                }
                else {   // line along y direction
                    xcoor.insert(xi[this->conductorIn[indi].x[indj]]);
                }
            }
            if (this->conductorIn[indi].x[indj] != this->conductorIn[indi].x[0]) {   // line spanning across x
                if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[0]) {
                    x1 = xi[this->conductorIn[indi].x[indj]];   // smaller x
                    x2 = xi[this->conductorIn[indi].x[0]];    // larger x
                }
                else {
                    x1 = xi[this->conductorIn[indi].x[0]];
                    x2 = xi[this->conductorIn[indi].x[indj]];
                }
                for (indl = x1; indl <= x2; indl++) {
                    xcoor.insert(indl);
                }
            }
            else {    // line along y direction
                xcoor.insert(xi[this->conductorIn[indi].x[indj]]);
            }

            for (auto xcoori : xcoor) {
                xcoorv.push_back(xcoori);
            }
            xrange_max = xcoorv.back();   // the maximal x coordinate

            for (indj = 0; indj < xcoorv.size() - 1; indj++) {
                mark1 = 1;    // the x coordinates are more than 1
                xrange[xcoorv[indj]] = xcoorv[indj + 1];
            }

            if (xcoorv.size() == 1) {    // If it has only one value
                xrange[xcoorv[0]] = xcoorv[0];
            }

            for (indj = 0; indj < this->conductorIn[indi].numVert - 1; indj++) {

                if (this->conductorIn[indi].y[indj] == this->conductorIn[indi].y[indj + 1]) {   // line along x direction
                    if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[indj + 1]) {
                        ss = xi[this->conductorIn[indi].x[indj]];
                        ee = xi[this->conductorIn[indi].x[indj + 1]];
                        if (ss == ee && mark1 == 1) {
                            continue;
                        }

                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is the last point
                            xcoory[ss].push_back(yi[this->conductorIn[indi].y[indj]]);
                            if (mark1 == 0) {   // xrange only has one value, break
                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                    else if (this->conductorIn[indi].x[indj] > this->conductorIn[indi].x[indj + 1]) {
                        ss = xi[this->conductorIn[indi].x[indj + 1]];
                        ee = xi[this->conductorIn[indi].x[indj]];
                        if (ss == ee && mark1 == 1) {
                            continue;
                        }
                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                            xcoory[ss].push_back(yi[this->conductorIn[indi].y[indj]]);
                            if (mark1 == 0) {    // xrange only has one value, break
                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                }
                else if (this->conductorIn[indi].y[indj] != this->conductorIn[indi].y[indj + 1] && this->conductorIn[indi].x[indj] != this->conductorIn[indi].x[indj + 1]) {   // this edge is with area slope
                    if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[indj + 1]) {
                        ss = xi[this->conductorIn[indi].x[indj]];
                        ee = xi[this->conductorIn[indi].x[indj + 1]];
                        if (ss == ee && mark1 == 1) {
                            continue;
                        }
                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                            y3 = (this->conductorIn[indi].x[indj + 1] - this->xn[ss]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[ss] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[indj + 1] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj + 1];
                            if (this->conductorIn[indi].y[indj] > this->conductorIn[indi].y[indj + 1]) {
                                y1 = yi[this->conductorIn[indi].y[indj + 1]];
                                y2 = yi[this->conductorIn[indi].y[indj]];
                            }
                            else {
                                y1 = yi[this->conductorIn[indi].y[indj]];
                                y2 = yi[this->conductorIn[indi].y[indj + 1]];
                            }
                            mini = DOUBLEMAX;
                            mini_k = -1;
                            for (indl = y1; indl <= y2; indl++) {    // find the closest y to y3
                                if (mini > abs(this->yn[indl] - y3)) {    // if there is a yn that is closer to y3
                                    mini = abs(this->yn[indl] - y3);
                                    mini_k = indl;
                                }
                            }
                            xcoory[ss].push_back(mini_k);
                            if (mark1 == 0) {   // if this range is the rightmost one, or xrange only has one value, break

                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                    else if (this->conductorIn[indi].x[indj] > this->conductorIn[indi].x[indj + 1]) {
                        ss = xi[this->conductorIn[indi].x[indj + 1]];
                        ee = xi[this->conductorIn[indi].x[indj]];
                        if (ss == ee && mark1 == 1) {
                            continue;
                        }
                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {
                            y3 = (this->conductorIn[indi].x[indj] - this->xn[ss]) / (this->conductorIn[indi].x[indj] - this->conductorIn[indi].x[indj + 1]) * this->conductorIn[indi].y[indj + 1] + (this->xn[ss] - this->conductorIn[indi].x[indj + 1]) / (this->conductorIn[indi].x[indj] - this->conductorIn[indi].x[indj + 1]) * this->conductorIn[indi].y[indj];
                            if (this->conductorIn[indi].y[indj] > this->conductorIn[indi].y[indj + 1]) {
                                y1 = yi[this->conductorIn[indi].y[indj + 1]];
                                y2 = yi[this->conductorIn[indi].y[indj]];
                            }
                            else {
                                y1 = yi[this->conductorIn[indi].y[indj]];
                                y2 = yi[this->conductorIn[indi].y[indj + 1]];
                            }
                            mini = DOUBLEMAX;
                            mini_k = -1;
                            for (indl = y1; indl <= y2; indl++) {    // find the closest y to y3
                                if (mini > abs(this->yn[indl] - y3)) {
                                    mini = abs(this->yn[indl] - y3);
                                    mini_k = indl;
                                }
                            }
                            xcoory[ss].push_back(mini_k);
                            if (mark1 == 0) {

                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                }
            }

            if (this->conductorIn[indi].y[indj] == this->conductorIn[indi].y[0]) {
                if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[0]) {
                    ss = xi[this->conductorIn[indi].x[indj]];
                    ee = xi[this->conductorIn[indi].x[0]];
                    if (!(ss == ee && mark1 == 1)) {

                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                            xcoory[ss].push_back(yi[this->conductorIn[indi].y[indj]]);
                            if (mark1 == 0) {

                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                }
                else if (this->conductorIn[indi].x[indj] > this->conductorIn[indi].x[0]) {
                    ss = xi[this->conductorIn[indi].x[0]];
                    ee = xi[this->conductorIn[indi].x[indj]];
                    if (!(ss == ee && mark1 == 1)) {

                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                            xcoory[ss].push_back(yi[this->conductorIn[indi].y[indj]]);
                            if (mark1 == 0) {

                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                }
            }
            else if (this->conductorIn[indi].y[indj] != this->conductorIn[indi].y[0] && this->conductorIn[indi].x[indj] != this->conductorIn[indi].x[0]) {   // this edge is with area slope
                if (this->conductorIn[indi].x[indj] < this->conductorIn[indi].x[0]) {
                    ss = xi[this->conductorIn[indi].x[indj]];
                    ee = xi[this->conductorIn[indi].x[0]];
                    if (!(ss == ee && mark1 == 1)) {

                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {   // if ss is not the last point, loop
                            y3 = (this->conductorIn[indi].x[0] - this->xn[ss]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[indj] + (this->xn[ss] - this->conductorIn[indi].x[indj]) / (this->conductorIn[indi].x[0] - this->conductorIn[indi].x[indj]) * this->conductorIn[indi].y[0];
                            if (this->conductorIn[indi].y[indj] > this->conductorIn[indi].y[0]) {
                                y1 = yi[this->conductorIn[indi].y[0]];
                                y2 = yi[this->conductorIn[indi].y[indj]];
                            }
                            else {
                                y1 = yi[this->conductorIn[indi].y[indj]];
                                y2 = yi[this->conductorIn[indi].y[0]];
                            }
                            mini = DOUBLEMAX;
                            mini_k = -1;
                            for (indl = y1; indl <= y2; indl++) {    // find the closest y to y3
                                if (mini > abs(this->yn[indl] - y3)) {
                                    mini = abs(this->yn[indl] - y3);
                                    mini_k = indl;
                                }
                            }
                            xcoory[ss].push_back(mini_k);
                            if (mark1 == 0) {   // if this range is the rightmost one, or xrange only has one value, break
                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                }
                else if (this->conductorIn[indi].x[indj] > this->conductorIn[indi].x[0]) {
                    ss = xi[this->conductorIn[indi].x[0]];
                    ee = xi[this->conductorIn[indi].x[indj]];
                    if (!(ss == ee && mark1 == 1)) {

                        while (xrange.find(ss) != xrange.end() && xrange[ss] <= ee) {    // if ss is not the last point, loop
                            y3 = (this->conductorIn[indi].x[indj] - this->xn[ss]) / (this->conductorIn[indi].x[indj] - this->conductorIn[indi].x[0]) * this->conductorIn[indi].y[0] + (this->xn[ss] - this->conductorIn[indi].x[0]) / (this->conductorIn[indi].x[indj] - this->conductorIn[indi].x[0]) * this->conductorIn[indi].y[indj];
                            if (this->conductorIn[indi].y[indj] > this->conductorIn[indi].y[0]) {
                                y1 = yi[this->conductorIn[indi].y[0]];
                                y2 = yi[this->conductorIn[indi].y[indj]];
                            }
                            else {
                                y1 = yi[this->conductorIn[indi].y[indj]];
                                y2 = yi[this->conductorIn[indi].y[0]];
                            }
                            mini = DOUBLEMAX;
                            mini_k = -1;
                            for (indl = y1; indl <= y2; indl++) {    // find the closest y to y3
                                if (mini > abs(this->yn[indl] - y3)) {
                                    mini = abs(this->yn[indl] - y3);
                                    mini_k = indl;
                                }
                            }
                            xcoory[ss].push_back(mini_k);
                            if (mark1 == 0) {

                                break;
                            }
                            ss = xrange[ss];
                        }
                    }
                }

            }


            if (mark1 == 0) {    // only has one x
                for (auto xcooryi : xcoory) {    // only one xcooryi
                    vector<myint> xcooryiv = xcooryi.second;
                    sort(xcooryiv.begin(), xcooryiv.end());
                    mark = 0;
                    for (auto xrangey : xcooryiv) {
                        mark++;
                        if (mark % 2 == 1) {
                            y1 = xrangey;
                            continue;
                        }
                        if (mark % 2 == 0) {
                            y2 = xrangey;
                        }

                        indj = xcooryi.first;
                        for (indl = zi[this->conductorIn[indi].zmin]; indl <= zi[this->conductorIn[indi].zmax]; indl++) {
                            for (indk = y1; indk < y2; indk++) {
                                this->markNode[indl * this->N_node_s + indj * (this->N_cell_y + 1) + indk] = (myint)1;
                                this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + indj * (this->N_cell_y) + indk] = indi + (myint)1;

                                if (indl != zi[this->conductorIn[indi].zmax]) {
                                    this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                                }
                            }
                            this->markNode[indl * this->N_node_s + indj * (this->N_cell_y + 1) + indk] = (myint)1;

                            if (indl != zi[this->conductorIn[indi].zmax]) {
                                this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                            }
                        }

                    }
                }
                continue;
            }
            for (auto xcooryi : xcoory) {
                vector<myint> xcooryiv = xcooryi.second;
                sort(xcooryiv.begin(), xcooryiv.end());
                mark = 0;
                x1 = xcooryi.first;
                x2 = xrange[xcooryi.first];
                //cout << endl << endl;
                //cout << "x are " << this->xn[x1] << " " << this->xn[x2] << " : " << endl;
                mark = 0;
                for (auto xrangey : xcooryiv) {
                    mark++;
                    if (mark % 2 == 1) {
                        y1 = xrangey;
                        continue;
                    }
                    else if (mark % 2 == 0) {
                        y2 = xrangey;
                        //cout << this->yn[y1] << " " << this->yn[y2] << " ";
                        //cout << this->conductorIn[indi].zmax << " " << this->conductorIn[indi].zmin;
                        for (indj = x1; indj <= x2; indj++) {
                            for (indl = zi[this->conductorIn[indi].zmin]; indl <= zi[this->conductorIn[indi].zmax]; indl++) {
                                for (indk = y1; indk < y2; indk++) {
                                    if (indj < x2){
                                        this->markNode[indl * this->N_node_s + indj * (this->N_cell_y + 1) + indk] = (myint)1;
                                        this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + indj * (this->N_cell_y) + indk] = indi + (myint)1;
                                        this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                                        if (indl != zi[this->conductorIn[indi].zmax]) {
                                            this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                                        }
                                    }
                                    else {    // If this range is the rightmost range, include the rightmost point
                                        this->markNode[indl * this->N_node_s + indj * (this->N_cell_y + 1) + indk] = (myint)1;
                                        this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + indj * (this->N_cell_y) + indk] = indi + (myint)1;
                                        if (indl != zi[this->conductorIn[indi].zmax]) {
                                            this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                                        }
                                    }
                                }
                                if (indj < x2){
                                    this->markNode[indl * this->N_node_s + indj * (this->N_cell_y + 1) + indk] = (myint)1;
                                    this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                                    if (indl != zi[this->conductorIn[indi].zmax]) {
                                        this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                                    }
                                }
                                else {    // If this range is the rightmost range, include the rightmost point
                                    this->markNode[indl * this->N_node_s + indj * (this->N_cell_y + 1) + indk] = (myint)1;
                                    if (indl != zi[this->conductorIn[indi].zmax]) {
                                        this->markEdge[indl * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indj * (this->N_cell_y + 1) + indk] = indi + (myint)1;
                                    }
                                }
                            }

                        }
                    }
                }
            }

        }

    }

    /* Try to resolve the issue when the port node is outside the conductor */
    int findPortNode(int step, int x1n, int y1n, int z1n, int x2n, int y2n, int z2n, myint indi, myint indk) {
        /* return 1 means find the port point
        return 0 means not find the port point
        indi : the port #
        indk : this port's multiplicity*/
        if (this->markNode[z1n * this->N_node_s + x1n * (this->N_cell_y + 1) + y1n] == 0) {   // the first port node is not inside the conductor
            int i = 0;
            if (x1n < x2n) {   // x1n goes further left
                while (i < step && x1n > 0) {
                    x1n -= 1;
                    if (this->markNode[z1n * this->N_node_s + x1n * (this->N_cell_y + 1) + y1n]) {
                        this->portCoor[indi].x1[indk] = this->xn[x1n];
                        break;
                    }
                    i++;
                }
                if (i == step || x1n == 0) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (y1n < y2n) {   // y1n goes further front
                while (i < step && y1n > 0) {
                    y1n -= 1;
                    if (this->markNode[z1n * this->N_node_s + x1n * (this->N_cell_y + 1) + y1n]) {
                        this->portCoor[indi].y1[indk] = this->yn[y1n];
                        break;
                    }
                    i++;
                }
                if (i == step || y1n == 0) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (z1n < z2n) {   // z1n goes further down
                while (i < step && z1n > 0) {
                    z1n -= 1;
                    if (this->markNode[z1n * this->N_node_s + x1n * (this->N_cell_y + 1) + y1n]) {
                        this->portCoor[indi].z1[indk] = this->zn[z1n];
                        break;
                    }
                    i++;
                }
                if (i == step || z1n == 0) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (x1n > x2n) {   // x1n goes further right
                while (i < step && x1n < this->nx - 1) {
                    x1n += 1;
                    if (this->markNode[z1n * this->N_node_s + x1n * (this->N_cell_y + 1) + y1n]) {
                        this->portCoor[indi].x1[indk] = this->xn[x1n];
                        break;
                    }
                    i++;
                }
                if (i == step || x1n == this->nx - 1) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (y1n > y2n) {   // y1n goes further back
                while (i < step && y1n < this->ny - 1) {
                    y1n += 1;
                    if (this->markNode[z1n * this->N_node_s + x1n * (this->N_cell_y + 1) + y1n]) {
                        this->portCoor[indi].y1[indk] = this->yn[y1n];
                        break;
                    }
                    i++;
                }
                if (i == step || y1n == this->ny - 1) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (z1n > z2n) {   // z1n goes further up
                while (i < step && z1n < this->nz - 1) {
                    z1n += 1;
                    if (this->markNode[z1n * this->N_node_s + x1n * (this->N_cell_y + 1) + y1n]) {
                        this->portCoor[indi].z1[indk] = this->zn[z1n];
                        break;
                    }
                    i++;
                }
                if (i == step || z1n == this->nz - 1) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
        }
        if (this->markNode[z2n * this->N_node_s + x2n * (this->N_cell_y + 1) + y2n] == 0) {   // the second port node is not inside the conductor
            int i = 0;
            if (x1n < x2n) {   // x2n goes further right
                while (i < step && x2n < this->nx - 1) {
                    x2n += 1;
                    if (this->markNode[z2n * this->N_node_s + x2n * (this->N_cell_y + 1) + y2n]) {
                        this->portCoor[indi].x2[indk] = this->xn[x2n];
                        break;
                    }
                    i++;
                }
                if (i == step || x2n == this->nx - 1) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (y1n < y2n) {   // y2n goes further back
                while (i < step && y2n < this->ny - 1) {
                    y2n += 1;
                    if (this->markNode[z2n * this->N_node_s + x2n * (this->N_cell_y + 1) + y2n]) {
                        this->portCoor[indi].y2[indk] = this->yn[y2n];
                        break;
                    }
                    i++;
                }
                if (i == step || y2n == this->ny - 1) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (z1n < z2n) {   // z2n goes further up
                while (i < step && z2n < this->nz - 1) {
                    z2n += 1;
                    if (this->markNode[z2n * this->N_node_s + x2n * (this->N_cell_y + 1) + y2n]) {
                        this->portCoor[indi].z2[indk] = this->zn[z2n];
                        break;
                    }
                    i++;
                }
                if (i == step || z2n == this->nz - 1) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (x1n > x2n) {   // x2n goes further left
                while (i < step && x2n > 0) {
                    x2n -= 1;
                    if (this->markNode[z2n * this->N_node_s + x2n * (this->N_cell_y + 1) + y2n]) {
                        this->portCoor[indi].x2[indk] = this->xn[x2n];
                        break;
                    }
                    i++;
                }
                if (i == step || x2n == 0) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (y1n > y2n) {   // y2n goes further front
                while (i < step && y2n > 0) {
                    y2n -= 1;
                    if (this->markNode[z2n * this->N_node_s + x2n * (this->N_cell_y + 1) + y2n]) {
                        this->portCoor[indi].y2[indk] = this->yn[y2n];
                        break;
                    }
                    i++;
                }
                if (i == step || y2n == 0) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
            else if (z1n > z2n) {   // z1n goes further down
                while (i < step && z2n > 0) {
                    z2n -= 1;
                    if (this->markNode[z2n * this->N_node_s + x2n * (this->N_cell_y + 1) + y2n]) {
                        this->portCoor[indi].z2[indk] = this->zn[z2n];
                        break;
                    }
                    i++;
                }
                if (i == step || z2n == 0) {   // doesn't find a point that is inside the conductor
                    return 0;
                }
            }
        }
        return 1;
    }

    /* Find this conductor node's corresponding conductorIn */
    void findCond2CondIn(int inx, int iny, int inz){
        /* inx : this conductor node's x index
        iny : this conductor node's y index
        inz : this conductor node's z index */
        myint eno;

        if (inz != 0) {    // this node is not on the bottom plane
            eno = (inz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + inx * (this->N_cell_y + 1) + iny;    // the node's lower edge
            if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
                this->cond2condIn.insert(this->markEdge[eno]);
            }
        }
        if (inz != this->nz - 1) {    // this node is not on the upper plane
            eno = inz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + inx * (this->N_cell_y + 1) + iny;    // the node's lower edge
            if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
                this->cond2condIn.insert(this->markEdge[eno]);
            }
        }
        if (inx != 0) {    // this node is not on the left plane
            eno = inz *(this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (inx - 1) * (this->N_cell_y + 1) + iny;
            if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
                this->cond2condIn.insert(this->markEdge[eno]);
            }
        }
        if (inx != this->nx - 1) {    // this node is not on the right plane
            eno = inz *(this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + inx * (this->N_cell_y + 1) + iny;
            if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
                this->cond2condIn.insert(this->markEdge[eno]);
            }
        }
        if (iny != 0) {    // this node is not on the front plane
            eno = inz *(this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny - 1;
            if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
                this->cond2condIn.insert(this->markEdge[eno]);
            }
        }
        if (iny != this->ny - 1) {    // this node is not on the back plane
            eno = inz *(this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny;
            if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
                this->cond2condIn.insert(this->markEdge[eno]);
            }
        }
    }

    /* Find boundary nodes and edges */
    void findBoundNodeEdge(int inx, int iny, int inz) {
        /* inx : this conductor node's x index
        iny : this conductor node's y index
        inz : this conductor node's z index */
        myint eno;

        if (inz == 0) {
            this->lbdn.insert(inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny);   // this node is lower PEC boundary

            if (inx != 0) {    // this node is not on the left plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (inx - 1) * (this->N_cell_y + 1) + iny;   // left edge
                if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
                    this->lbde.insert(eno);
                }
            }
            if (inx != this->nx - 1) {    // this node is not on the right plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + inx * (this->N_cell_y + 1) + iny;   // right edge
                if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
                    this->lbde.insert(eno);
                }
            }
            if (iny != 0) {    // this node is not on the front plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny - 1;   // front edge
                if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
                    this->lbde.insert(eno);
                }
            }
            if (iny != this->ny - 1) {    // this node is not on the back plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny;   // back edge
                if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
                    this->lbde.insert(eno);
                }
            }
        }
        else if (inz == this->nz - 1) {
            this->ubdn.insert(inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny);   // this node is upper PEC boundary
            if (inx != 0) {    // this node is not on the left plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (inx - 1) * (this->N_cell_y + 1) + iny;    // left edge
                if (this->markEdge[eno] != 0 && this->ubde.find(eno) == this->ubde.end()) {
                    this->ubde.insert(eno);
                }
            }
            if (inx != this->nx - 1) {    // this node is not on the right plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + inx * (this->N_cell_y + 1) + iny;    // right edge
                if (this->markEdge[eno] != 0 && this->ubde.find(eno) == this->ubde.end()) {
                    this->ubde.insert(eno);
                }
            }
            if (iny != 0) {    // this node is not on the front plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny - 1;    // front edge
                if (this->markEdge[eno] != 0 && this->ubde.find(eno) == this->ubde.end()) {
                    this->ubde.insert(eno);
                }
            }
            if (iny != this->ny - 1) {    // this node is not on the back plane
                eno = inz *(this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny;    // back edge
                if (this->markEdge[eno] != 0 && this->ubde.find(eno) == this->ubde.end()) {
                    this->ubde.insert(eno);
                }
            }
        }
    }

    /* Map the original edges to the new edge # with PEC boundary */
    void setMapEdge() {
        /* if this edge is removed, value = -1
        else value = new edge # */
        this->mapEdge = new myint[this->N_edge]();
        this->mapEdgeR = new myint[this->N_edge - this->bden]();
        int ind, count = 0;    // count : how many boundary edges already identified

        for (ind = 0; ind < this->N_edge; ind++) {
            if (ind >= 0 && ind < this->N_edge_s && this->lbde.find(ind) != this->lbde.end()) {   // if this edge is on the lower PEC boundary
                this->mapEdge[ind] = -1;
                count++;
            }
            else if (ind >= this->N_edge - this->N_edge_s && ind < this->N_edge && this->ubde.find(ind) != this->ubde.end()) {   // if this edge is on the upper PEC boundary
                this->mapEdge[ind] = -1;
                count++;
            }
            else {    // if this edge is not on the PEC boundaries
                this->mapEdge[ind] = ind - count;
                this->mapEdgeR[ind - count] = ind;
            }
        }
    }

    /* Generate V0 and V0a */
    //void merge_V0(myint &v0num, myint &leng_v0, myint &v0anum, myint &leng_v0a, myint* map) {
    //    v0num = 0;
    //    leng_v0 = 0;
    //    v0anum = 0;
    //    leng_v0a = 0;
    //    myint count = 1, nno, eno;
    //    int* visited;
    //    int indx, indy;

    //    for (int iz = 0; iz < this->nz; iz++) {    // generate on each layer
    //        for (int ix = 0; ix < this->nx; ix++) {
    //            for (int iy = 0; iy < this->ny; iy++) {
    //                nno = iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy;
    //                if ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(nno) == this->lbdn.end() && this->ubdn.find(nno) == this->ubdn.end())) {   // if this node is not boundary node
    //                    if (visited[ix * (this->N_cell_y + 1) + iy] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] == 0) {

    //                        map[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;

    //                        if (iz != 0) {    // this node is not on the bottom plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (iz != this->nz - 1) {   // this node is not on the top plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indx != 0) {    // this node is not on the left plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indx != this->nx - 1) {    // this node is not on the right plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indy != 0) {    // this node is not on the front plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indy != this->ny - 1) {   // this node is not on the back plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                    }

    //                    count++;
    //                }
    //            }
    //        }
    //    }
    //    for (auto node : this->ubdn) {
    //        iz = sys->nz - 1;
    //        indx = (node % this->N_node_s) / (this->N_cell_y + 1);
    //        indy = (node % this->N_node_s) % (this->N_cell_y + 1);
    //        map[node] = count;
    //        if (iz != 0) {    // this node is not on the bottom plane
    //            v0d1num++;
    //            v0d1anum++;
    //        }
    //        //if (iz != this->nz - 1) {   // this node is not on the top plane
    //        //    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
    //        //    compute_edgelink(eno, node1, node2);
    //        //    if (node1 != ndi) {
    //        //        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
    //        //            v0d1num++;
    //        //            v0d1anum++;
    //        //        }
    //        //    }
    //        //    else if (node2 != ndi) {
    //        //        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
    //        //            v0d1num++;
    //        //            v0d1anum++;
    //        //        }
    //        //    }
    //        //}
    //        if (indx != 0) {    // this node is not on the left plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
    //            compute_edgelink(eno, node1, node2);
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }

    //        }
    //        if (indx != this->nx - 1) {    // this node is not on the right plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //        }
    //        if (indy != 0) {    // this node is not on the front plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
    //            compute_edgelink(eno, node1, node2);
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //        }
    //        if (indy != this->ny - 1) {   // this node is not on the back plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
    //            compute_edgelink(eno, node1, node2);
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //        }
    //    }
    //    

    //    this->V0RowId = (myint*)malloc(v0num, sizeof(myint));
    //    this->v0ColId = (myint*)malloc(v0num, sizeof(myint));
    //    this->v0val = (double*)malloc(v0num, sizeof(myint));
    //    this->v0aval = (double*)malloc(v0anum, sizeof(double));
    //    double lx_whole_avg = 0;
    //    double ly_whole_avg = 0;
    //    double lz_whole_avg = 0;
    //    lx_whole_avg = (this->xn[this->nx - 1] - this->xn[0]) / (this->nx - 1);
    //    ly_whole_avg = (this->yn[this->ny - 1] - this->yn[0]) / (this->ny - 1);
    //    lz_whole_avg = (this->zn[this->nz - 1] - this->zn[0]) / (this->nz - 1);
    //    v0num = 0;
    //    v0anum = 0;

    //    for (int iz = 0; iz < this->nz; iz++) {    // generate on each layer
    //        for (int ix = 0; ix < this->nx; ix++) {
    //            for (int iy = 0; iy < this->ny; iy++) {
    //                nno = iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy;
    //                if ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(nno) == this->lbdn.end() && this->ubdn.find(nno) == this->ubdn.end())) {   // if this node is not boundary node
    //                    if (visited[ix * (this->N_cell_y + 1) + iy] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] == 0) {


    //                        avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
    //                        if (iz != 0) {    // this node is not on the bottom plane
    //                            eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
    //                            this->v0RowId[v0num] = eno;
    //                            this->v0ColId[v0num] = leng_v0;
    //                            this->V0val[v0num] = -1 / (this->zn[iz] - this->zn[iz - 1]);
    //                            this->v0aval[v0anum] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);;
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (iz != this->nz - 1) {   // this node is not on the top plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indx != 0) {    // this node is not on the left plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indx != this->nx - 1) {    // this node is not on the right plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indy != 0) {    // this node is not on the front plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                        if (indy != this->ny - 1) {   // this node is not on the back plane
    //                            v0num++;
    //                            v0anum++;
    //                        }
    //                    }

    //                    count++;
    //                }
    //            }
    //        }
    //    }
    //    for (auto node : this->ubdn) {
    //        iz = sys->nz - 1;
    //        indx = (node % this->N_node_s) / (this->N_cell_y + 1);
    //        indy = (node % this->N_node_s) % (this->N_cell_y + 1);
    //        if (iz != 0) {    // this node is not on the bottom plane
    //            v0d1num++;
    //            v0d1anum++;
    //        }
    //        //if (iz != this->nz - 1) {   // this node is not on the top plane
    //        //    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
    //        //    compute_edgelink(eno, node1, node2);
    //        //    if (node1 != ndi) {
    //        //        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
    //        //            v0d1num++;
    //        //            v0d1anum++;
    //        //        }
    //        //    }
    //        //    else if (node2 != ndi) {
    //        //        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
    //        //            v0d1num++;
    //        //            v0d1anum++;
    //        //        }
    //        //    }
    //        //}
    //        if (indx != 0) {    // this node is not on the left plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
    //            compute_edgelink(eno, node1, node2);
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }

    //        }
    //        if (indx != this->nx - 1) {    // this node is not on the right plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //        }
    //        if (indy != 0) {    // this node is not on the front plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
    //            compute_edgelink(eno, node1, node2);
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //        }
    //        if (indy != this->ny - 1) {   // this node is not on the back plane
    //            eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
    //            compute_edgelink(eno, node1, node2);
    //            if (node1 != node) {
    //                if (this->ubdn.find(node1) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //            else if (node2 != node) {
    //                if (this->ubdn.find(node2) == this->ubdn.end()) {
    //                    v0num++;
    //                    v0anum++;
    //                }
    //            }
    //        }
    //    }
    //}

    /* Generate V0d: both V0d1 and V0d2 are put into V0d1 */
    void merge_v0d1(double block1_x, double block1_y, double block2_x, double block2_y, double block3_x, double block3_y, myint &v0d1num, myint &leng_v0d1, myint &v0d1anum, myint &leng_v0d1a, myint *map, double sideLen){
        int *visited;
        clock_t t1;
        double t, ta;
        double ratio;
        double startx, starty;    // the start coordinates of each block
        queue<int> st;    // dfs stack
        //vector<int> ind;
        int indsize;
        myint indi = 0;
        int indx, indy;
        int mark;
        int indnum;
        int *markLayerNode = (int*)calloc(this->N_node_s, sizeof(int));
        /* Mark layer nodes from port sides */
        for (int indPort = 0; indPort < this->numPorts; indPort++) {

            for (int indPortSide = 0; indPortSide < this->portCoor[indPort].multiplicity; indPortSide++) {
                myint indCdt = this->portCoor[indPort].portCnd[indPortSide] - 1; // Conductor index for this port side

                for (int indCdtNode = 0; indCdtNode < this->cdtNumNode[indCdt]; indCdtNode++) {
                    markLayerNode[this->conductor[indCdt].node[indCdtNode] % (this->N_node_s)] = 1;
                }
            }
        }

        leng_v0d1 = 0;
        leng_v0d1a = 0;
        v0d1num = 0;
        v0d1anum = 0;

        /* First assign a larger number of storage, don't need to calculate the entries twice */
        //myint *v0d1RowId = (myint*)malloc(2 * this->outedge * sizeof(myint));
        //myint *v0d1ColId = (myint*)malloc(2 * this->outedge * sizeof(myint));
        //double *v0d1val = (double*)malloc(2 * this->outedge * sizeof(double));
        ///*myint *v0d1aRowId = (myint*)malloc(this->N_edge * sizeof(myint));
        //myint *v0d1aColId = (myint*)malloc(this->N_edge * sizeof(myint));*/
        //double *v0d1aval = (double*)malloc(2 * this->outedge * sizeof(double));

        /* V0d1 generation */
        int count = 1;    /* count which box it is */
        clock_t t2 = clock();
        unordered_map<myint, double> va, v;
        vector<set<myint>> node_group;
        set<myint> base;
        myint eno, nno;
        double lx_avg, ly_avg, lz_avg;
        t = 0.;
        ta = 0.;
        myint node1, node2;
        int nodegs;   // node group #
        for (int iz = 0; iz < this->nz; iz++) {    // merge on each layer
            visited = (int*)calloc(this->nx * this->ny, sizeof(int));
            for (int ix = 0; ix < this->nx; ix++) {
                for (int iy = 0; iy < this->ny; iy++) {
                    nno = iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy;
                    if ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(nno) == this->lbdn.end() && this->ubdn.find(nno) == this->ubdn.end())) {   // if this node is not boundary node
                        if (visited[ix * (this->N_cell_y + 1) + iy] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] == 0) {
                            if (markLayerNode[ix * (this->N_cell_y + 1) + iy] == 0 && !this->markProSide[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this point is not visited and it is outside the conductor, not in the projection of the excited conductor
                                //if (!ind.empty())
                                //    ind.clear();
                                startx = this->xn[ix];
                                starty = this->yn[iy];
                                node_group.push_back(base);
                                nodegs = node_group.size() - 1;
                                node_group[nodegs].insert(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);
                                map[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
                                //ind.push_back(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);
                                st.push(ix * (this->N_cell_y + 1) + iy);
                                visited[ix * (this->N_cell_y + 1) + iy] = 1;

                                while (!st.empty()) {
                                    indx = (st.front()) / (this->N_cell_y + 1);
                                    indy = st.front() % (this->N_cell_y + 1);
                                    if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
                                        if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block1_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block1_y) {   // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() + this->N_cell_y + 1] == 0
                                                && visited[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
                                                && markLayerNode[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
                                                && !this->markProSide[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited, not in the markLayerNode (Projection of the excited conductors), not in the projection side, not among the boundary nodes
                                                st.push((indx + 1)*(this->N_cell_y + 1) + indy);
                                                visited[(indx + 1)*(this->N_cell_y + 1) + indy] = 1;
                                                map[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);
                                            }
                                        }
                                    }
                                    if (indx != 0) {    // it must have a left x edge, thus left x node
                                        if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block1_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block1_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() - this->N_cell_y - 1] == 0
                                                && visited[(indx - 1) * (this->N_cell_y + 1) + indy] == 0
                                                && markLayerNode[(indx - 1) * (this->N_cell_y + 1) + indy] == 0 && !this->markProSide[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited, not in the markLayerNode (Projection of the excited conductors), not in the projection side, not among the boundary nodes

                                                st.push((indx - 1)*(this->N_cell_y + 1) + indy);
                                                visited[(indx - 1)*(this->N_cell_y + 1) + indy] = 1;
                                                map[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
                                            }
                                        }
                                    }
                                    if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
                                        if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block1_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block1_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() + 1] == 0
                                                && visited[indx * (this->N_cell_y + 1) + indy + 1] == 0
                                                && markLayerNode[indx * (this->N_cell_y + 1) + indy + 1] == 0
                                                && !this->markProSide[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx)*(this->N_cell_y + 1) + indy + 1);
                                                visited[(indx)*(this->N_cell_y + 1) + indy + 1] = 1;
                                                map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1);
                                            }
                                        }
                                    }
                                    if (indy != 0) {    // it must have a closer y edge, thus closer y node
                                        if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block1_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block1_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
                                                && visited[(indx)* (this->N_cell_y + 1) + indy - 1] == 0
                                                && markLayerNode[(indx)* (this->N_cell_y + 1) + indy - 1] == 0
                                                && !this->markProSide[iz * this->N_node_s + (indx)* (this->N_cell_y + 1) + indy - 1]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx)*(this->N_cell_y + 1) + indy - 1);
                                                visited[(indx)*(this->N_cell_y + 1) + indy - 1] = 1;
                                                map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1);
                                            }
                                        }
                                    }
                                    st.pop();
                                }

                                for (auto ndi : node_group[nodegs]) {
                                    indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
                                    indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
                                    avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
                                    if (iz != 0) {    // this node is not on the bottom plane
                                        eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (iz != this->nz - 1) {   // this node is not on the top plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indx != 0) {    // this node is not on the left plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }

                                    }
                                    if (indx != this->nx - 1) {    // this node is not on the right plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indy != 0) {    // this node is not on the front plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indy != this->ny - 1) {   // this node is not on the back plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                }

                                count++;
                            }

                            else if (markLayerNode[ix * (this->N_cell_y + 1) + iy] == 1 && !this->markProSide[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {//&& this->exciteCdtLayer[iz] == 1) {    // this point is not visited and it is outside the conductor, in the projection of the excited conductor
                                startx = this->xn[ix];
                                starty = this->yn[iy];

                                map[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
                                st.push(ix * (this->N_cell_y + 1) + iy);
                                visited[ix * (this->N_cell_y + 1) + iy] = 1;
                                node_group.push_back(base);
                                nodegs = node_group.size() - 1;
                                node_group[nodegs].insert(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);
                                while (!st.empty()) {
                                    mark = 0;
                                    indx = (st.front()) / (this->N_cell_y + 1);
                                    indy = st.front() % (this->N_cell_y + 1);

                                    if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
                                        if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() + this->N_cell_y + 1] == 0
                                                && visited[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
                                                && markLayerNode[(indx + 1) * (this->N_cell_y + 1) + indy] == 1
                                                && !this->markProSide[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited and not the boundary nodes

                                                st.push((indx + 1)*(this->N_cell_y + 1) + indy);
                                                visited[(indx + 1)*(this->N_cell_y + 1) + indy] = 1;
                                                map[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);
                                            }
                                        }
                                    }

                                    if (indx != 0) {    // it must have a left x edge, thus left x node
                                        if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() - this->N_cell_y - 1] == 0
                                                && visited[(indx - 1) * (this->N_cell_y + 1) + indy] == 0
                                                && markLayerNode[(indx - 1) * (this->N_cell_y + 1) + indy] == 1
                                                && !this->markProSide[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx - 1)*(this->N_cell_y + 1) + indy);
                                                visited[(indx - 1)*(this->N_cell_y + 1) + indy] = 1;
                                                map[iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy);
                                            }
                                        }
                                    }
                                    if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
                                        if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block2_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() + 1] == 0
                                                && visited[indx * (this->N_cell_y + 1) + indy + 1] == 0
                                                && markLayerNode[indx * (this->N_cell_y + 1) + indy + 1] == 1
                                                && !this->markProSide[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx)*(this->N_cell_y + 1) + indy + 1);
                                                visited[(indx)*(this->N_cell_y + 1) + indy + 1] = 1;
                                                map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1);
                                            }
                                        }
                                    }
                                    if (indy != 0) {    // it must have a closer y edge, thus closer y node
                                        if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block2_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
                                                && visited[(indx)* (this->N_cell_y + 1) + indy - 1] == 0
                                                && markLayerNode[(indx)* (this->N_cell_y + 1) + indy - 1] == 1
                                                && !this->markProSide[iz * this->N_node_s + (indx)* (this->N_cell_y + 1) + indy - 1]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx)*(this->N_cell_y + 1) + indy - 1);
                                                visited[(indx)*(this->N_cell_y + 1) + indy - 1] = 1;
                                                map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1);
                                            }
                                        }
                                    }
                                    st.pop();
                                }

                                for (auto ndi : node_group[nodegs]) {
                                    indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
                                    indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
                                    avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
                                    if (iz != 0) {    // this node is not on the bottom plane
                                        eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (iz != this->nz - 1) {   // this node is not on the top plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indx != 0) {    // this node is not on the left plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indx != this->nx - 1) {    // this node is not on the right plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indy != 0) {    // this node is not on the front plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indy != this->ny - 1) {   // this node is not on the back plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                }


                                count++;

                            }

                            else {
                                startx = this->xn[ix];
                                starty = this->yn[iy];

                                map[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
                                st.push(ix * (this->N_cell_y + 1) + iy);
                                visited[ix * (this->N_cell_y + 1) + iy] = 1;
                                node_group.push_back(base);
                                nodegs = node_group.size() - 1;
                                node_group[nodegs].insert(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);

                                while (!st.empty()) {
                                    mark = 0;
                                    indx = (st.front()) / (this->N_cell_y + 1);
                                    indy = st.front() % (this->N_cell_y + 1);

                                    if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
                                        if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block3_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block3_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() + this->N_cell_y + 1] == 0
                                                && visited[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
                                                && this->markProSide[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in the sideLen && this node is not among the boundary nodes

                                                st.push((indx + 1)*(this->N_cell_y + 1) + indy);
                                                visited[(indx + 1)*(this->N_cell_y + 1) + indy] = 1;
                                                map[iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy);
                                            }
                                        }
                                    }
                                    if (indx != 0) {    // it must have a left x edge, thus left x node
                                        if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block3_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block3_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() - this->N_cell_y - 1] == 0
                                                && visited[(indx - 1) * (this->N_cell_y + 1) + indy] == 0
                                                && this->markProSide[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx - 1)*(this->N_cell_y + 1) + indy);
                                                visited[(indx - 1)*(this->N_cell_y + 1) + indy] = 1;
                                                map[iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy);
                                            }
                                        }
                                    }
                                    if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
                                        if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block3_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block3_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() + 1] == 0
                                                && visited[indx * (this->N_cell_y + 1) + indy + 1] == 0
                                                && this->markProSide[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx)*(this->N_cell_y + 1) + indy + 1);
                                                visited[(indx)*(this->N_cell_y + 1) + indy + 1] = 1;
                                                map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1);
                                            }
                                        }
                                    }
                                    if (indy != 0) {    // it must have a closer y edge, thus closer y node
                                        if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block3_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block3_y) {    // this node is within the block area
                                            if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
                                                && visited[(indx)* (this->N_cell_y + 1) + indy - 1] == 0
                                                && this->markProSide[iz * this->N_node_s + (indx)* (this->N_cell_y + 1) + indy - 1]
                                                && ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

                                                st.push((indx)*(this->N_cell_y + 1) + indy - 1);
                                                visited[(indx)*(this->N_cell_y + 1) + indy - 1] = 1;
                                                map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1] = count;
                                                node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1);
                                            }
                                        }
                                    }
                                    st.pop();
                                }
                                for (auto ndi : node_group[nodegs]) {
                                    indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
                                    indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
                                    avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
                                    if (iz != 0) {    // this node is not on the bottom plane
                                        eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (iz != this->nz - 1) {   // this node is not on the top plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indx != 0) {    // this node is not on the left plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indx != this->nx - 1) {    // this node is not on the right plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indy != 0) {    // this node is not on the front plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                    if (indy != this->ny - 1) {   // this node is not on the back plane
                                        eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
                                        compute_edgelink(eno, node1, node2);
                                        if (node1 != ndi) {
                                            if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                        else if (node2 != ndi) {
                                            if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                                                v0d1num++;
                                                v0d1anum++;
                                            }
                                        }
                                    }
                                }

                                count++;
                            }
                        }
                    }
                }
            }
            free(visited); visited = NULL;
        }

        /* V0d2 generation */
        myint indj;
        myint inz, inx, iny;
        myint iz, ix, iy;
        queue<myint> qu;
        visited = (int*)calloc(this->N_node, sizeof(int));
        t2 = clock();


        for (indi = 0; indi < this->numCdt; indi++) {
            //cout << this->conductor[indi].markPort << " ";
            if (this->conductor[indi].markPort == -1) {   // if this conductor is the reference conductor, no V0d2 corresponding to it
                continue;
            }
            else {
                mark = 0;    // if mark = 0 it means that no V0d2 for this conductor, leng_v0d doesn't increase by 1
                //v.clear();
                //va.clear();
                for (indj = 0; indj < this->cdtNumNode[indi]; indj++) {
                    map[this->conductor[indi].node[indj]] = count;
                    iz = this->conductor[indi].node[indj] / this->N_node_s;
                    ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
                    iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
                    avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
                    if (iz != 0) {    // this node is not on the bottom plane
                        eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
                        if (this->markEdge[eno] == 0 && this->markNode[(iz - 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this edge is in the dielectric
                            v0d1num++;
                            v0d1anum++;
                            mark = 1;
                        }
                    }
                    if (iz != this->nz - 1) {   // this node is not on the top plane
                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
                        if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                            v0d1num++;
                            v0d1anum++;
                            mark = 1;
                        }
                    }
                    if (ix != 0) {    // this node is not on the left plane
                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
                        if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                            v0d1num++;
                            v0d1anum++;
                            mark = 1;
                        }
                    }
                    if (ix != this->nx - 1) {    // this node is not on the right plane
                        eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
                        if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                            v0d1num++;
                            v0d1anum++;
                            mark = 1;
                        }
                    }
                    if (iy != 0) {    // this node is not on the front plane
                        eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
                        if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                            v0d1num++;
                            v0d1anum++;
                            mark = 1;
                        }
                    }
                    if (iy != this->ny - 1) {   // this node is not on the back plane
                        eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
                        if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                            v0d1num++;
                            v0d1anum++;
                            mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
                        }
                    }
                }
                if (mark == 1) {
                    count++;
                }

            }
        }

        /* Sparse matrix construction for V0d1 */
        this->v0d1RowId = (myint*)malloc(v0d1num * sizeof(myint));
        this->v0d1ColId = (myint*)malloc(v0d1num * sizeof(myint));
        this->v0d1val = (double*)malloc(v0d1num * sizeof(double));
        this->v0d1aval = (double*)malloc(v0d1anum * sizeof(double));

        double lx_whole_avg = 0;
        double ly_whole_avg = 0;
        double lz_whole_avg = 0;
        lx_whole_avg = (this->xn[this->nx - 1] - this->xn[0]) / (this->nx - 1);
        ly_whole_avg = (this->yn[this->ny - 1] - this->yn[0]) / (this->ny - 1);
        lz_whole_avg = (this->zn[this->nz - 1] - this->zn[0]) / (this->nz - 1);
        leng_v0d1 = 0;
        leng_v0d1a = 0;
        v0d1num = 0;
        v0d1anum = 0;
        for (nodegs = 0; nodegs < node_group.size(); nodegs++) {
            for (auto ndi : node_group[nodegs]) {
                iz = ndi / (this->N_node_s);
                indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
                indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
                avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
                if (iz != 0) {    // this node is not on the bottom plane
                    eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = -1 / (this->zn[iz] - this->zn[iz - 1]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = -1 / (this->zn[iz] - this->zn[iz - 1]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                }
                if (iz != this->nz - 1) {   // this node is not on the top plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = 1 / (this->zn[iz + 1] - this->zn[iz]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = 1 / (this->zn[iz + 1] - this->zn[iz]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                }
                if (indx != 0) {    // this node is not on the left plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = -1 / (this->xn[indx] - this->xn[indx - 1]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = -1 / (this->xn[indx] - this->xn[indx - 1]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                }
                if (indx != this->nx - 1) {    // this node is not on the right plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = 1 / (this->xn[indx + 1] - this->xn[indx]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = 1 / (this->xn[indx + 1] - this->xn[indx]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                }
                if (indy != 0) {    // this node is not on the front plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = -1 / (this->yn[indy] - this->yn[indy - 1]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = -1 / (this->yn[indy] - this->yn[indy - 1]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                }
                if (indy != this->ny - 1) {   // this node is not on the back plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = 1 / (this->yn[indy + 1] - this->yn[indy]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0d1RowId[v0d1num] = eno;
                            this->v0d1ColId[v0d1num] = leng_v0d1;
                            this->v0d1val[v0d1num] = 1 / (this->yn[indy + 1] - this->yn[indy]);
                            v0d1num++;
                            this->v0d1aval[v0d1anum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0d1anum++;
                        }
                    }
                }
            }
            leng_v0d1++;
            leng_v0d1a++;
        }
        for (indi = 0; indi < this->numCdt; indi++) {
            if (this->conductor[indi].markPort == -1) {
                continue;
            }
            for (indj = 0; indj < this->cdtNumNode[indi]; indj++) {
                iz = this->conductor[indi].node[indj] / this->N_node_s;
                ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
                iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
                avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
                if (iz != 0) {    // this node is not on the bottom plane
                    eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
                    if (this->markEdge[eno] == 0 && this->markNode[(iz - 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this edge is in the dielectric
                        this->v0d1RowId[v0d1num] = eno;
                        this->v0d1ColId[v0d1num] = leng_v0d1;
                        this->v0d1val[v0d1num] = -1 / (this->zn[iz] - this->zn[iz - 1]);
                        v0d1num++;
                        this->v0d1aval[v0d1anum] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                        v0d1anum++;
                        mark = 1;
                    }
                }
                if (iz != this->nz - 1) {   // this node is not on the top plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
                    if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                        this->v0d1RowId[v0d1num] = eno;
                        this->v0d1ColId[v0d1num] = leng_v0d1;
                        this->v0d1val[v0d1num] = 1 / (this->zn[iz + 1] - this->zn[iz]);
                        v0d1num++;
                        this->v0d1aval[v0d1anum] = lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                        v0d1anum++;
                        mark = 1;
                    }
                }
                if (ix != 0) {    // this node is not on the left plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
                    if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                        this->v0d1RowId[v0d1num] = eno;
                        this->v0d1ColId[v0d1num] = leng_v0d1;
                        this->v0d1val[v0d1num] = -1 / (this->xn[ix] - this->xn[ix - 1]);
                        v0d1num++;
                        this->v0d1aval[v0d1anum] = -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                        v0d1anum++;
                        mark = 1;
                    }
                }
                if (ix != this->nx - 1) {    // this node is not on the right plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
                    if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                        this->v0d1RowId[v0d1num] = eno;
                        this->v0d1ColId[v0d1num] = leng_v0d1;
                        this->v0d1val[v0d1num] = 1 / (this->xn[ix + 1] - this->xn[ix]);
                        v0d1num++;
                        this->v0d1aval[v0d1anum] = ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                        v0d1anum++;
                        mark = 1;
                    }
                }
                if (iy != 0) {    // this node is not on the front plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
                    if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                        this->v0d1RowId[v0d1num] = eno;
                        this->v0d1ColId[v0d1num] = leng_v0d1;
                        this->v0d1val[v0d1num] = -1 / (this->yn[iy] - this->yn[iy - 1]);
                        v0d1num++;
                        this->v0d1aval[v0d1anum] = -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                        v0d1anum++;
                        mark = 1;
                    }
                }
                if (iy != this->ny - 1) {   // this node is not on the back plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
                    if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
                        this->v0d1RowId[v0d1num] = eno;
                        this->v0d1ColId[v0d1num] = leng_v0d1;
                        this->v0d1val[v0d1num] = 1 / (this->yn[iy + 1] - this->yn[iy]);
                        v0d1num++;
                        this->v0d1aval[v0d1anum] = lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                        v0d1anum++;
                        mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
                    }
                }
            }
            if (mark == 1) {
                leng_v0d1++;
                leng_v0d1a++;
            }
        }
        //ofstream out;
        //out.open("V0d.txt", std::ofstream::out | std::ofstream::trunc);
        //for (indi = 0; indi < v0d1num; indi++) {
        //    out << this->v0d1RowId[indi] << " " << this->v0d1ColId[indi] << " ";
        //    out << std::setprecision(15) << this->v0d1val[indi] << endl;
        //}
        //out.close();
        //out.open("V0da.txt", std::ofstream::out | std::ofstream::trunc);
        //for (indi = 0; indi < v0d1num; indi++) {
        //    out << this->v0d1RowId[indi] << " " << this->v0d1ColId[indi] << " ";
        //    out << std::setprecision(15) << this->v0d1aval[indi] << endl;
        //}
        //out.close();
        node_group.clear();
    }


    /* Generate Ad */
    void generateAd(myint *map, myint v0d1num, myint v0d1anum, myint leng_v0d1, myint& leng_Ad){
        unordered_map<myint, unordered_map<myint, double>> Ad1;
        myint indi, inz, inx, iny, node1, node2;

        for (indi = 0; indi < v0d1anum; indi++) {
            if (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // this edge is along z axis
                inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
                inx = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) / (this->N_cell_y + 1);
                iny = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) % (this->N_cell_y + 1);
                node1 = inz * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
                node2 = (inz + 1) * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
                if (map[node1] != this->v0d1ColId[indi] + 1 && map[node1] != 0) {
                    Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                }
                else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
                    Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                }
                else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0);
                }
            }
            else if (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
                inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
                inx = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
                iny = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
                node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
                node2 = inz * this->N_node_s + (inx + 1) * (this->N_cell_y + 1) + iny;
                if (map[node1] != this->v0d1ColId[indi] + 1 && map[node1] != 0) {
                    Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                }
                else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
                    Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                }
                else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0);
                }
            }
            else{    // this edge is along y axis
                inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
                inx = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
                iny = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
                node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
                node2 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny + 1;
                if (map[node1] != this->v0d1ColId[indi] + 1 && map[node1] != 0) {
                    Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                }
                else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
                    Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
                }
                else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
                    Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0);
                }
            }
        }


        for (indi = 0; indi < leng_v0d1; indi++) {
            leng_Ad += Ad1[indi].size();
        }

        this->AdRowId = (myint*)calloc(leng_Ad, sizeof(myint));
        this->AdColId = (myint*)calloc(leng_Ad, sizeof(myint));
        this->Adval = (double*)calloc(leng_Ad, sizeof(double));
        myint indj = 0;

        for (indi = 0; indi < leng_v0d1; indi++) {
            vector<pair<myint, double>> v(Ad1[indi].begin(), Ad1[indi].end());
            sort(v.begin(), v.end());
            for (auto adi : v) {
                if (abs(adi.second) > 1e-8) {
                    this->AdRowId[indj] = indi;
                    this->AdColId[indj] = adi.first;
                    this->Adval[indj] = adi.second;
                    indj++;
                }
            }
            v.clear();
        }
        Ad1.clear();

    }

    /* Generate V0c */
    void merge_v0c(double block_x, double block_y, double block2_x, double block2_y, myint &v0cnum, myint &leng_v0c, myint &v0canum, myint &leng_v0ca, myint *map) {

        int *visited;
        double ratio;
        double startx, starty;    // the start coordinates of each block
        queue<int> st;    // dfs stack
        //vector<int> ind;
        int indsize;
        int indx, indy;
        int mark;
        int markcond;
        int count = 0;
        leng_v0c = 0;
        leng_v0ca = 0;
        v0cnum = 0;
        v0canum = 0;
        int indnum;
        int ix, iy, iz;
        int n;
        int i, j;
        int map_count = 1;
        myint eno;
        double lx_avg, ly_avg, lz_avg;
        myint node1, node2;


        unordered_map<myint, double> v, va;
        visited = (int*)calloc(this->N_node, sizeof(int));
        vector<set<myint>> node_group;
        set<myint> base;
        int nodegs;

        for (int ic = 0; ic < this->numCdt; ic++) {
            if (this->conductor[ic].markPort <= 0) {    // not excited conductors
                markcond = ic + 1;
                //visited = (int*)calloc(this->N_node, sizeof(int));
                n = this->cdtNumNode[ic] - 1;
                if (this->conductor[ic].markPort <= -1)
                    n = this->cdtNumNode[ic];
                for (int jc = 0; jc < n; jc++) {
                    if (visited[this->conductor[ic].node[jc]] == 0 && this->ubdn.find(this->conductor[ic].node[jc]) == this->ubdn.end() && this->lbdn.find(this->conductor[ic].node[jc]) == this->lbdn.end()) {   // this node is not visited and it is not in the boundary


                        //if (!ind.empty())
                        //    ind.clear();
                        iz = this->conductor[ic].node[jc] / (this->N_node_s);
                        ix = (this->conductor[ic].node[jc] - iz * this->N_node_s) / (this->N_cell_y + 1);
                        iy = this->conductor[ic].node[jc] % (this->N_cell_y + 1);
                        startx = this->xn[ix];
                        starty = this->yn[iy];
                        //ind.push_back(this->conductor[ic].node[jc]);
                        st.push(ix * (this->N_cell_y + 1) + iy);
                        visited[this->conductor[ic].node[jc]] = 1;
                        map[this->conductor[ic].node[jc]] = map_count;
                        node_group.push_back(base);
                        nodegs = node_group.size() - 1;
                        node_group[nodegs].insert(this->conductor[ic].node[jc]);
                        while (!st.empty()) {

                            mark = 0;
                            indx = (st.front()) / (this->N_cell_y + 1);
                            indy = st.front() % (this->N_cell_y + 1);

                            if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
                                if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block_y) {    // this node is within the block area
                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (this->N_cell_x + 1) * indx + indy] == markcond
                                        && visited[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] == 0
                                        && (this->conductor[ic].markPort == 0 && iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy != this->conductor[ic].node[this->cdtNumNode[ic] - 1]) || (this->conductor[ic].markPort == -1 && iz != 0 && iz != this->nz - 1)) {    // this node is in conductor and this node is not visited, markPort == 0 then this node is not the last node, markPort == -1 then this node is not on the boundary


                                        st.push((indx + 1)*(this->N_cell_y + 1) + indy);
                                        visited[iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy] = 1;
                                        map[iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy);

                                    }
                                }
                            }
                            if (indx != 0) {    // it must have a left x edge, thus left x node
                                if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block_y) {    // this node is within the block area
                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (this->N_cell_x + 1) * (indx - 1) + indy] == markcond
                                        && visited[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] == 0
                                        && (this->conductor[ic].markPort == 0 && iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy != this->conductor[ic].node[this->cdtNumNode[ic] - 1]) || (this->conductor[ic].markPort == -1 && iz != 0 && iz != this->nz - 1)) {    // this node is in conductor and this node is not visited


                                        st.push((indx - 1)*(this->N_cell_y + 1) + indy);
                                        visited[iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy] = 1;
                                        map[iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy);
                                        //mark = 1;

                                        //continue;
                                    }
                                }
                            }
                            if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
                                if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block_y) {    // this node is within the block area
                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * indx + indy] == markcond
                                        && visited[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1] == 0
                                        && (this->conductor[ic].markPort == 0 && iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1 != this->conductor[ic].node[this->cdtNumNode[ic] - 1]) || (this->conductor[ic].markPort == -1 && iz != 0 && iz != this->nz - 1)) {    // this node is in conductor and this node is not visited


                                        st.push((indx)*(this->N_cell_y + 1) + indy + 1);
                                        visited[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1] = 1;
                                        map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1);
                                        //mark = 1;

                                        //continue;
                                    }
                                }
                            }
                            if (indy != 0) {    // it must have a closer y edge, thus closer y node
                                if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block_y) {    // this node is within the block area

                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * indx + indy - 1] == markcond
                                        && visited[iz * this->N_node_s + (indx)* (this->N_cell_y + 1) + indy - 1] == 0
                                        && (this->conductor[ic].markPort == 0 && iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1 != this->conductor[ic].node[this->cdtNumNode[ic] - 1]) || (this->conductor[ic].markPort == -1 && iz != 0 && iz != this->nz - 1)) {    // this node is in conductor and this node is not visited

                                        st.push((indx)*(this->N_cell_y + 1) + indy - 1);
                                        visited[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1] = 1;
                                        map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1);
                                        //mark = 1;

                                        //continue;
                                    }
                                }
                            }
                            //if (mark == 0) {
                            st.pop();
                            //}
                        }
                        for (auto ndi : node_group[nodegs]) {
                            indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
                            indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
                            avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
                            if (iz != 0) {    // this node is not on the bottom plane
                                eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (iz != this->nz - 1) {   // this node is not on the top plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indx != 0) {    // this node is not on the left plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indx != this->nx - 1) {    // this node is not on the right plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indy != 0) {    // this node is not on the front plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indy != this->ny - 1) {   // this node is not on the back plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                        }

                        map_count++;

                    }
                }

                if (leng_v0c > this->acu_cnno.back())
                    this->acu_cnno.push_back(leng_v0c);

                //free(visited); visited = NULL;
            }
            else{   // markPort > 0 then the conductor is not reference conductor, no need to judge the boundary nodes


                markcond = ic + 1;

                //visited = (int*)calloc(this->N_node, sizeof(int));
                n = this->cdtNumNode[ic] - 1;
                for (int jc = 0; jc < n; jc++) {
                    if (visited[this->conductor[ic].node[jc]] == 0) {

                        iz = this->conductor[ic].node[jc] / (this->N_node_s);
                        ix = (this->conductor[ic].node[jc] - iz * this->N_node_s) / (this->N_cell_y + 1);
                        iy = this->conductor[ic].node[jc] % (this->N_cell_y + 1);
                        startx = this->xn[ix];
                        starty = this->yn[iy];
                        //ind.push_back(this->conductor[ic].node[jc]);
                        st.push(ix * (this->N_cell_y + 1) + iy);
                        visited[this->conductor[ic].node[jc]] = 1;
                        map[this->conductor[ic].node[jc]] = map_count;
                        node_group.push_back(base);
                        nodegs = node_group.size() - 1;
                        node_group[nodegs].insert(this->conductor[ic].node[jc]);

                        while (!st.empty()) {

                            mark = 0;
                            indx = (st.front()) / (this->N_cell_y + 1);
                            indy = st.front() % (this->N_cell_y + 1);

                            if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
                                if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (this->N_cell_x + 1) * indx + indy] == markcond && visited[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] == 0 && (iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy != this->conductor[ic].node[this->cdtNumNode[ic] - 1] || this->conductor[ic].markPort == -1)) {    // this node is in conductor and this node is not visited


                                        st.push((indx + 1)*(this->N_cell_y + 1) + indy);
                                        visited[iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy] = 1;
                                        map[iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx + 1)*(this->N_cell_y + 1) + indy);
                                        //mark = 1;

                                        //continue;
                                    }
                                }
                            }
                            if (indx != 0) {    // it must have a left x edge, thus left x node
                                if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (this->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] == 0 && (iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy != this->conductor[ic].node[this->cdtNumNode[ic] - 1] || this->conductor[ic].markPort == -1)) {    // this node is in conductor and this node is not visited


                                        st.push((indx - 1)*(this->N_cell_y + 1) + indy);
                                        visited[iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy] = 1;
                                        map[iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx - 1)*(this->N_cell_y + 1) + indy);
                                        //mark = 1;

                                        //continue;
                                    }
                                }
                            }
                            if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
                                if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block2_y) {    // this node is within the block area
                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * indx + indy] == markcond && visited[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1] == 0 && (iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1 != this->conductor[ic].node[this->cdtNumNode[ic] - 1] || this->conductor[ic].markPort == -1)) {    // this node is in conductor and this node is not visited


                                        st.push((indx)*(this->N_cell_y + 1) + indy + 1);
                                        visited[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1] = 1;
                                        map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy + 1);
                                        //mark = 1;

                                        //continue;
                                    }
                                }
                            }
                            if (indy != 0) {    // it must have a closer y edge, thus closer y node
                                if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block2_y) {    // this node is within the block area
                                    if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * indx + indy - 1] == markcond && visited[iz * this->N_node_s + (indx)* (this->N_cell_y + 1) + indy - 1] == 0 && (iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1 != this->conductor[ic].node[this->cdtNumNode[ic] - 1] || this->conductor[ic].markPort == -1)) {    // this node is in conductor and this node is not visited

                                        st.push((indx)*(this->N_cell_y + 1) + indy - 1);
                                        visited[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1] = 1;
                                        map[iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1] = map_count;
                                        node_group[nodegs].insert(iz * this->N_node_s + (indx)*(this->N_cell_y + 1) + indy - 1);
                                        //mark = 1;

                                        //continue;
                                    }
                                }
                            }
                            //if (mark == 0) {

                            st.pop();

                            //}
                        }
                        for (auto ndi : node_group[nodegs]) {
                            indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
                            indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
                            avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
                            if (iz != 0) {    // this node is not on the bottom plane
                                eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (iz != this->nz - 1) {   // this node is not on the top plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indx != 0) {    // this node is not on the left plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indx != this->nx - 1) {    // this node is not on the right plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indy != 0) {    // this node is not on the front plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                            if (indy != this->ny - 1) {   // this node is not on the back plane
                                eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
                                compute_edgelink(eno, node1, node2);
                                if (node1 != ndi) {
                                    if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                                else if (node2 != ndi) {
                                    if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {

                                        v0cnum++;
                                        v0canum++;
                                    }
                                }
                            }
                        }
                        map_count++;


                    }
                }

                if (leng_v0c > this->acu_cnno.back())
                    this->acu_cnno.push_back(leng_v0c);

                //free(visited); visited = NULL;

            }
        }
        this->v0cRowId = (myint*)malloc(v0cnum * sizeof(myint));
        this->v0cColId = (myint*)malloc(v0cnum * sizeof(myint));
        this->v0cval = (double*)malloc(v0cnum * sizeof(double));
        //this->v0caRowId = (myint*)malloc(v0canum * sizeof(myint));
        //this->v0caColId = (myint*)malloc(v0canum * sizeof(myint));
        this->v0caval = (double*)malloc(v0canum * sizeof(double));
        v0cnum = 0;
        v0canum = 0;
        leng_v0c = 0;
        leng_v0ca = 0;
        double lx_whole_avg = 0;
        double ly_whole_avg = 0;
        double lz_whole_avg = 0;
        lx_whole_avg = (this->xn[this->nx - 1] - this->xn[0]) / (this->nx - 1);
        ly_whole_avg = (this->yn[this->ny - 1] - this->yn[0]) / (this->ny - 1);
        lz_whole_avg = (this->zn[this->nz - 1] - this->zn[0]) / (this->nz - 1);

        for (nodegs = 0; nodegs < node_group.size(); nodegs++) {
            for (auto ndi : node_group[nodegs]) {
                iz = ndi / (this->N_node_s);
                indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
                indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
                avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
                if (iz != 0) {    // this node is not on the bottom plane
                    eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = -1 / (this->zn[iz] - this->zn[iz - 1]);
                            v0cnum++;
                            this->v0caval[v0canum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = -1 / (this->zn[iz] - this->zn[iz - 1]);
                            v0cnum++;
                            this->v0caval[v0canum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                }
                if (iz != this->nz - 1) {   // this node is not on the top plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = 1 / (this->zn[iz + 1] - this->zn[iz]);
                            v0cnum++;
                            this->v0caval[v0canum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = 1 / (this->zn[iz + 1] - this->zn[iz]);
                            v0cnum++;
                            this->v0caval[v0canum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                }
                if (indx != 0) {    // this node is not on the left plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = -1 / (this->xn[indx] - this->xn[indx - 1]);
                            v0cnum++;
                            this->v0caval[v0canum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = -1 / (this->xn[indx] - this->xn[indx - 1]);
                            v0cnum++;
                            this->v0caval[v0canum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                }
                if (indx != this->nx - 1) {    // this node is not on the right plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = 1 / (this->xn[indx + 1] - this->xn[indx]);
                            v0cnum++;
                            this->v0caval[v0canum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = 1 / (this->xn[indx + 1] - this->xn[indx]);
                            v0cnum++;
                            this->v0caval[v0canum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                }
                if (indy != 0) {    // this node is not on the front plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = -1 / (this->yn[indy] - this->yn[indy - 1]);
                            v0cnum++;
                            this->v0caval[v0canum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = -1 / (this->yn[indy] - this->yn[indy - 1]);
                            v0cnum++;
                            this->v0caval[v0canum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                }
                if (indy != this->ny - 1) {   // this node is not on the back plane
                    eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
                    compute_edgelink(eno, node1, node2);
                    if (node1 != ndi) {
                        if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = 1 / (this->yn[indy + 1] - this->yn[indy]);
                            v0cnum++;
                            this->v0caval[v0canum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                    else if (node2 != ndi) {
                        if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
                            this->v0cRowId[v0cnum] = eno;
                            this->v0cColId[v0cnum] = leng_v0c;
                            this->v0cval[v0cnum] = 1 / (this->yn[indy + 1] - this->yn[indy]);
                            v0cnum++;
                            this->v0caval[v0canum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
                            v0canum++;
                        }
                    }
                }
            }

            leng_v0c++;
            leng_v0ca++;
        }
        //ofstream out;
        //out.open("V0c.txt", std::ofstream::out | std::ofstream::trunc);
        //for (int indi = 0; indi < v0cnum; indi++) {
        //    out << this->v0cRowId[indi] << " " << this->v0cColId[indi] << " ";
        //    out << std::setprecision (15) << this->v0cval[indi] << endl;
        //}
        //out.close();
        //out.open("V0ca.txt", std::ofstream::out | std::ofstream::trunc);
        //for (int indi = 0; indi < v0cnum; indi++) {
        //    out << this->v0cRowId[indi] << " " << this->v0cColId[indi] << " ";
        //    out << std::setprecision (15) << this->v0caval[indi] << endl;
        //}
        //out.close();
    }

    void generateAc(myint *map, myint v0cnum, myint v0canum, myint leng_v0c, myint& leng_Ac) {
        int indi, indj, inx, iny, inz;
        myint node1, node2;
        unordered_map<myint, unordered_map<myint, double>> Ac;

        for (indi = 0; indi < v0canum; indi++) {    // the upper and lower planes are PEC
            if (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // this edge is along z axis
                inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
                inx = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) / (this->N_cell_y + 1);
                iny = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) % (this->N_cell_y + 1);
                node1 = inz * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
                node2 = (inz + 1) * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
                if (map[node1] != this->v0cColId[indi] + 1 && map[node1] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][map[node1] - 1] += this->v0caval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;
                }
                else if (map[node2] != this->v0cColId[indi] + 1 && map[node2] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][map[node2] - 1] += this->v0caval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;
                }
                else if (this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += abs(this->v0caval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * SIGMA);
                }
            }
            else if (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
                inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
                inx = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
                iny = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
                node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
                node2 = inz * this->N_node_s + (inx + 1) * (this->N_cell_y + 1) + iny;
                if (map[node1] != this->v0cColId[indi] + 1 && map[node1] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][map[node1] - 1] += this->v0caval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;
                }
                else if (map[node2] != this->v0cColId[indi] + 1 && map[node2] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][map[node2] - 1] += this->v0caval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;
                }
                else if (this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += abs(this->v0caval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * SIGMA);
                }
            }
            else{    // this edge is along y axis
                inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
                inx = (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
                iny = (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
                node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
                node2 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny + 1;
                if (map[node1] != this->v0cColId[indi] + 1 && map[node1] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][map[node1] - 1] += this->v0caval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;
                }
                else if (map[node2] != this->v0cColId[indi] + 1 && map[node2] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][map[node2] - 1] += this->v0caval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;
                }
                else if (this->markEdge[this->v0cRowId[indi]] != 0) {
                    Ac[this->v0cColId[indi]][this->v0cColId[indi]] += abs(this->v0caval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * SIGMA);
                }
            }
        }

        for (indi = 0; indi < leng_v0c; indi++) {
            leng_Ac += Ac[indi].size();
        }
        this->AcRowId = (myint*)calloc(leng_Ac, sizeof(myint));
        this->AcColId = (myint*)calloc(leng_Ac, sizeof(myint));
        this->Acval = (double*)calloc(leng_Ac, sizeof(double));
        indj = 0;
        int k = 1;
        for (indi = 0; indi < leng_v0c; indi++) {
            vector<pair<myint, double>> v(Ac[indi].begin(), Ac[indi].end());
            sort(v.begin(), v.end());
            for (auto aci : v) {
                if (abs(aci.second) > 1e5) {
                    this->AcRowId[indj] = indi;
                    this->AcColId[indj] = aci.first;
                    this->Acval[indj] = aci.second;
                    if (this->AcRowId[indj] >= this->acu_cnno[k]) {
                        this->cindex.push_back(indj - 1);
                        k++;
                    }
                    indj++;
                }
            }
            v.clear();
        }
        this->cindex.push_back(indj - 1);
        Ac.clear();
    }

    /* Find Vh by using Arnoldi method */
    void find_Vh_Arnoldi(int k) {
        /* sourcePort : The port No. now considered
        k : Arnoldi run k steps
        */
        if ((this->N_edge - this->bden) * 2 < k) {
            k = (this->N_edge - this->bden) * 2;
        }
        // Arnoldi method
        double scale = 1.e+14;    // balance the matrix
        int i, j, indi, start;
        double* V = new double[(k + 1) * 2 * (this->N_edge - this->bden)]();   // Arnoldi orthogonal vectors
        double* w = new double[2 * (this->N_edge - this->bden)];
        double* H = new double[k * k]();   // Hessenberg matrix with (k + 1) rows and k columns, stored in column major and initial with 0
        double temp, cri_zero = 1e-3;

        for (i = 0; i < 2 * (this->N_edge - this->bden); i++){
            V[i] = 1 / sqrt(2 * (this->N_edge - this->bden));   // arbitrary starting vector
        }
        for (j = 1; j <= k; j++) {
            indi = 0;
            while (indi < this->leng_S) {
                // w = B^-1*A*v_j-1
                start = this->SRowId[indi];
                w[this->N_edge - this->bden + start] = 0;
                w[start] = V[(j - 1) * 2 * (this->N_edge - this->bden) + (this->N_edge - this->bden) + start];
                while (indi < this->leng_S && this->SRowId[indi] == start) {
                    w[this->N_edge - this->bden + start] += -this->Sval[indi] * V[(j - 1) * 2 * (this->N_edge - this->bden) + this->SColId[indi]] / (this->stackEpsn[(this->mapEdgeR[start] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0 * pow(scale, 2));
                    indi++;
                }
                if (this->markEdge[this->mapEdgeR[start]] != 0)
                    w[this->N_edge - this->bden + start] += -V[(j - 1) * 2 * (this->N_edge - this->bden) + (this->N_edge - this->bden) + start] * SIGMA / (this->stackEpsn[(this->mapEdgeR[start] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0 * scale);
            }
            for (i = 0; i <= j - 1; i++) {
                for (int in = 0; in < (this->N_edge - this->bden) * 2; in++) {
                    H[(j - 1) * k + i] += V[i * ((this->N_edge - this->bden) * 2) + in] * w[in];
                }
                for (int in = 0; in < (this->N_edge - this->bden) * 2; in++) {
                    w[in] -= H[(j - 1) * k + i] * V[i * ((this->N_edge - this->bden) * 2) + in];
                }
            }
            if (j == k) {
                temp = 0;
                for (i = 0; i < (this->N_edge - this->bden) * 2; i++) {
                    temp += pow(w[i], 2);
                }
                temp = sqrt(temp);
            }
            else {
                for (i = 0; i < (this->N_edge - this->bden) * 2; i++) {
                    H[(j - 1) * k + j] += pow(w[i], 2);
                }
                H[(j - 1) * k + j] = sqrt(H[(j - 1) * k + j]);
                temp = H[(j - 1) * k + j];
            }
            for (i = 0; i < (this->N_edge - this->bden) * 2; i++) {
                V[j * (this->N_edge - this->bden) * 2 + i] = w[i] / temp;
            }
            //if (temp < cri_zero) {    // if H[(j - 1) * k + j] is very small, Arnoldi terminates
            //    break;
            //}
        }


        double* H1 = new double[k * k]();   // copy of H, because dhseqr will change H
        for (i = 0; i < k; i++){
            for (j = 0; j < k; j++){
                if (i <= j + 1){
                    H1[j * k + i] = H[j * k + i];
                }
            }
        }
        ofstream out;
        //out.open("H.txt", std::ofstream::out | std::ofstream::trunc);
        //for (i = 0; i < k; i++){
        //    for (j = 0; j < k; j++){
        //        out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1) << H1[j * k + i] << " ";
        //    }
        //    out << endl;
        //}
        //out.close();

        //out.open("V.txt", std::ofstream::out | std::ofstream::trunc);
        //for (i = 0; i < 2 * (this->N_edge - this->bden); i++) {
        //    for (j = 0; j < k + 1; j++) {
        //        out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1) << V[j * 2 * (this->N_edge - this->bden) + i] << " ";
        //    }
        //    out << endl;
        //}
        //out.close();


        /* LAPACKE_dgeev to find the eigenvalues and eigenvectors of H */
        //// make the matrix balanced
        //int matrix_layout = LAPACK_COL_MAJOR;
        //lapack_int info;
        //char job = 'B';
        //lapack_int n = k;    // the order of matrix H
        //lapack_int lda = n;    // leading dimension of H
        //lapack_int ilo, ihi;
        //double* scale; 
        //scale = new double[n];
        //info = LAPACKE_dgebal(matrix_layout, job, n, H, lda, &ilo, &ihi, scale);



        //char jobvl = 'N';    // left eigenvectors of H are not computed
        //char jobvr = 'V';    // right eigenvectors of H are computed
        //double* wr, * wi;    // real and imaginary eigenvalues
        //double* vl, * vr;    // eigenvectors
        //lapack_int ldvl = n, ldvr = n;

        //wr = new double[n];
        //wi = new double[n];
        //vr = new double[n * n];

        //
        //info = LAPACKE_dgeev(matrix_layout, jobvl, jobvr, n, H, lda, wr, wi, vl, ldvl, vr, ldvr);

        //if (info != 0) {
        //	cout << "LAPACKE_dgeev is not successful and info is " << info << endl;
        //	return;
        //}
        //lapack_complex_double* v = new lapack_complex_double[n * n];
        //for (i = 0; i < n; i++) {
        //	for (j = 0; j < n; j++) {
        //		if (wi[j] == 0) {
        //			v[j * n + i].real = vr[j * n + i];
        //		}
        //		else {
        //			v[j * n + i].real = vr[j * n + i];
        //			v[j * n + i].imag = vr[(j + 1) * n + i];
        //			v[(j + 1) * n + i].real = vr[j * n + i];
        //			v[(j + 1) * n + i].imag = -vr[(j + 1) * n + i];
        //			j++;
        //		}
        //	}
        //}

        //// transform eigenvectors of the balanced matrix to the original
        ///*char side = 'R';
        //lapack_int m = n;
        //lapack_int ldv = n;
        //info = LAPACKE_zgebak(matrix_layout, job, side, n, ilo, ihi, scale, m, v, ldv);*/

        

        //out.open("v.txt", std::ofstream::out | std::ofstream::trunc);
        //for (i = 0; i < n; i++) {
        //	for (j = 0; j < n; j++) {
        //		out << v[j * n + i].real << " " << v[j * n + i].imag << " ";
        //	}
        //	out << endl;
        //}


        //delete[] wr;
        //delete[] wi;
        //delete[] vr;
        //      delete[] v;
        //delete[] scale;

        /* mkl to find the eigenvalues and eigenvectors of H (not accurate) */
        int matrix_layout = LAPACK_COL_MAJOR;
        char job = 'E';   // only eigenvalues are required
        char compz = 'N';   // no Schur vectors are computed
        lapack_int n = k;    // the order of the matrix
        lapack_int ilo = 1, ihi = n;
        double* z;    // compz = 'N' and z need not to be set
        lapack_int ldh = n;    // the leading dimension of H
        double* wr, *wi;    // real and imaginary parts of eigenvalues, initialize ?
        lapack_int ldz = n;    // if compz = 'N' then ldz >= 1
        double eps = 1e-6;    // esp * ma is considered to be zero eigenvalues

        wr = new double[n];
        wi = new double[n];
        lapack_int info;
        info = LAPACKE_dhseqr(matrix_layout, job, compz, n, ilo, ihi, H, ldh, wr, wi, z, ldz);    // calculate all the eigenvalues of H



        char side = 'R';    // only right eigenvectors are computed
        char eigsrc = 'Q';    // the eigenvalues of H are from hseqr
        char initv = 'N';   // no initial estimates for eigenvectors
        lapack_logical* select = new lapack_logical[n];   // specify which eigenvectors are computed
        double* vl, *vr;   // initv = 'N' need not to be set
        lapack_int ldvl = n;   // leading dimension of vl
        lapack_int ldvr = n;    // leading dimension of vr
        lapack_int mm = 0;
        double* wrc, *wic;    // copy of wr and wc because after dhsein wr will be modified, eigenvalues of the original system
        wrc = new double[n];
        wic = new double[n];
        for (i = 0; i < n; i++) {
            wrc[i] = wr[i] * scale;
            wic[i] = wi[i] * scale;
            if (abs(wi[i]) > eps){
                select[i] = 1;    // calculate the non-zero eigenvalues' eigenvectors
                //cout << wr[i] << " " << wi[i] << endl;
                mm++;
            }
            else{
                select[i] = 0;
                //cout << wr[i] << " " << wi[i] << endl;
            }
        }
        out.open("w.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < n; i++) {
            out << wr[i] * scale << " " << wi[i] * scale << endl;
        }
        out.close();
        lapack_int* m = new lapack_int(mm);
        lapack_int* ifaill = new lapack_int[mm];   // 0 for converge
        lapack_int* ifailr = new lapack_int[mm];
        vr = new double[ldvr * mm];
        /*out.open("H1.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < k; i++){
        for (j = 0; j < k; j++){
        out << H1[j * k + i] << " ";
        }
        out << endl;
        }
        out.close();*/
        info = LAPACKE_dhsein(matrix_layout, side, eigsrc, initv, select, n, H1, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, ifaill, ifailr);

        for (i = 0; i < mm; i++) {
            if (ifailr[i] > 0)
                cout << "The " << i << "th column failed to converge and eigenvalue is " << wr[i] << " + 1i * " << wi[i] << "!" << endl;
        }

        //out.open("vr.txt", std::ofstream::out | std::ofstream::trunc);
        lapack_complex_double* v = new lapack_complex_double[k * mm];
        for (i = 0; i < n; i++) {
            int temp = 0;;
            for (j = 0; j < n; j++) {
                if (select[j] == 1){   // always generate the complex conjugate pairs eigenvectors
                    //out << vr[temp * n + i] << " " << vr[(temp + 1) * n + i] << " ";
                    //out << vr[temp * n + i] << " " << -vr[(temp + 1) * n + i] << " ";
                    //if (i == 0)
                    //    cout << wrc[j] << " " << wic[j] << endl;
                    v[temp * n + i].real = vr[temp * n + i];
                    v[temp * n + i].imag = vr[(temp + 1) * n + i];
                    v[(temp + 1) * n + i].real = vr[temp * n + i];
                    v[(temp + 1) * n + i].imag = -vr[(temp + 1) * n + i];
                    temp++;
                    temp++;
                    j++;
                    //if (i == 0)
                    //    cout << wrc[j] << " " << wic[j] << endl;
                }
            }
            //out << endl;
        }
        //out.close();
        delete[] H1;
        delete[] wr;
        delete[] wi;
        delete[] wrc;
        delete[] wic;
        delete[] vl;
        delete[] vr;
        delete[] select;
        delete[] m;
        delete[] ifaill;
        delete[] ifailr;

        /* Generate Vh = V * vr */
        this->Vh = new lapack_complex_double[(this->N_edge - this->bden) * mm];
        for (i = 0; i < mm; i++){
            for (j = 0; j < (this->N_edge - this->bden); j++){
                this->Vh[i * (this->N_edge - this->bden) + j].real = 0;
                this->Vh[i * (this->N_edge - this->bden) + j].imag = 0;
                for (int in = 0; in < k; in++){
                    this->Vh[i * (this->N_edge - this->bden) + j].real += V[in * (this->N_edge - this->bden) * 2 + j] * v[i * k + in].real;
                    this->Vh[i * (this->N_edge - this->bden) + j].imag += V[in * (this->N_edge - this->bden) * 2 + j] * v[i * k + in].imag;
                }
            }
        }
        this->leng_Vh = mm;
        cout << "There are " << this->leng_Vh << " number of columns in Vh!" << endl;
        delete[] v;



        delete[] w;
        delete[] V;
        delete[] H;

    }

    /* Reference eigenvectors for the generalized eigenvalue problem [-S, 0; 0, D_eps]v = \lambda [D_sig, D_eps; D_eps, 0]v */
    void find_reference_Vh() {
        double scale = 1e+14;
        double* A, *B;
        myint indi = 0, indj;

        A = new double[4 * (this->N_edge - this->bden) * (this->N_edge - this->bden)]();
        B = new double[4 * (this->N_edge - this->bden) * (this->N_edge - this->bden)]();

        // set A
        while (indi < this->leng_S) {
            A[this->SColId[indi] * (this->N_edge - this->bden) * 2 + this->SRowId[indi]] = -this->Sval[indi];
            indi++;
        }
        indi = 0;
        while (indi < (this->N_edge - this->bden)) {
            A[(this->N_edge - this->bden + indi) * (this->N_edge - this->bden) * 2 + this->N_edge - this->bden + indi] = pow(scale, 2) * this->stackEpsn[(this->mapEdgeR[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
            indi++;
        }

        // set B
        indi = 0;
        while (indi < (this->N_edge - this->bden)) {
            if (this->markEdge[this->mapEdgeR[indi]] != 0) {
                B[indi * (this->N_edge - this->bden) * 2 + indi] = scale * SIGMA;
            }
            B[(this->N_edge - this->bden + indi) * (this->N_edge - this->bden) * 2 + indi] = pow(scale, 2) * this->stackEpsn[(this->mapEdgeR[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
            B[indi * (this->N_edge - this->bden) * 2 + this->N_edge - this->bden + indi] = pow(scale, 2) * this->stackEpsn[(this->mapEdgeR[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0;
            indi++;
        }
        
        lapack_int n;
        char jobvl, jobvr;
        lapack_int lda, ldb, ldvl, ldvr;
        double *alphar, *alphai, *beta, *vl, *vr, *lscale, *rscale;
        n = (this->N_edge - this->bden) * 2;
        jobvl = 'N'; jobvr = 'V';
        lda = n; ldb = n;
        ldvl = n; ldvr = n;
        alphar = new double[n];
        alphai = new double[n];
        beta = new double[n];
        vl = new double[n * n];
        vr = new double[n * n];
        lscale = new double[n];
        rscale = new double[n];

        cout << "Begin to solve the generalized eigenvalue problem!\n";
        lapack_int info = LAPACKE_dggev(LAPACK_COL_MAJOR, jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr);
        if (info != 0) {
            cout << "Generalized eigenvalue problem ddgev is failed!\n";
            return;
        }
        cout << "Finish solving the generalized eigenvalue problem!\n";
        double eps = 1e-4;
        this->leng_Vh = 0;
        ofstream out;
        //out.open("ww.txt", std::ofstream::out | std::ofstream::trunc);
        for (indi = 0; indi < n; indi++) {
            if (abs(alphai[indi] / beta[indi]) > eps && abs(alphar[indi] / beta[indi]) < 1) {   // if the imaginary part of Vh is not so small
                //out << alphar[indi] / beta[indi] << " " << alphai[indi] / beta[indi] << endl;
                this->leng_Vh++;
            }
        }
        //out.close();
        out.open("vv.txt", std::ofstream::out | std::ofstream::trunc);
        this->Vh = new lapack_complex_double[(this->N_edge - this->bden) * this->leng_Vh];
        indi = 0; indj = 0;
        while (indi < n) {
            if (abs(alphai[indi] / beta[indi]) > eps && abs(alphar[indi] / beta[indi]) < 1) {    // the eigenvalue should have a small real part and a non-zero imaginary part
                for (int indk = 0; indk < this->N_edge - this->bden; indk++) {
                    this->Vh[indj * (this->N_edge - this->bden) + indk].real = vr[indi * (this->N_edge - this->bden) * 2 + indk];
                    this->Vh[indj * (this->N_edge - this->bden) + indk].imag = vr[(indi + 1) * (this->N_edge - this->bden) * 2 + indk];
                    this->Vh[(indj + 1) * (this->N_edge - this->bden) + indk].real = vr[indi * (this->N_edge - this->bden) * 2 + indk];
                    this->Vh[(indj + 1) * (this->N_edge - this->bden) + indk].imag = -vr[(indi + 1) * (this->N_edge - this->bden) * 2 + indk];
                }
                indj++;
                indj++;
                indi++;
            }
            indi++;
        }
        for (indi = 0; indi < this->N_edge - this->bden; indi++) {
            for (indj = 0; indj < this->leng_Vh; indj++) {
                out << this->Vh[indj * (this->N_edge - this->bden) + indi].real << " " << this->Vh[indj * (this->N_edge - this->bden) + indi].imag << " ";
            }
            out << endl;
        }
        out.close();
        cout << "Vh's number of columns is " << this->leng_Vh << endl;
        delete[] A;
        delete[] B;
        delete[] alphar;
        delete[] alphai;
        delete[] beta;
        delete[] vl;
        delete[] vr;
        delete[] lscale;
        delete[] rscale;
    }

	/* Generate backward difference left hand matrix */
	void backDiffLeftMatrix(myint* rowId, myint* colId, double* val_source, int leng, double* val, double dt) {
		/* (D_eps + dt * D_sig + dt^2 * S)
		   val : this matrix's val
		   dt : the time step
		   RowId is S's rowId, ColId is S's colId */
		int i = 0;
        
		while (i < leng) {
			val[i] = val_source[i] * pow(dt, 2);
			if (rowId[i] == colId[i]) {
				val[i] += this->stackEpsn[(this->mapEdgeR[rowId[i]] + this->N_edge_v) / (this->N_edge_v + this->N_edge_s)] * EPSILON0;
				if (this->markEdge[this->mapEdgeR[rowId[i]]]) {
					val[i] += dt * SIGMA;
				}
			}
            i++;
		}
	}

    /* Calculate the reference */
    void reference1(int freqNo, int sourcePort, complex<double>* xr) {
        double freq = this->freqNo2freq(freqNo);
        myint size = this->N_edge - this->bden;
        myint *RowId1 = (myint*)malloc((size + 1) * sizeof(myint));
        int count = 0;
        int indi = 0;
        int k = 0;
        complex<double> *valc;
        valc = (complex<double>*)calloc(this->leng_S, sizeof(complex<double>));
        complex<double> *J;
        J = (complex<double>*)calloc((this->N_edge - this->bden), sizeof(complex<double>));
        int indz, indy, temp;

        for (int sourcePortSide = 0; sourcePortSide < this->portCoor[sourcePort].multiplicity; sourcePortSide++) {
            for (int inde = 0; inde < this->portCoor[sourcePort].portEdge[sourcePortSide].size(); inde++){
                J[this->mapEdge[this->portCoor[sourcePort].portEdge[sourcePortSide][inde]]] = 0. - (1i) * (this->portCoor[sourcePort].portDirection[sourcePortSide] * 1.0) * freq * 2. * M_PI;
                //cout << "SourcePort " << sourcePort << " sourcePortSide " << sourcePortSide << " is " << this->portCoor[sourcePort].portDirection[sourcePortSide] << endl;

            }
        }
    
        

        /* Used in plasma2D for upper and lower excitation */
        /*myint current_edge = this->portEdge[sourcePort][indi - 1] + (this->N_edge_s + this->N_edge_v);

        while (current_edge < this->N_edge - this->N_edge_s){
        if (this->markEdge[current_edge] == 0){
        J[current_edge - this->N_edge_s] = 0. + (1i) * this->portCoor[sourcePort].portDirection * freq * 2. * M_PI;
        }
        current_edge = current_edge + (this->N_edge_s + this->N_edge_v);
        }*/
        /* end of Used in plasma2D for upper and lower excitation */

        RowId1[k] = 0;
        k++;
        myint start;
        myint nnz = this->leng_S;
        //cout << "Start to generate CSR form for S!\n";
        indi = 0;
        //ofstream out;
        //out.open("A1.txt", std::ofstream::out | std::ofstream::trunc);
        while (indi < nnz){
            start = this->SRowId[indi];
            while (indi < nnz && this->SRowId[indi] == start) {
                valc[indi] += this->Sval[indi]; // val[indi] is real
                if (this->SRowId[indi] == this->SColId[indi]){
                    if (this->markEdge[this->mapEdgeR[this->SRowId[indi]]] != 0) {
                        complex<double> addedPart(-(2. * M_PI * freq) * this->stackEpsn[(this->mapEdgeR[this->SRowId[indi]] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0, SIGMA);
                        valc[indi] += (2. * M_PI * freq) * addedPart;
                    }
                    else {
                        complex<double> addedPart(-(2. * M_PI * freq) * this->stackEpsn[(this->mapEdgeR[this->SRowId[indi]] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] * EPSILON0, 0);
                        valc[indi] += (2. * M_PI * freq) * addedPart;
                    }
                }
                //out << valc[indi].real() << " " << valc[indi].imag() << endl;
                count++;
                indi++;
            }
            RowId1[k] = (count);
            k++;
        }
        //out.close();
        //cout << "(-w^2*D_eps+iw*D_sig+S) is generated!\n" << endl;
        myint mtype = 13;    /* Real complex unsymmetric matrix */
        myint nrhs = 1;    /* Number of right hand sides */
        void *pt[64];

        /* Pardiso control parameters */
        myint iparm[64];
        myint maxfct, mnum, phase, error, msglvl, solver;
        double dparm[64];
        int v0csin;
        myint perm;

        /* Auxiliary variables */
        char *var;

        msglvl = 0;    /* print statistical information */
        solver = 0;    /* use sparse direct solver */
        error = 0;
        maxfct = 1;
        mnum = 1;
        phase = 13;

        pardisoinit(pt, &mtype, iparm);
        iparm[38] = 1;
        iparm[34] = 1;    // 0-based indexing
        iparm[3] = 2;    // number of processors
        //iparm[59] = 2;    // out of core version to solve very large problem
        //iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */

        //cout << "Begin to solve (-w^2*D_eps+iwD_sig+S)x=-iwJ\n";
        complex<double> * ddum;

        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, valc, RowId1, this->SColId, &perm, &nrhs, iparm, &msglvl, J, xr, &error);
        if (error != 0){
            printf("\nERROR during numerical factorization: %d", error);
            exit(2);
        }

        phase = -1;     // Release internal memory
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, &ddum, RowId1, this->SColId, &perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

        free(RowId1); RowId1 = NULL;
        free(valc); valc = NULL;
        free(J); J = NULL;

    }

    /* Calculate the averaged length */
    void avg_length(int iz, int iy, int ix, double &lx, double &ly, double &lz) {    // given a node, we can know its averaged lengths along x, y, z directions
        if (iz == 0) {
            lz = this->zn[1] - this->zn[0];
        }
        else if (iz == this->nz - 1) {
            lz = this->zn[iz] - this->zn[iz - 1];
        }
        else {
            lz = (this->zn[iz + 1] - this->zn[iz - 1]) / 2;
        }

        if (iy == 0) {
            ly = this->yn[1] - this->yn[0];
        }
        else if (iy == this->ny - 1) {
            ly = this->yn[iy] - this->yn[iy - 1];
        }
        else {
            ly = (this->yn[iy + 1] - this->yn[iy - 1]) / 2;
        }

        if (ix == 0) {
            lx = this->xn[1] - this->xn[0];
        }
        else if (ix == this->nx - 1) {
            lx = this->xn[ix] - this->xn[ix - 1];
        }
        else {
            lx = (this->xn[ix + 1] - this->xn[ix - 1]) / 2;
        }

    }

    /* Compute edgelink */
    void compute_edgelink(myint eno, myint &node1, myint &node2) {
        myint inz, inx, iny;
        if (eno % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // this edge is along z axis
            inz = eno / (this->N_edge_s + this->N_edge_v);
            inx = ((eno % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) / (this->N_cell_y + 1);
            iny = ((eno % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) % (this->N_cell_y + 1);
            node1 = inz * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
            node2 = (inz + 1) * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
        }
        else if (eno % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
            inz = eno / (this->N_edge_s + this->N_edge_v);
            inx = ((eno % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
            iny = ((eno % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
            node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
            node2 = inz * this->N_node_s + (inx + 1) * (this->N_cell_y + 1) + iny;
        }
        else {    // this edge is along y axis
            inz = eno / (this->N_edge_s + this->N_edge_v);
            inx = (eno % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
            iny = (eno % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
            node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
            node2 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny + 1;
        }
    }

    /* Compute one node's indx, indy, indz */
    void compute_node_index(myint node, int &indx, int &indy, int &indz) {
        /* node : the node number before adding the boundary condition
           indx : the node's x index 
           indy : the node's y index
           indz : the node's z index */
        indz = node / this->N_node_s;
        indx = (node % this->N_node_s) / (this->N_cell_y + 1);
        indy = (node % this->N_node_s) % (this->N_cell_y + 1);


    }

    /* Construct Z parameters with V0 and Vh */
    void Construct_Z_V0_Vh(complex<double> *x, int freqNo, int sourcePort){
        /* x: field distribution
        freqNo: frequency no.
        sourcePort: port no.
        Note: portDirection is the relative position of the port to the ground. E.g., ground is on the top, then portDirection = -1 */
        int inz, inx, iny;
        double leng;

        /* Only divide matrix entry by current at end of response port calculation */
        double sourceCurrent = 0.; // In-phase current from unit source port edge current densities into supply point (A)
        for (int sourcePortSide = 0; sourcePortSide < this->portCoor[sourcePort].multiplicity; sourcePortSide++)
        {
            sourceCurrent += this->portCoor[sourcePort].portArea[sourcePortSide];
        }

        for (int indPort = 0; indPort < this->numPorts; indPort++) {
            int indPortSide = 0; // Only deal with first port side to get response edge line integral
            for (int indEdge = 0; indEdge < this->portCoor[indPort].portEdge[indPortSide].size(); indEdge++) {
                myint thisEdge = this->portCoor[indPort].portEdge[indPortSide][indEdge];
                if (thisEdge % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // This edge is along the z-axis
                    inz = thisEdge / (this->N_edge_s + this->N_edge_v);
                    leng = this->zn[inz + 1] - this->zn[inz];
                }
                else if (thisEdge % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // This edge is along the x-axis
                    inx = ((thisEdge % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
                    leng = this->xn[inx + 1] - this->xn[inx];
                }
                else {    // This edge is along the y-axis
                    iny = (thisEdge % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
                    leng = this->yn[iny + 1] - this->yn[iny];
                }

                this->x[freqNo * (this->numPorts * this->numPorts) + indPort + this->numPorts * sourcePort] -= x[this->mapEdge[thisEdge]] * leng * (this->portCoor[indPort].portDirection[indPortSide] * 1.0); // Accumulating responses due to each response edge line integral (V)
            }


            this->x[freqNo * (this->numPorts * this->numPorts) + indPort + this->numPorts * sourcePort] /= sourceCurrent; // Final matrix entry (ohm)
        }
    }

    /* Construct Z parameters with V0 */
    void Construct_Z_V0(complex<double> *x, int sourcePort){
        /* x: field distribution from V0 solution
        sourcePort: port no. */
        myint inz, inx, iny;
        double leng;

        /* Z_ij = V_i / I_j = - integ[(grad V_i - part A_i / part t) * dl_i] / iinteg[(J_j) * dS_j] = - sum[(yd + yh)_respedge * leng_respedge]_oneside / sum[1 * area_source_side] */
        double sourceCurrent = 0.; // In-phase current from unit source port edge current densities into supply point (A)
        for (int sourcePortSide = 0; sourcePortSide < this->portCoor[sourcePort].multiplicity; sourcePortSide++)
        {
            sourceCurrent += this->portCoor[sourcePort].portArea[sourcePortSide];
        }

        for (int indPort = 0; indPort < this->numPorts; indPort++) {
            int indPortSide = 0; // Only deal with first port side to get response edge line integral
            for (int indEdge = 0; indEdge < this->portCoor[indPort].portEdge[indPortSide].size(); indEdge++) {
                myint thisEdge = this->portCoor[indPort].portEdge[indPortSide][indEdge];
                if (thisEdge % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // This edge is along the z-axis
                    inz = thisEdge / (this->N_edge_s + this->N_edge_v);
                    leng = this->zn[inz + 1] - this->zn[inz];
                }
                else if (thisEdge % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // This edge is along the x-axis
                    inx = ((thisEdge % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
                    leng = this->xn[inx + 1] - this->xn[inx];
                }
                else {    // This edge is along the y-axis
                    iny = (thisEdge % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
                    leng = this->yn[iny + 1] - this->yn[iny];
                }

                /*leng = pow((this->nodepos[this->edgelink[thisEdge * 2] * 3] - this->nodepos[this->edgelink[thisEdge * 2 + 1] * 3]), 2);
                leng += pow((this->nodepos[this->edgelink[thisEdge * 2] * 3 + 1] - this->nodepos[this->edgelink[thisEdge * 2 + 1] * 3 + 1]), 2);
                leng += pow((this->nodepos[this->edgelink[thisEdge * 2] * 3 + 2] - this->nodepos[this->edgelink[thisEdge * 2 + 1] * 3 + 2]), 2);
                leng = sqrt(leng);*/
                this->x[indPort + this->numPorts * sourcePort] -= x[thisEdge] * leng * (this->portCoor[indPort].portDirection[indPortSide] * 1.0); // Accumulating responses due to each response edge line integral (V)
            }

            /* Only divide matrix entry by current at end of response port calculation */
            //cout << "  leng = " << leng << ", first side portArea = " << this->portCoor[sourcePort].portArea[0] << " m^2, first side portDirection = " << this->portCoor[sourcePort].portDirection[0] << endl;
            //cout << "  Response port voltage = " << this->x[indPort + this->numPorts * xcol] << " V, Source port current = " << sourceCurrent << " A" << endl;
            this->x[indPort + this->numPorts * sourcePort] /= sourceCurrent; // Final matrix entry (ohm)
        }
    }

    /* print Z for both V0 and Vh */
    void print_z_V0_Vh(){
        int indi, inde, indj;
        double freq;

        for (indi = 0; indi < this->nfreq; indi++){
            // this point's frequency

            if (this->nfreq == 1) {    // to avoid (sys->nfreq - 1)
                freq = this->freqStart * this->freqUnit;
            }
            else {
                if (this->freqScale == 1) {
                    freq = (this->freqStart + indi * (this->freqEnd - this->freqStart) / (this->nfreq - 1)) * this->freqUnit;
                }
                else {
                    freq = this->freqStart * this->freqUnit * pow(this->freqEnd / this->freqStart, (indi * 1.0 / (this->nfreq - 1)));
                }
            }
            cout << "Z-parameters at frequency " << freq << " Hz:" << endl;

            for (inde = 0; inde < this->numPorts; inde++){
                for (indj = 0; indj < this->numPorts; indj++){   // One port excitation, different ports response as one row
                    cout << this->x[indi * (this->numPorts * this->numPorts) + inde * this->numPorts + indj] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }

    /* print Z for just V0 */
    void print_z_V0(){
        int indi, indj;
        double freq;
        complex<double> Zresult;

        if (this->nfreq > 1) {

            for (int id = 0; id < this->nfreq; id++) {
                if (id == 0)
                {
                    // First frequency in sweep
                    freq = this->freqStart * this->freqUnit;

                    // Report the saved result
                    cout << "Z-parameters at frequency " << (this->freqStart + id * (this->freqEnd - this->freqStart) / (this->nfreq - 1)) * this->freqUnit << " Hz:" << endl;
                    for (indi = 0; indi < this->numPorts; indi++) {
                        for (indj = 0; indj < this->numPorts; indj++) {
                            Zresult = this->x[indj + indi*this->numPorts];
                            cout << "  " << Zresult;
                        }
                        cout << endl;
                    }
                    continue;
                }
                else if (id == this->nfreq - 1)
                {
                    // Last frequency in sweep
                    freq = this->freqEnd * this->freqUnit;
                }
                else
                {
                    // All other frequencies in sweep
                    if (this->freqScale == 1)
                    {
                        // Linear interpolation of frequency sweep
                        freq = (this->freqStart + id * (this->freqEnd - this->freqStart) / (this->nfreq - 1)) * this->freqUnit;
                    }
                    else
                    {
                        // Logarithmic interpolation of frequency sweep
                        freq = this->freqStart * this->freqUnit * pow(this->freqEnd / this->freqStart, (id * 1.0 / (this->nfreq - 1))); // Should be most numerically stable calculated like this
                        //cout << "Log freq interp: " << freq << " Hz" << endl;
                    }
                }

                // Report the results beyond the first and append to storage object
                cout << "Z-parameters at frequency " << freq << " Hz:" << endl;
                for (indi = 0; indi < this->numPorts; indi++) {
                    for (indj = 0; indj < this->numPorts; indj++) {
                        Zresult = this->x[indj + indi*this->numPorts].real() + (1i) * this->x[indj + indi*this->numPorts].imag() * this->freqStart * this->freqUnit / freq;
                        cout << "  " << Zresult;
                        this->x[id * (this->numPorts * this->numPorts) + indi + indj * this->numPorts] = Zresult;
                    }
                    cout << endl;
                }
            }
        }
        else{
            cout << "Z-parameters at single frequency " << (this->freqStart) * this->freqUnit << " Hz:" << endl;
            for (indi = 0; indi < this->numPorts; indi++) {
                for (indj = 0; indj < this->numPorts; indj++) {
                    Zresult = this->x[indj + indi*this->numPorts];
                    cout << Zresult << " ";
                    //cout << Zresult.real() << "+ 1i* " << Zresult.imag() << " "; // Alternative for copying and pasting into MATLAB
                }
                cout << endl;
            }
        }
    }

    /* freqNo to freq */
    double freqNo2freq(int freqNo){
        double freq;

        if (this->nfreq == 1) {    // to avoid (sys->nfreq - 1)
            freq = this->freqStart * this->freqUnit;
        }
        else {
            if (this->freqScale == 1) {
                freq = (this->freqStart + freqNo * (this->freqEnd - this->freqStart) / (this->nfreq - 1)) * this->freqUnit;
            }
            else {
                freq = this->freqStart * this->freqUnit * pow(this->freqEnd / this->freqStart, (freqNo * 1.0 / (this->nfreq - 1)));
            }
        }
        return freq;
    }

	void generateLaplacianLeft(myint* LrowId, myint* LcolId, double* Lval, myint Lleng, double dt) {
		/* Generate the matrix [V0a'*D*V0, V0a'*D; D*V0, D+L] */
		myint i;
		int indx, indy, indz, mark;
		leng_Ld = 0;
		set<myint> remNode;   // the node which are removed due to PEC
		if (lbdn.size() == 0 && ubdn.size() == 0)
			remNode.insert(N_node);
		else if (lbdn.size() != 0)
			remNode = lbdn;
		else
			remNode = ubdn;
		

		for (i = 0; i < N_node; ++i) {
			if (remNode.find(i) != remNode.end())
				continue;
			compute_node_index(i, indx, indy, indz);
			leng_Ld++;   // with itself
			if (indz > 0) {   // has a lower node
				leng_Ld++;
				leng_Ld++;    // V0a'*D has one non-zero edge
				leng_Ld++;    // D*V0 has one non-zero edge
			}
			if (indz < nz - 1) {   // has an upper node
				leng_Ld++;
				leng_Ld++;
				leng_Ld++;
			}
			if (indx > 0) {   // has a left node
				leng_Ld++;
				leng_Ld++;
				leng_Ld++;
			}
			if (indx < nx - 1) {    // has a right node
				leng_Ld++;
				leng_Ld++;
				leng_Ld++;
			}
			if (indy > 0) {    // has a front node
				leng_Ld++;
				leng_Ld++;
				leng_Ld++;
			}
			if (indy < ny - 1) {    // has a back node
				leng_Ld++;
				leng_Ld++;
				leng_Ld++;
			}
		}
		leng_Ld += Lleng;    // add the non-zero entries in (D_eps+dt*D_sig+dt^2*L)

		// balance the norm of the matrices, normalize the V0 vectors, columnwise matrix
		LdrowId = (myint*)malloc(leng_Ld * sizeof(myint));
		LdcolId = (myint*)malloc(leng_Ld * sizeof(myint));
		Ldval = (double*)malloc(leng_Ld * sizeof(double));
		leng_Ld = 0;
		for (i = 0; i < N_node; ++i) {
			compute_node_index(i, indx, indy, indz);
			if (indz > 0) {
				LdrowId[leng_Ld] = i - N_node_s;
				LdcolId[leng_Ld] = i;

			}
		}
	}

    /* Destructor */
    ~fdtdMesh(){
        // Set some important numbers to zero
        this->numCdtRow = (myint)0;
        this->leng_Vh = (myint)0;
        this->leng_S = (myint)0;

        // Free all the pointers
        free(this->xn);
        free(this->yn);
        free(this->zn);
        free(this->xnu);
        free(this->ynu);
        free(this->znu);
        //free(this->nodepos);
        free(this->Epoints);
        //free(this->edgelink);
        free(this->Hpoints);
        free(this->bd_node1);
        free(this->bd_node2);
        free(this->bd_edge);
        free(this->stackCdtMark);
        free(this->markEdge);
        free(this->markCell);
        free(this->cdtNumNode);
        free(this->sig);
        free(this->conductor);
        free(this->markNode);
        free(this->exciteCdtLayer);
        free(this->patch);
        free(this->bound);
        free(this->v0cRowId);
        free(this->v0cColId);
        free(this->v0cColIdo);
        free(this->v0cval);
        free(this->v0cvalo);
        free(this->v0caRowId);
        free(this->v0caColId);
        free(this->v0caColIdo);
        free(this->v0caval);
        free(this->v0cavalo);
        free(this->v0c2y0c2);
        free(this->v0c2y0c2o);
        free(this->yc);
        free(this->v0cy0c);
        free(this->AcRowId);
        free(this->AcRowId1);
        free(this->AcColId);
        free(this->Acval);
        free(this->AdRowId);
        free(this->AdRowId1);
        free(this->AdColId);
        free(this->Adval);
        free(this->crhs);
        free(this->v0d1RowId);
        free(this->v0d1ColId);
        free(this->v0d1ColIdo);
        free(this->v0d1val);
        free(this->v0d1valo);
        free(this->v0d1aRowId);
        free(this->v0d1aColId);
        free(this->v0d1aColIdo);
        free(this->v0d1aval);
        free(this->v0d1avalo);
        free(this->v0d2RowId);
        free(this->v0d2ColId);
        free(this->v0d2ColIdo);
        free(this->v0d2val);
        free(this->v0d2valo);
        free(this->v0d2aRowId);
        free(this->v0d2aColId);
        free(this->v0d2aColIdo);
        free(this->v0d2aval);
        free(this->v0d2avalo);
        free(this->yd);
        delete[] this->Vh;
        free(this->SRowId);
        free(this->SColId);
        free(this->Sval);
        free(this->y);
        free(this->J);
        free(this->v0csJ);
        free(this->Y);
    }
};




int meshAndMark(fdtdMesh* sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi, unordered_set<double> *portCoorx, unordered_set<double> *portCoory);
int compute_edgelink(fdtdMesh *sys, myint eno, myint &node1, myint &node2);
int parameterConstruction(fdtdMesh* sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi);
void freePara(fdtdMesh *sys);
int matrixConstruction(fdtdMesh *sys);
int portSet(fdtdMesh *sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi);
int mklMatrixMulti(fdtdMesh *sys, int &leng_A, int *aRowId, int *aColId, double *aval, int arow, int acol, int *bRowId, int *bColId, double *bval, int mark);
// The first is read row by row, and the second one is read column by column
int COO2CSR(vector<int>& rowId, vector<int>& ColId, vector<double>& val);
int mvMulti(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int>& bRowId, vector<int>& bColId, vector<double>& bval, double *index_val, int size);
int nodeAdd(int *index, int size, int total_size, fdtdMesh *sys, int &v0d2num, int &leng_v0d2, int mark);
int nodeAddLarger(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int *RowId, int *ColId, double *Val);
int nodeAdd_count(int *index, int size, int total_size, fdtdMesh *sys, int &v0d2num, int &leng_v0d2);
int nodeAddAvg(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int mark);
int nodeAddAvgLarger(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int* RowId, int *ColId, double *Val);
int nodeAddAvg_count(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng);
int interativeSolver(int N, int nrhs, double *rhs, int *ia, int *ja, double *a, int *ib, int *jb, double *b, double *solution, fdtdMesh *sys);
int output(fdtdMesh *sys);
int paraGenerator(fdtdMesh *sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi);
int yParaGenerator(fdtdMesh *sys);
int solveV0dSystem(fdtdMesh *sys, double *dRhs, double *y0d, int leng_v0d1);
int pardisoSolve(fdtdMesh *sys, double *rhs, double *solution, int leng_v0d1);
int COO2CSR_malloc(myint *rowId, myint *ColId, double *val, myint totalnum, myint leng, myint *rowId1);
int generateStiff(fdtdMesh *sys);
int merge_v0d1(fdtdMesh *sys, double block1_x, double block1_y, double block2_x, double block2_y, double block3_x, double block3_y, myint &v0d1num, myint &leng_v0d1, myint &v0d1anum, myint &leng_v0d1a, myint *map, double sideLen);
int merge_v0c(fdtdMesh *sys, double block_x, double block_y, double block2_x, double block2_y, myint &v0cnum, myint &leng_v0c, myint &v0canum, myint &leng_v0ca, myint *map);
int setsideLen(int node, double sideLen, int *markLayerNode, int *markProSide, fdtdMesh *sys);
int generateStiff(fdtdMesh *sys);
int mklMatrixMulti_nt(fdtdMesh *sys, myint &leng_A, myint *aRowId, myint *aColId, double *aval, myint arow, myint acol, myint *bRowId, myint *bColId, double *bval);
int find_Vh_central(fdtdMesh *sys, lapack_complex_double *u0, lapack_complex_double *u0a, int sourcePort);
int find_Vh_back(fdtdMesh* sys, int sourcePort);
int matrix_multi(char operation, lapack_complex_double *a, myint arow, myint acol, lapack_complex_double *b, myint brow, myint bcol, lapack_complex_double *tmp3);
int reference(fdtdMesh *sys, int freqNo, myint *RowId, myint *ColId, double *val);
int plotTime(fdtdMesh *sys, int sourcePort, double *u0d, double *u0c);
int avg_length(fdtdMesh *sys, int iz, int iy, int ix, double &lx, double &ly, double &lz);
myint generateLaplacian_count(fdtdMesh* sys);    // count how many nnz in L, and return the number
int generateLaplacian(fdtdMesh* sys, myint* rowId, myint* colId, double* val);    // generate the Laplacian matrix
#endif