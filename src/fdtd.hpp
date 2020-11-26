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
#include <map>
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
#include <mkl_spblas.h>
#include <mkl.h>


// Manipulate namespace
using namespace std;

// Fundamental physical constant macros
#define MU (4*M_PI*1.e-7)
#define CSPED (299792458.)
#define EPSILON0 (1./(CSPED*CSPED*MU))
#define MAXC (256*6)

// Solver control macros
#define SIGMA (5.8e+7)  // the load edge sigma //(3.508) // Default conductivity for conductors is copper (S/m)
#define DOUBLEMAX (1.e+30)
#define DOUBLEMIN (-1.e+30)
#define MINDISFRACX (0.001) // Fraction setting minimum discretization retained in x-directions after node merging in terms of smaller of x-extent
#define MINDISFRACY (0.001) // Fraction setting minimum discretization retained in y-directions after node merging in terms of smaller of y-extent
#define MINDISFRACZ (0.05) // Fraction setting minimum discretization retained in z-direction after node merging in terms of distance between closest layers
#define MAXDISFRACX (0.005) // Fraction setting largest discretization in x-direction in terms of x-extent
#define MAXDISFRACY (0.005) // Fraction setting largest discretization in y-direction in terms of y-extent
#define MAXDISLAYERZ (2.)// Largest discretization in z-direction represented as fewest nodes placed between closest layers (1. = distance between closest layers, 2. = half distance between closest layers)
#define DT (1.e-11) // Time step for finding high-frequency modes (s)

// Debug testing macros (comment out if not necessary)
#define UPPER_BOUNDARY_PEC
//#define LOWER_BOUNDARY_PEC
//#define CONDUCTOR_PEC    // all conductor edges are removed
#define PRINT_NODE_COORD
#define PRINT_DIS_COUNT (1)
//#define SKIP_MARK_CELL
#define PRINT_VERBOSE_TIMING // Terminal output has extra runtime clock information
//#define PRINT_PORT_SET
//#define PRINT_V0D_BLOCKS
#define PRINT_V0_Z_PARAM
//#define PRINT_V0_Vh_Z_PARAM
#define SKIP_PARDISO // Remove PARDISO solver code
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
	double* x;
	double* y;
	double* z;
	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;
	double sigma;
	myint* cdtInNode;
	int* markCdtInNode;
	int layer;
};

class fdtdCdt {
public:
	myint* node;
	int markPort;
	int* portNode;
	int portind;
	myint cdtNodeind;
};

class fdtdBound {
public:
	int numBound;
	int numPeri; /*number of open surrounding boundaries */
	int topNum; /*top open boundary index*/
	int botNum; /*bot open boundary index*/
	int* p1;
	int* p2;
	int* bound;
	int* ukwnBegEdge; /*starting edge of each bound*/
	int* ukwnBegNode; /*starting node of each bound*/
	int* ukwnNumEdge; /*number of edges of each bound*/
	int* ukwnNumNode;  /*number of nodes of each bound*/
	int top;
	int bot;
	double** norm;/*bound norm direction*/
	int* numPort;
	int** port;
};

class fdtdPort {
	/* Port information */
public:
	double* x, *y, *z;
	int multiplicity;                        // Number of port sides to excite
	vector<double> x1, x2, y1, y2, z1, z2;   // Supply/return coordinates for each port side
	vector<myint> portCnd;                   // Conductor index of supply point for each port side
	vector<vector<myint>> portEdge;          // Edge indices from return to supply with nonzero current density for each port side
	vector<double> portArea;                 // Area of each port side (m^2)
	vector<int> portDirection;               // Sign of current density component for each port side
	int* node;

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
	int* node;
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

	double* xn, *yn, *zn;    // coordinates of the nodes along x, y, z
	double* xnu, *ynu, *znu;

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
	double* Epoints;
	//myint *edgelink;
	double* Hpoints;
	vector<vector<pair<myint, double>>> nodeEdge;    // for each node which edge is connected with this node
	vector<vector<pair<myint, double>>> nodeEdgea;    // for each node which edge is connected with this node

													  /* Upper and lower PEC */
	int* bd_node1;   //lower PEC
	int* bd_node2;   //upper PEC
	int* bd_edge;
	set<myint> ubde, lbde;    // upper boundary edge and lower boundary edge
	set<myint> ubdn, lbdn;    // upper boundary node and lower boundary node
	myint* mapEdge;   // map the original edges to the new edge # with upper and lower PEC boundaries
	myint* mapEdgeR;    // map the new edge # to the original edges
	myint* mapNode;
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
	double* stackCdtMark;

	/* Conductor parameters */
	vector<fdtdOneCondct> conductorIn;
	myint numCdtRow;                      // how many input rows
	myint numCdt;                         // number of isolated conductors in design
	myint* markEdge;                      // mark if this edge is inside a conductor
	myint* markCell;
	myint* cdtNumNode;                    // number of nodes along each isolated conductor
	double* sig;
	fdtdCdt* conductor;                   // information about isolated conductors
	myint* markNode;                      // mark this node if it is inside the conductor
	vector<vector<int>> edgeCell;         // for each cell which edge is around it
	vector<vector<double>> edgeCellArea;  // for each cell the area of the perpendicular rectangle
	vector<int> acu_cnno;                 // accumulated conductor number of nodes
	vector<int> cindex;
	int* exciteCdtLayer;
	unordered_set<int> cond2condIn;       // put the active conductors' corresponding conductorIn #
	vector<bool> markProSide;             // mark nodes near excited conductors for less aggressive node merging in rapid field change area
	int* markCdt;                         // mark 1 if this conductor cannot be formed to V0g by V0b

										  /* Patch information */
	fdtdPatch* patch;

	/* Boundary information */
	fdtdBound* bound;

	/* V0c row, column, value */
	myint* v0cRowId;
	myint* v0cColId;
	myint* v0cColIdo;
	double* v0cval;
	double* v0cvalo;

	myint* v0caRowId;
	myint* v0caColId;
	myint* v0caColIdo;
	double* v0caval;
	double* v0cavalo;
	myint leng_v0c, leng_v0ca, v0cnum, v0canum;
	myint* mapc;

	/* norm of V0d, V0da, V0c, V0ca */
	double* v0dn, *v0dan, *v0cn, *v0can;

	/* V0c2 and yc row, column, value */
	double* v0c2y0c2;
	double* v0c2y0c2o;
	double* yc;
	double* v0cy0c;

	/* V0c'*D_sig*V0c row, column, value */
	myint* AcRowId;
	myint* AcRowId1;
	myint* AcColId;
	double* Acval;
	myint leng_Ac;
	myint* AdRowId;
	myint* AdRowId1;
	myint* AdColId;
	double* Adval;
	myint leng_Ad;

	myint* MdRowId, *MdColId;
	double* Mdval;
	myint leng_Md;

	double* crhs;

	/* V0d row, column, value */
	myint* v0d1RowId;
	myint* v0d1ColId;
	myint* v0d1ColIdo;
	double* v0d1val;
	double* v0d1valo;

	myint* v0d1aRowId;
	myint* v0d1aColId;
	myint* v0d1aColIdo;
	double* v0d1aval;
	double* v0d1aval1;
	double* v0d1avalo;

	myint* v0d2RowId;
	myint* v0d2ColId;
	myint* v0d2ColIdo;
	double* v0d2val;
	double* v0d2valo;

	myint* v0d2aRowId;
	myint* v0d2aColId;
	myint* v0d2aColIdo;
	double* v0d2aval;
	double* v0d2avalo;
	myint leng_v0d1, leng_v0d1a, v0d1num, v0d1anum, v0dnum;
	myint leng_v0d, leng_v0dd;   // leng_v0d is the column number for V0d in V0db, leng_v0dd is the column number for V0d in V0d1
	myint* v0dbRowId, *v0dbColId, *v0dbColIdo;
	double* v0dbval, *v0dbaval;
	myint leng_v0db, v0dbnum, leng_v0db1;
	myint* mapd;

	myint* XgRowId, *XgColId;
	double* Xgval;
	myint leng_Xg;

	double* yd;

	/* Vh */
	lapack_complex_double* Vh;
	myint leng_Vh;

	/* Se and Sh */
	myint* SRowId;
	myint* SColId;
	double* Sval;
	myint leng_S;

	/* Laplacian decomposed v0 xh left matrix */
	myint* LlRowId;
	myint* LlRowIdo;
	myint* LlColId;
	double* Llval;
	myint leng_Ll;

	/* outside conductor Laplacian rowId, colId, val */
	myint* LooRowId, *LooRowId1, *LooColId;
	double* Looval;
	myint leng_Loo;

	/* map the edges to the groups inside or outside */
	myint* mapio, *mapioR;
	myint inside, outside;    // size ranges of inside and outside


							  /* Solution storage */
	complex<double>* y;
	vector<complex<double>> x;    // the solution involving all the sourcePorts

								  /* Port information */
	int numPorts;
	vector<fdtdPort> portCoor;

	/* Current sources */
	double* J;

	/* Current V0c,s^T*I matrix */
	complex<double>* v0csJ;
	complex<double>* Y;

	/* Default Constructor */
	fdtdMesh() {
		// Set some important numbers to zero
		this->outedge = (myint)0;
		this->inedge = (myint)0;
		this->numCdtRow = (myint)0;
		this->leng_Vh = (myint)0;
		this->leng_S = (myint)0;
		inside = 0;
		outside = 0;
		leng_Loo = 0;

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
		this->v0d1aval1 = NULL;
		this->v0d1avalo = NULL;
		this->v0dbRowId = NULL;
		this->v0dbColId = NULL;
		this->v0dbColIdo = NULL;
		this->v0dbval = NULL;
		this->v0dbaval = NULL;
		this->XgRowId = NULL;
		this->XgColId = NULL;
		this->Xgval = NULL;
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
		this->v0dn = NULL;
		this->v0dan = NULL;
		this->v0cn = NULL;
		this->v0can = NULL;
		this->yd = NULL;
		this->Vh = NULL;
		this->SRowId = NULL;
		this->SColId = NULL;
		this->Sval = NULL;
		this->y = NULL;
		this->J = NULL;
		this->v0csJ = NULL;
		this->Y = NULL;
		this->mapd = NULL;
		this->mapc = NULL;
		this->LlRowId = NULL;
		this->LlColId = NULL;
		this->Llval = NULL;
		this->mapio = NULL;
		this->mapioR = NULL;
		this->LooRowId = NULL;
		this->LooRowId1 = NULL;
		this->LooColId = NULL;
		this->Looval = NULL;
		this->mapEdge = NULL;
		this->mapEdgeR = NULL;
		this->markCdt = NULL;


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

	/* Print dx min */
	void printDxMin() {
		double res = DOUBLEMAX;

		for (int i = 1; i < nx; ++i) {
			if (xn[i] - xn[i - 1] < res) {
				res = xn[i] - xn[i - 1];
			}
		}
		cout << "dx min is " << res << endl;
	}

	/* Print dy min */
	void printDyMin() {
		double res = DOUBLEMAX;

		for (int i = 1; i < ny; ++i) {
			if (yn[i] - yn[i - 1] < res) {
				res = yn[i] - yn[i - 1];
			}
		}
		cout << "dy min is " << res << endl;
	}

	/* Print dz min */
	void printDzMin() {
		double res = DOUBLEMAX;

		for (int i = 1; i < nz; ++i) {
			if (zn[i] - zn[i - 1] < res) {
				res = zn[i] - zn[i - 1];
			}
		}
		cout << "dz min is " << res << endl;
	}

	/* Print conductorIn */
	void printConductorIn() {
		int i, j;

		cout << "Print conductorIn information: " << endl;
		for (i = 0; i < this->numCdtRow; i++) {
			for (j = 0; j < this->conductorIn[i].numVert - 1; j++) {
				if (this->conductorIn[i].layer == 5) {
					cout << this->conductorIn[i].x[j] << " " << this->conductorIn[i].y[j] << " " << this->conductorIn[i].zmin << " " << this->conductorIn[i].zmax << endl;
				}
			}
		}
	}

	void setCurrent(int sourcePort, double* J) {
		for (int sourcePortSide = 0; sourcePortSide < this->portCoor[sourcePort].multiplicity; sourcePortSide++) {
			for (int indEdge = 0; indEdge < this->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
				J[this->mapEdge[this->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]]] = this->portCoor[sourcePort].portDirection[sourcePortSide];
			}
		}
	}

	/* Check whether the point's markNode */
	void checkPoint(double x, double y, double z, unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi) {
		int inx, iny, inz;

		inx = xi[x];
		iny = yi[y];
		inz = zi[z];
		cout << "inx " << inx << " iny " << iny << " inz " << inz << endl;
		cout << "x " << this->xn[inx] << " y " << this->yn[iny] << " z " << this->zn[inz] << endl;
		cout << "This node's mark is " << this->markNode[inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny] << endl;
	}

	/* Get the eps for a particular edge */
	double getEps(myint eno) {
		/* eno : edge number from the original mesh */
		if (eno >= N_edge_s && eno < N_edge - N_edge_s && eno % (N_edge_s + N_edge_v) < N_edge_s) {   // this edge is not on the top and bottom planes and this edge are on the interface between two layers
			return (stackEpsn[(eno + N_edge_v) / (N_edge_s + N_edge_v)] + stackEpsn[(eno + N_edge_s + N_edge_v * 2) / (N_edge_s + N_edge_v)]) / 2 * EPSILON0;
		}
		else {
			return stackEpsn[(eno + N_edge_v) / (N_edge_s + N_edge_v)] * EPSILON0;
		}
	}

	/* Is point (x,y) within the polygon? */
	bool polyIn(double x, double y, int inPoly) {
		int npol;
		myint indi = 0, indj = 0;
		bool isCond = false;
		double disMin = 1.e-10;

		npol = this->conductorIn[inPoly].numVert;

		for (indi = 0, indj = npol - 1; indi < npol; indj = indi++) {
			//if ((abs(y - this->conductorIn[inPoly].y[indj]) < disMin && abs(y - this->conductorIn[inPoly].y[indi]) < disMin &&
			//	((x >= this->conductorIn[inPoly].x[indj] && x <= this->conductorIn[inPoly].x[indi]) ||
			//	(x >= this->conductorIn[inPoly].x[indi] && x <= this->conductorIn[inPoly].x[indj])))) {    // on x direction edge
			//	return true;
			//}
			//else if (abs(x - this->conductorIn[inPoly].x[indj]) < disMin && abs(x - this->conductorIn[inPoly].x[indi]) < disMin &&
			//	((y >= this->conductorIn[inPoly].y[indj] && y <= this->conductorIn[inPoly].y[indi]) ||
			//	(y >= this->conductorIn[inPoly].y[indi] && y <= this->conductorIn[inPoly].y[indj]))) {    // on y direction edge
			//	return true;
			//}

			//else if ((abs(this->conductorIn[inPoly].y[indi] - this->conductorIn[inPoly].y[indj]) > disMin &&
			//	(((this->conductorIn[inPoly].y[indi] <= y) && (y < this->conductorIn[inPoly].y[indj])) ||
			//	((this->conductorIn[inPoly].y[indj] <= y) && (y < this->conductorIn[inPoly].y[indi])))) &&
			//		(x < (this->conductorIn[inPoly].x[indj] - this->conductorIn[inPoly].x[indi]) * (y - this->conductorIn[inPoly].y[indi]) /
			//	(this->conductorIn[inPoly].y[indj] - this->conductorIn[inPoly].y[indi]) + this->conductorIn[inPoly].x[indi])) {
			//	isCond = !isCond;
			//}
			if ((((this->conductorIn[inPoly].y[indi] <= y) && (y < this->conductorIn[inPoly].y[indj])) ||
				((this->conductorIn[inPoly].y[indj] <= y) && (y < this->conductorIn[inPoly].y[indi]))) &&
				(x < (this->conductorIn[inPoly].x[indj] - this->conductorIn[inPoly].x[indi]) * (y - this->conductorIn[inPoly].y[indi]) /
				(this->conductorIn[inPoly].y[indj] - this->conductorIn[inPoly].y[indi]) + this->conductorIn[inPoly].x[indi])) {
				isCond = !isCond;
			}

		}
		return isCond;
	}

	/* Map the after removed PEC edges to the inside outside edges */
	void mapEdgeInsideOutside() {
		myint edge = N_edge - bden;
		mapio = (myint*)malloc(edge * sizeof(myint));
		mapioR = (myint*)malloc(edge * sizeof(myint));
		inside = 0;
		outside = 0;

		/* count how many edges are outside the conductors */
		for (int ind = 0; ind < edge; ++ind) {
			if (markEdge[mapEdgeR[ind]] == 0)
			{
				inside++;
			}
		}
		/* map the edges to inside and outside edges */
		for (int ind = 0; ind < edge; ++ind) {
			if (markEdge[mapEdgeR[ind]] == 0) {
				mapio[ind] = outside;
				mapioR[outside] = ind;
				outside++;
			}
			else {
				mapio[ind] = inside;
				mapioR[inside] = ind;
				inside++;
			}
		}
		return;
	}

	//dj new algorithms for finding conducting edges and nodes

	void findInsideCondNew(unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi) {

		double xmin;
		double xmax;
		double ymin;
		double ymax;
		double zmin, zmax;
		ofstream out;
		ofstream out2;
		int indi, indk, ip = 0, iq = 0, nip = 2649;
		myint node, node1, node2, idi, idj, idk;
		myint i1, i2, j1, j2, k1, k2, eno;

		//cout << "Pls enter a conductor " << endl;
		//cin >> nip;

		for (indi = 0; indi < this->numCdtRow; indi++) {
			xmin = this->xlim2*this->lengthUnit;
			xmax = this->xlim1*this->lengthUnit;
			ymin = this->ylim2*this->lengthUnit;
			ymax = this->ylim1*this->lengthUnit;
			//cout << "xmin max, ymin, ymax" << xmin << " " << xmax <<" "<<ymin<<" "<<ymax<< endl;

			//if (this->conductorIn[indi].zmax == this->conductorIn[indi].zmin)
			//continue;
			//if (this->conductorIn[indi].zmin >= this->zlim1*this->lengthUnit && this->conductorIn[indi].zmin <= this->zlim2*this->lengthUnit) {
			//if (this->conductorIn[indi].zmax >= this->zlim1*this->lengthUnit && this->conductorIn[indi].zmax <= this->zlim2*this->lengthUnit) {
			if (this->conductorIn[indi].zmax != this->conductorIn[indi].zmin) {

				ip++;
				//if (indi == nip) {
				//	out.open("Cnt1.txt", std::ofstream::out | std::ofstream::trunc);
				//	out.precision(17);
				//	out2.open("Cnt1_new.txt", std::ofstream::out | std::ofstream::trunc);
				//	out2.precision(17);
				//}
				for (indk = 0; indk < this->conductorIn[indi].numVert; indk++) {
					//if (this->conductorIn[indi].x[indk] < this->xlim1 * lengthUnit || this->conductorIn[indi].x[indk] >= this->xlim2 * lengthUnit || this->conductorIn[indi].y[indk] < this->ylim1 * lengthUnit || this->conductorIn[indi].y[indk] >= this->ylim2 * lengthUnit)
					//	continue;
					xmin = fmin(xmin, this->conductorIn[indi].x[indk]);
					xmax = fmax(xmax, this->conductorIn[indi].x[indk]);
					ymin = fmin(ymin, this->conductorIn[indi].y[indk]);
					ymax = fmax(ymax, this->conductorIn[indi].y[indk]);
					//cout << "indi " << indi << "indk " << indk << "" << this->conductorIn[indi].x[indk] << " " << this->conductorIn[indi].y[indk] << endl;
					//if (indi == nip) {
					//	out << this->conductorIn[indi].x[indk] << " " << this->conductorIn[indi].y[indk] << endl;
					//}
				}
				//if (indi == nip) {
				//	out.close();
				//}
				//cout << "new xmin max, ymin, ymax" << xmin << " " << xmax << " " << ymin << " " << ymax << endl;
				zmin = this->conductorIn[indi].zmin;
				zmax = this->conductorIn[indi].zmax;
				//cout << "zmin, zmax" << zmin << " " << zmax << endl;
				i1 = xi[xmin];//get to know the bounding box;
				i2 = xi[xmax];
				j1 = yi[ymin];
				j2 = yi[ymax];
				k1 = zi[zmin];
				k2 = zi[zmax];
				//cout << "i1..." << i1 << " " << i2 << " " << j1 << " " << j2 << " " << k1 << " " << k2 << endl;

				//cin >> iq;
				for (idi = i1; idi <= i2; idi++) {
					for (idj = j1; idj <= j2; idj++) {
						if (polyIn(this->xn[idi], this->yn[idj], indi)) {
							//if (indi == nip) {
							//	out2 << this->xn[idi] << " " << this->yn[idj] << " " << k1 << " " << k2 << endl;
							//}
							//cout << idi << " " << idj << endl;
							for (idk = k1; idk <= k2; idk++) {
								node = idk * this->N_node_s + (this->N_cell_y + 1) * idi + idj;
								this->markNode[node] = indi + 1;
							}
						}
					}
				}
				//if (indi == nip) {
				//	out2.close();
				//	//cin >> iq;
				//}
				//dj test
				for (idi = i1; idi < i2; idi++) {
					for (idj = j1; idj <= j2; idj++) {
						for (idk = k1; idk <= k2; idk++) {// x-direction
							eno = idk * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + idi * (this->N_cell_y + 1) + idj;
							node1 = idk * this->N_node_s + (this->N_cell_y + 1) * idi + idj;
							node2 = idk * this->N_node_s + (this->N_cell_y + 1) * (idi + 1) + idj;
							if (this->markNode[node1] == (indi + 1) && this->markNode[node2] == (indi + 1)) {
								this->markEdge[eno] = indi + 1;
							}
						}
					}
				}
				for (idi = i1; idi <= i2; idi++) {
					for (idj = j1; idj < j2; idj++) {
						for (idk = k1; idk <= k2; idk++) {// y-direction
							eno = idk * (this->N_edge_s + this->N_edge_v) + idi * (this->N_cell_y) + idj;
							node1 = idk * this->N_node_s + (this->N_cell_y + 1) * idi + idj;
							node2 = idk * this->N_node_s + (this->N_cell_y + 1) * idi + idj + 1;
							if (this->markNode[node1] == (indi + 1) && this->markNode[node2] == (indi + 1)) {
								this->markEdge[eno] = indi + 1;
							}
						}
					}
				}
				for (idi = i1; idi <= i2; idi++) {
					for (idj = j1; idj <= j2; idj++) {
						for (idk = k1; idk < k2; idk++) {// z-direction
							eno = idk*(this->N_edge_s + this->N_edge_v) + this->N_edge_s + idi * (this->N_cell_y + 1) + idj;
							node1 = idk * this->N_node_s + (this->N_cell_y + 1) * idi + idj;
							node2 = (idk + 1) * this->N_node_s + (this->N_cell_y + 1) * idi + idj;
							if (this->markNode[node1] == (indi + 1) && this->markNode[node2] == (indi + 1)) {
								this->markEdge[eno] = indi + 1;
							}
						}
					}
				}
				//dj


			}
			//}
			//}
		}
		cout << "old conduct number " << this->numCdtRow << "new number " << ip << endl;
		//cin >> indi;
		// the following part cannot be used because the node conductor can get changed during the assignment
		// So I rewrote the above x-y-z section to get markEdge.
		/*for (indi = 0; indi < this->N_edge; indi++) {
		compute_edgelink(indi, node1, node2);
		if (this->markNode[node1] == this->markNode[node2] && this->markNode[node1] > 0) {
		this->markEdge[indi] = this->markNode[node1];
		}
		}*/
	}


	/* Find nodes inside conductors by judging whether the point in each small window is inside the polygon or not */
	void findInsideCond(unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi) {
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
			if (!(conductorIn[indi].zmin >= zlim1 * lengthUnit && conductorIn[indi].zmin <= zlim2 * lengthUnit && conductorIn[indi].zmax >= zlim1 * lengthUnit && conductorIn[indi].zmax <= zlim2 * lengthUnit))
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
						if (xcoory[ss][0] > min(y1, y2)) {
							xcoory[ss][0] = min(y1, y2);
						}
						if (xcoory[ss][1] < max(y1, y2)) {
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
						if (xcoory[ss][0] > min(y1, y2)) {
							xcoory[ss][0] = min(y1, y2);
						}
						if (xcoory[ss][1] < max(y1, y2)) {
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
						if (xcoory[ss][0] > min(y1, y2)) {
							xcoory[ss][0] = min(y1, y2);
						}
						if (xcoory[ss][1] < max(y1, y2)) {
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
						if (xcoory[ss][0] > min(y1, y2)) {
							xcoory[ss][0] = min(y1, y2);
						}
						if (xcoory[ss][1] < max(y1, y2)) {
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
				for (auto xcooryi : xcoory) {
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
	void findInsideCond_xrangey(unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi) {
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
			if (!(conductorIn[indi].zmin >= zlim1 * lengthUnit && conductorIn[indi].zmin <= zlim2 * lengthUnit && conductorIn[indi].zmax >= zlim1 * lengthUnit && conductorIn[indi].zmax <= zlim2 * lengthUnit))
				continue;

			xrange.clear();
			xcoor.clear();
			xcoorv.clear();
			xcoory.clear();
			mark1 = 0;
			//cout << endl;
			//cout << endl;
			//cout << this->conductorIn[indi].numVert << endl;
			for (indj = 0; indj < this->conductorIn[indi].numVert - 1; indj++) {
				//cout << this->conductorIn[indi].x[indj] << " " << this->conductorIn[indi].y[indj] << endl;
				if ((this->conductorIn[indi].x[indj] >= xlim1 && this->conductorIn[indi].x[indj] <= xlim2) ||
					(this->conductorIn[indi].x[indj + 1] >= xlim1 && this->conductorIn[indi].x[indj + 1] <= xlim2)) {
					// if one of the nodes of the range is inside the total structure
					if (this->conductorIn[indi].x[indj] < xlim1 || this->conductorIn[indi].x[indj] > xlim2) {
						if (this->conductorIn[indi].x[indj] < xlim1) {
							xi[this->conductorIn[indi].x[indj]] = 0;
						}
						else if (this->conductorIn[indi].x[indj] > xlim2) {
							xi[this->conductorIn[indi].x[indj]] = nx - 1;
						}
					}
					else if (this->conductorIn[indi].x[indj + 1] < xlim1 || this->conductorIn[indi].x[indj + 1] > xlim2) {
						if (this->conductorIn[indi].x[indj + 1] < xlim1) {
							xi[this->conductorIn[indi].x[indj + 1]] = 0;
						}
						else if (this->conductorIn[indi].x[indj + 1] > xlim2) {
							xi[this->conductorIn[indi].x[indj + 1]] = nx - 1;
						}
					}



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
			}
			//cout << this->conductorIn[indi].x[indj] << " " << this->conductorIn[indi].y[indj] << endl;
			if ((this->conductorIn[indi].x[indj] >= xlim1 && this->conductorIn[indi].x[indj] <= xlim2) ||
				(this->conductorIn[indi].x[0] >= xlim1 && this->conductorIn[indi].x[0] <= xlim2))
				// if the x range is outside the total size
			{

				if (this->conductorIn[indi].x[indj] < xlim1 || this->conductorIn[indi].x[indj] > xlim2) {
					if (this->conductorIn[indi].x[indj] < xlim1) {
						xi[this->conductorIn[indi].x[indj]] = 0;
					}
					else if (this->conductorIn[indi].x[indj] > xlim2) {
						xi[this->conductorIn[indi].x[indj]] = nx - 1;
					}
				}
				else if (this->conductorIn[indi].x[0] < xlim1 || this->conductorIn[indi].x[0] > xlim2) {
					if (this->conductorIn[indi].x[0] < xlim1) {
						xi[this->conductorIn[indi].x[0]] = 0;
					}
					else if (this->conductorIn[indi].x[0] > xlim2) {
						xi[this->conductorIn[indi].x[0]] = nx - 1;
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
			}
			if (xcoor.size() == 0) {
				continue;
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

				if ((this->conductorIn[indi].x[indj] >= xlim1 && this->conductorIn[indi].x[indj] <= xlim2) ||
					(this->conductorIn[indi].x[indj + 1] >= xlim1 && this->conductorIn[indi].x[indj + 1] <= xlim2)) {
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
			}


			if ((this->conductorIn[indi].x[indj] >= xlim1 && this->conductorIn[indi].x[indj] <= xlim2) ||
				(this->conductorIn[indi].x[0] >= xlim1 && this->conductorIn[indi].x[0] <= xlim2))
				// if the x range is outside the total size
			{
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

						if ((y1 == y2 && y1 == 0))   // if y are all on the boundary = 0 we neglect because they may be outside the total structure
							continue;
						else if (y1 > y2 && y2 == 0) {
							y2 = ny - 1;
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
						if (y1 == y2 && y1 == 0)    // if y is on the boundary = 0 this may be they are all outside the total structure
							continue;
						else if (y1 > y2 && y2 == 0) {   // if y2 is outside then put y2 as the max position
							y2 = ny - 1;
						}

						for (indj = x1; indj <= x2; indj++) {
							for (indl = zi[this->conductorIn[indi].zmin]; indl <= zi[this->conductorIn[indi].zmax]; indl++) {
								for (indk = y1; indk < y2; indk++) {
									if (indj < x2) {
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
								if (indj < x2) {
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
	void findCond2CondIn(int inx, int iny, int inz) {
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
			eno = inz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (inx - 1) * (this->N_cell_y + 1) + iny;
			if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
				this->cond2condIn.insert(this->markEdge[eno]);
			}
		}
		if (inx != this->nx - 1) {    // this node is not on the right plane
			eno = inz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + inx * (this->N_cell_y + 1) + iny;
			if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
				this->cond2condIn.insert(this->markEdge[eno]);
			}
		}
		if (iny != 0) {    // this node is not on the front plane
			eno = inz * (this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny - 1;
			if (this->markEdge[eno] != 0 && this->cond2condIn.find(this->markEdge[eno]) == this->cond2condIn.end()) {
				this->cond2condIn.insert(this->markEdge[eno]);
			}
		}
		if (iny != this->ny - 1) {    // this node is not on the back plane
			eno = inz * (this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny;
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
				eno = inz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (inx - 1) * (this->N_cell_y + 1) + iny;   // left edge
				if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
					this->lbde.insert(eno);
				}
			}
			if (inx != this->nx - 1) {    // this node is not on the right plane
				eno = inz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + inx * (this->N_cell_y + 1) + iny;   // right edge
				if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
					this->lbde.insert(eno);
				}
			}
			if (iny != 0) {    // this node is not on the front plane
				eno = inz * (this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny - 1;   // front edge
				if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
					this->lbde.insert(eno);
				}
			}
			if (iny != this->ny - 1) {    // this node is not on the back plane
				eno = inz * (this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny;   // back edge
				if (this->markEdge[eno] != 0 && this->lbde.find(eno) == this->lbde.end()) {
					this->lbde.insert(eno);
				}
			}
		}
		else if (inz == this->nz - 1) {
			this->ubdn.insert(inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny);   // this node is upper PEC boundary
			if (inx != 0) {    // this node is not on the left plane
				eno = inz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (inx - 1) * (this->N_cell_y + 1) + iny;    // left edge
				if (this->markEdge[eno] != 0 && this->ubde.find(eno) == this->ubde.end()) {
					this->ubde.insert(eno);
				}
			}
			if (inx != this->nx - 1) {    // this node is not on the right plane
				eno = inz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + inx * (this->N_cell_y + 1) + iny;    // right edge
				if (this->markEdge[eno] != 0 && this->ubde.find(eno) == this->ubde.end()) {
					this->ubde.insert(eno);
				}
			}
			if (iny != 0) {    // this node is not on the front plane
				eno = inz * (this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny - 1;    // front edge
				if (this->markEdge[eno] != 0 && this->ubde.find(eno) == this->ubde.end()) {
					this->ubde.insert(eno);
				}
			}
			if (iny != this->ny - 1) {    // this node is not on the back plane
				eno = inz * (this->N_edge_s + this->N_edge_v) + inx * this->N_cell_y + iny;    // back edge
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
		this->mapEdge = (myint*)calloc(this->N_edge, sizeof(myint));
		this->mapEdgeR = (myint*)calloc(this->N_edge - this->bden, sizeof(myint));
		int ind, count = 0;    // count : how many boundary edges already identified
#ifdef CONDUCTOR_PEC
		for (int ind = 0; ind < this->N_edge; ++ind) {
			if (this->markEdge[ind]) {
				this->mapEdge[ind] = -1;
				count++;
			}
			else {
				this->mapEdge[ind] = ind - count;
				this->mapEdgeR[ind - count] = ind;
			}
		}
#else
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
		//for (ind = 0; ind < N_edge; ++ind) {
		//	cout << mapEdge[ind] << " ";
		//}
#endif
	}

	void setMapNode() {
		/* map the node number before PEC to after PEC */

		mapNode = (myint*)calloc(N_node, sizeof(myint));
		int ind, count = -1;

		for (ind = 0; ind < N_node; ++ind) {
			if (ind >= 0 && ind < N_node_s && lbdn.find(ind) != lbdn.end()) {
				continue;
			}
			else if (ind >= N_node - N_node_s && ind < N_node && ubdn.find(ind) != ubdn.end()) {
				continue;
			}
			else {
				count++;
				mapNode[ind] = count;
			}
		}

		if (ubdn.size() > 0) {    // first upper and then lower because lower is removed to make the indexes continuous
			count++;
		}
		for (auto ubdni : ubdn) {

			mapNode[ubdni] = count;    // all the nodes in the upper boundary have 1 common map value
		}

		if (lbdn.size() > 0) {
			count++;
		}
		for (auto lbdni : lbdn) {

			mapNode[lbdni] = count;    // all the nodes in the lower boundary have 1 common map value
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
	void merge_v0d1(double block1_x, double block1_y, double block2_x, double block2_y, double block3_x, double block3_y, double sideLen) {
		int* visited;
		clock_t t1;
		double t, ta;
		double ratio;
		double startx, starty;    // the start coordinates of each block
		queue<int> st;    // dfs stack
		int indsize;
		myint indi = 0;
		int indx, indy;
		int mark;
		int indnum;
		int* markLayerNode = (int*)calloc(this->N_node_s, sizeof(int));
		/* Mark layer nodes from port sides */
		for (int indPort = 0; indPort < this->numPorts; indPort++) {

			for (int indPortSide = 0; indPortSide < this->portCoor[indPort].portCnd.size(); indPortSide++) {
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
								startx = this->xn[ix];
								starty = this->yn[iy];
								node_group.push_back(base);
								nodegs = node_group.size() - 1;
								node_group[nodegs].insert(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);
								mapd[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
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
												st.push((indx + 1) * (this->N_cell_y + 1) + indy);
												visited[(indx + 1) * (this->N_cell_y + 1) + indy] = 1;
												mapd[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
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

												st.push((indx - 1) * (this->N_cell_y + 1) + indy);
												visited[(indx - 1) * (this->N_cell_y + 1) + indy] = 1;
												mapd[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = count;
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

												st.push((indx) * (this->N_cell_y + 1) + indy + 1);
												visited[(indx) * (this->N_cell_y + 1) + indy + 1] = 1;
												mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
											}
										}
									}
									if (indy != 0) {    // it must have a closer y edge, thus closer y node
										if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block1_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block1_y) {    // this node is within the block area
											if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
												&& visited[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
												&& markLayerNode[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
												&& !this->markProSide[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1]
												&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

												st.push((indx) * (this->N_cell_y + 1) + indy - 1);
												visited[(indx) * (this->N_cell_y + 1) + indy - 1] = 1;
												mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
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
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (iz != this->nz - 1) {   // this node is not on the top plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indx != 0) {    // this node is not on the left plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}

									}
									if (indx != this->nx - 1) {    // this node is not on the right plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indy != 0) {    // this node is not on the front plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indy != this->ny - 1) {   // this node is not on the back plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
								}
								leng_v0d1++;

								count++;
							}

							else if (markLayerNode[ix * (this->N_cell_y + 1) + iy] == 1 && !this->markProSide[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {//&& this->exciteCdtLayer[iz] == 1) {    // this point is not visited and it is outside the conductor, in the projection of the excited conductor
								startx = this->xn[ix];
								starty = this->yn[iy];

								mapd[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
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

												st.push((indx + 1) * (this->N_cell_y + 1) + indy);
												visited[(indx + 1) * (this->N_cell_y + 1) + indy] = 1;
												mapd[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
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

												st.push((indx - 1) * (this->N_cell_y + 1) + indy);
												visited[(indx - 1) * (this->N_cell_y + 1) + indy] = 1;
												mapd[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
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

												st.push((indx) * (this->N_cell_y + 1) + indy + 1);
												visited[(indx) * (this->N_cell_y + 1) + indy + 1] = 1;
												mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
											}
										}
									}
									if (indy != 0) {    // it must have a closer y edge, thus closer y node
										if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block2_y) {    // this node is within the block area
											if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
												&& visited[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
												&& markLayerNode[(indx) * (this->N_cell_y + 1) + indy - 1] == 1
												&& !this->markProSide[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1]
												&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

												st.push((indx) * (this->N_cell_y + 1) + indy - 1);
												visited[(indx) * (this->N_cell_y + 1) + indy - 1] = 1;
												mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
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
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (iz != this->nz - 1) {   // this node is not on the top plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indx != 0) {    // this node is not on the left plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indx != this->nx - 1) {    // this node is not on the right plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indy != 0) {    // this node is not on the front plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indy != this->ny - 1) {   // this node is not on the back plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
								}
								leng_v0d1++;

								count++;

							}

							else {
								startx = this->xn[ix];
								starty = this->yn[iy];

								mapd[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
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

												st.push((indx + 1) * (this->N_cell_y + 1) + indy);
												visited[(indx + 1) * (this->N_cell_y + 1) + indy] = 1;
												mapd[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);
											}
										}
									}
									if (indx != 0) {    // it must have a left x edge, thus left x node
										if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block3_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block3_y) {    // this node is within the block area
											if (this->markNode[iz * this->N_node_s + st.front() - this->N_cell_y - 1] == 0
												&& visited[(indx - 1) * (this->N_cell_y + 1) + indy] == 0
												&& this->markProSide[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy]
												&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

												st.push((indx - 1) * (this->N_cell_y + 1) + indy);
												visited[(indx - 1) * (this->N_cell_y + 1) + indy] = 1;
												mapd[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
											}
										}
									}
									if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
										if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block3_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block3_y) {    // this node is within the block area
											if (this->markNode[iz * this->N_node_s + st.front() + 1] == 0
												&& visited[indx * (this->N_cell_y + 1) + indy + 1] == 0
												&& this->markProSide[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1]
												&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

												st.push((indx) * (this->N_cell_y + 1) + indy + 1);
												visited[(indx) * (this->N_cell_y + 1) + indy + 1] = 1;
												mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
											}
										}
									}
									if (indy != 0) {    // it must have a closer y edge, thus closer y node
										if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block3_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block3_y) {    // this node is within the block area
											if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
												&& visited[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
												&& this->markProSide[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1]
												&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

												st.push((indx) * (this->N_cell_y + 1) + indy - 1);
												visited[(indx) * (this->N_cell_y + 1) + indy - 1] = 1;
												mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = count;
												node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
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
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (iz != this->nz - 1) {   // this node is not on the top plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indx != 0) {    // this node is not on the left plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indx != this->nx - 1) {    // this node is not on the right plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indy != 0) {    // this node is not on the front plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
									if (indy != this->ny - 1) {   // this node is not on the back plane
										eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
										compute_edgelink(eno, node1, node2);
										if (node1 != ndi) {
											if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
										else if (node2 != ndi) {
											if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
												v0d1num++;
											}
										}
									}
								}
								leng_v0d1++;
								count++;
							}
						}
					}
				}
			}
			free(visited); visited = NULL;
		}
		leng_v0dd = leng_v0d1;
		myint indj;
		myint inz, inx, iny;
		myint iz, ix, iy;

		/* V0d2 generation */
		queue<myint> qu;
		//visited = (int*)calloc(this->N_node, sizeof(int));
		t2 = clock();

		int markref = 0;   // mark is one reference conductor is excluded in V0d2
		for (indi = 0; indi < this->numCdt; indi++) {
			//cout << this->conductor[indi].markPort << " ";
			if (this->conductor[indi].markPort <= -1 && markref == 0) {   // if this conductor is the reference conductor, no V0d2 corresponding to it
				markref = 1;
				continue;
			}
			else {
				mark = 0;    // if mark = 0 it means that no V0d2 for this conductor, leng_v0d doesn't increase by 1
				//v.clear();
				//va.clear();
				for (indj = 0; indj < this->cdtNumNode[indi]; indj++) {
					mapd[this->conductor[indi].node[indj]] = count;
					iz = this->conductor[indi].node[indj] / this->N_node_s;
					ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
					iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
					avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
					if (iz != 0) {    // this node is not on the bottom plane
						eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
						if (this->markEdge[eno] == 0 && this->markNode[(iz - 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this edge is in the dielectric
							v0d1num++;
							mark = 1;
						}
					}
					if (iz != this->nz - 1) {   // this node is not on the top plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
						if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0d1num++;
							mark = 1;
						}
					}
					if (ix != 0) {    // this node is not on the left plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0d1num++;
							mark = 1;
						}
					}
					if (ix != this->nx - 1) {    // this node is not on the right plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0d1num++;
							mark = 1;
						}
					}
					if (iy != 0) {    // this node is not on the front plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0d1num++;
							mark = 1;
						}
					}
					if (iy != this->ny - 1) {   // this node is not on the back plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0d1num++;
							mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
						}
					}
				}
				if (mark == 1) {
					leng_v0d1++;
					count++;
				}

			}
		}

		/* V0b generation */
		//visited = (int*)calloc(this->N_node, sizeof(int));
		//for (indi = 0; indi < this->numCdt; indi++) {
		//	for (indj = 0; indj < this->cdtNumNode[indi]; indj++) {
		//		iz = this->conductor[indi].node[indj] / this->N_node_s;
		//		ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
		//		iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
		//		avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
		//		mark = 0;

		//		if (iz != 0) {    // this node is not on the bottom plane
		//			eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
		//			if (this->markEdge[eno] == 0) {    // this edge is in the dielectric
		//				mark = 1;
		//			}
		//		}
		//		if (iz != this->nz - 1) {   // this node is not on the top plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (ix != 0) {    // this node is not on the left plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (ix != this->nx - 1) {    // this node is not on the right plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (iy != 0) {    // this node is not on the front plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (iy != this->ny - 1) {   // this node is not on the back plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}

		//		if (mark == 1) {
		//			if (iz != 0) {    // this node is not on the bottom plane
		//				eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
		//				if (this->markEdge[eno] == 0) {    // this edge is in the dielectric; the other node is outside the conductor
		//					v0d1num++;
		//				}
		//			}
		//			if (iz != this->nz - 1) {   // this node is not on the top plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
		//				if (this->markEdge[eno] == 0) {
		//					v0d1num++;
		//				}
		//			}
		//			if (ix != 0) {    // this node is not on the left plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
		//				if (this->markEdge[eno] == 0) {
		//					v0d1num++;
		//				}
		//			}
		//			if (ix != this->nx - 1) {    // this node is not on the right plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
		//				if (this->markEdge[eno] == 0) {
		//					v0d1num++;
		//				}
		//			}
		//			if (iy != 0) {    // this node is not on the front plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
		//				if (this->markEdge[eno] == 0) {
		//					v0d1num++;
		//				}
		//			}
		//			if (iy != this->ny - 1) {   // this node is not on the back plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
		//				if (this->markEdge[eno] == 0) {
		//					v0d1num++;
		//				}
		//			}
		//			leng_v0d1++;
		//		}
		//	}
		//}


		/* Sparse matrix construction for V0d1 */
		this->v0d1RowId = (myint*)malloc(v0d1num * sizeof(myint));
		this->v0d1ColId = (myint*)malloc(v0d1num * sizeof(myint));
		this->v0d1val = (double*)malloc(v0d1num * sizeof(double));
		this->v0d1aval = (double*)malloc(v0d1num * sizeof(double));
		this->v0d1aval1 = (double*)malloc(v0d1num * sizeof(double));
		this->v0dn = (double*)calloc(leng_v0d1, sizeof(double));
		this->v0dan = (double*)calloc(leng_v0d1, sizeof(double));

		double lx_whole_avg = 0;
		double ly_whole_avg = 0;
		double lz_whole_avg = 0;
		lx_whole_avg = (this->xn[this->nx - 1] - this->xn[0]) / (this->nx - 1);
		ly_whole_avg = (this->yn[this->ny - 1] - this->yn[0]) / (this->ny - 1);
		lz_whole_avg = (this->zn[this->nz - 1] - this->zn[0]) / (this->nz - 1);
		leng_v0d1 = 0;
		leng_v0d1a = 0;
		v0d1num = 0;
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
							this->v0dn[leng_v0d1] += pow(1 / (this->zn[iz] - this->zn[iz - 1]), 2);
							this->v0d1aval[v0d1num] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = -1 / lz_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lz_avg, 2);
							v0d1num++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0d1RowId[v0d1num] = eno;
							this->v0d1ColId[v0d1num] = leng_v0d1;
							this->v0d1val[v0d1num] = -1 / (this->zn[iz] - this->zn[iz - 1]);
							this->v0dn[leng_v0d1] += pow(1 / (this->zn[iz] - this->zn[iz - 1]), 2);
							this->v0d1aval[v0d1num] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = -1 / lz_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lz_avg, 2);
							v0d1num++;
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
							this->v0dn[leng_v0d1] += pow(1 / (this->zn[iz + 1] - this->zn[iz]), 2);
							this->v0d1aval[v0d1num] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = 1 / lz_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lz_avg, 2);
							v0d1num++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0d1RowId[v0d1num] = eno;
							this->v0d1ColId[v0d1num] = leng_v0d1;
							this->v0d1val[v0d1num] = 1 / (this->zn[iz + 1] - this->zn[iz]);
							this->v0dn[leng_v0d1] += pow(1 / (this->zn[iz + 1] - this->zn[iz]), 2);
							this->v0d1aval[v0d1num] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = 1 / lz_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lz_avg, 2);
							v0d1num++;
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
							this->v0dn[leng_v0d1] += pow(1 / (this->xn[indx] - this->xn[indx - 1]), 2);
							this->v0d1aval[v0d1num] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = -1 / lx_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lx_avg, 2);
							v0d1num++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0d1RowId[v0d1num] = eno;
							this->v0d1ColId[v0d1num] = leng_v0d1;
							this->v0d1val[v0d1num] = -1 / (this->xn[indx] - this->xn[indx - 1]);
							this->v0dn[leng_v0d1] += pow(1 / (this->xn[indx] - this->xn[indx - 1]), 2);
							this->v0d1aval[v0d1num] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = -1 / lx_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lx_avg, 2);
							v0d1num++;
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
							this->v0dn[leng_v0d1] += pow(1 / (this->xn[indx + 1] - this->xn[indx]), 2);
							this->v0d1aval[v0d1num] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = 1 / lx_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lx_avg, 2);
							v0d1num++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0d1RowId[v0d1num] = eno;
							this->v0d1ColId[v0d1num] = leng_v0d1;
							this->v0d1val[v0d1num] = 1 / (this->xn[indx + 1] - this->xn[indx]);
							this->v0dn[leng_v0d1] += pow(1 / (this->xn[indx + 1] - this->xn[indx]), 2);
							this->v0d1aval[v0d1num] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = 1 / lx_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / lx_avg, 2);
							v0d1num++;
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
							this->v0dn[leng_v0d1] += pow(1 / (this->yn[indy] - this->yn[indy - 1]), 2);
							this->v0d1aval[v0d1num] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = -1 / ly_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / ly_avg, 2);
							v0d1num++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0d1RowId[v0d1num] = eno;
							this->v0d1ColId[v0d1num] = leng_v0d1;
							this->v0d1val[v0d1num] = -1 / (this->yn[indy] - this->yn[indy - 1]);
							this->v0dn[leng_v0d1] += pow(1 / (this->yn[indy] - this->yn[indy - 1]), 2);
							this->v0d1aval[v0d1num] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = -1 / ly_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / ly_avg, 2);
							v0d1num++;
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
							this->v0dn[leng_v0d1] += pow(1 / (this->yn[indy + 1] - this->yn[indy]), 2);
							this->v0d1aval[v0d1num] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = 1 / ly_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / ly_avg, 2);
							v0d1num++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0d1RowId[v0d1num] = eno;
							this->v0d1ColId[v0d1num] = leng_v0d1;
							this->v0d1val[v0d1num] = 1 / (this->yn[indy + 1] - this->yn[indy]);
							this->v0dn[leng_v0d1] += pow(1 / (this->yn[indy + 1] - this->yn[indy]), 2);
							this->v0d1aval[v0d1num] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0d1aval1[v0d1num] = 1 / ly_avg / MU;
							this->v0dan[leng_v0d1] += pow(1 / ly_avg, 2);
							v0d1num++;
						}
					}
				}
			}
			leng_v0d1++;
			leng_v0d1a++;
		}
		int start;
		markref = 0;
		double scalar = 1;   // debug use : check with scalar is better
		/* Note: should make the mesh around the conductor with the same size in order to make sure V0da for each conductor is generated correctly */
		for (indi = 0; indi < this->numCdt; indi++) {
			if (this->conductor[indi].markPort <= -1 && markref == 0) {
				markref = 1;
				continue;
			}
			
			start = v0d1num;
			mark = 0;
			myint cnode = this->cdtNumNode[indi];
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
						this->v0d1aval[v0d1num] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / lz_avg;
						this->v0d1aval1[v0d1num] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / MU;// -1 / lz_avg / MU;
						v0d1num++;
						mark = 1;
					}
				}
				if (iz != this->nz - 1) {   // this node is not on the top plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
					if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
						this->v0d1RowId[v0d1num] = eno;
						this->v0d1ColId[v0d1num] = leng_v0d1;
						this->v0d1val[v0d1num] = 1 / (this->zn[iz + 1] - this->zn[iz]);
						this->v0d1aval[v0d1num] = lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / lz_avg;
						this->v0d1aval1[v0d1num] = lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / MU;// 1 / lz_avg / MU;
						v0d1num++;
						mark = 1;
					}
				}
				if (ix != 0) {    // this node is not on the left plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
						this->v0d1RowId[v0d1num] = eno;
						this->v0d1ColId[v0d1num] = leng_v0d1;
						this->v0d1val[v0d1num] = -1 / (this->xn[ix] - this->xn[ix - 1]);
						this->v0d1aval[v0d1num] = -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / lx_avg;
						this->v0d1aval1[v0d1num] = -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / MU;// -1 / lx_avg / MU;
						v0d1num++;
						mark = 1;
					}
				}
				if (ix != this->nx - 1) {    // this node is not on the right plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
						this->v0d1RowId[v0d1num] = eno;
						this->v0d1ColId[v0d1num] = leng_v0d1;
						this->v0d1val[v0d1num] = 1 / (this->xn[ix + 1] - this->xn[ix]);
						this->v0d1aval[v0d1num] = ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / lx_avg;
						this->v0d1aval1[v0d1num] = ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / MU;// 1 / lx_avg / MU;
						v0d1num++;
						mark = 1;
					}
				}
				if (iy != 0) {    // this node is not on the front plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
						this->v0d1RowId[v0d1num] = eno;
						this->v0d1ColId[v0d1num] = leng_v0d1;
						this->v0d1val[v0d1num] = -1 / (this->yn[iy] - this->yn[iy - 1]);
						this->v0d1aval[v0d1num] = -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / ly_avg;
						this->v0d1aval1[v0d1num] = -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / MU;// -1 / ly_avg / MU;
						v0d1num++;
						mark = 1;
					}
				}
				if (iy != this->ny - 1) {   // this node is not on the back plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
						this->v0d1RowId[v0d1num] = eno;
						this->v0d1ColId[v0d1num] = leng_v0d1;
						this->v0d1val[v0d1num] = 1 / (this->yn[iy + 1] - this->yn[iy]);
						this->v0d1aval[v0d1num] = lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / ly_avg;
						this->v0d1aval1[v0d1num] = lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / MU;// 1 / ly_avg / MU;
						v0d1num++;
						mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
					}
				}
			}
			if (mark == 1) {
				for (indj = start; indj < v0d1num; indj++) {
					this->v0dn[leng_v0d1] += pow(this->v0d1val[indj], 2);
					this->v0dan[leng_v0d1] += pow(this->v0d1aval[indj], 2);
				}
				leng_v0d1++;
				leng_v0d1a++;
			}
		}

		//for (indi = 0; indi < this->numCdt; indi++) {
		//	for (indj = 0; indj < this->cdtNumNode[indi]; indj++) {
		//		mark = 0;
		//		iz = this->conductor[indi].node[indj] / this->N_node_s;
		//		ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
		//		iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
		//		avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);

		//		if (iz != 0) {    // this node is not on the bottom plane
		//			eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
		//			if (this->markEdge[eno] == 0) {    // this edge is in the dielectric
		//				mark = 1;
		//			}
		//		}
		//		if (iz != this->nz - 1) {   // this node is not on the top plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (ix != 0) {    // this node is not on the left plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (ix != this->nx - 1) {    // this node is not on the right plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (iy != 0) {    // this node is not on the front plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}
		//		if (iy != this->ny - 1) {   // this node is not on the back plane
		//			eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
		//			if (this->markEdge[eno] == 0) {
		//				mark = 1;
		//			}
		//		}

		//		if (mark == 1) {
		//			if (iz != 0) {    // this node is not on the bottom plane
		//				eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
		//				if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
		//					this->v0d1RowId[v0d1num] = eno;
		//					this->v0d1ColId[v0d1num] = leng_v0d1;
		//					this->v0d1val[v0d1num] = -1 / (this->zn[iz] - this->zn[iz - 1]);
		//					this->v0d1aval[v0d1num] = -1 / lz_avg;
		//					v0d1num++;
		//				}
		//			}
		//			if (iz != this->nz - 1) {   // this node is not on the top plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
		//				if (this->markEdge[eno] == 0) {
		//					this->v0d1RowId[v0d1num] = eno;
		//					this->v0d1ColId[v0d1num] = leng_v0d1;
		//					this->v0d1val[v0d1num] = 1 / (this->zn[iz + 1] - this->zn[iz]);
		//					this->v0d1aval[v0d1num] = 1 / lz_avg;

		//					v0d1num++;
		//				}
		//			}
		//			if (ix != 0) {    // this node is not on the left plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
		//				if (this->markEdge[eno] == 0) {
		//					this->v0d1RowId[v0d1num] = eno;
		//					this->v0d1ColId[v0d1num] = leng_v0d1;
		//					this->v0d1val[v0d1num] = -1 / (this->xn[ix] - this->xn[ix - 1]);
		//					this->v0d1aval[v0d1num] = -1 / lx_avg;

		//					v0d1num++;
		//				}
		//			}
		//			if (ix != this->nx - 1) {    // this node is not on the right plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
		//				if (this->markEdge[eno] == 0) {
		//					this->v0d1RowId[v0d1num] = eno;
		//					this->v0d1ColId[v0d1num] = leng_v0d1;
		//					this->v0d1val[v0d1num] = 1 / (this->xn[ix + 1] - this->xn[ix]);
		//					this->v0d1aval[v0d1num] = 1 / lx_avg;

		//					v0d1num++;
		//				}
		//			}
		//			if (iy != 0) {    // this node is not on the front plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
		//				if (this->markEdge[eno] == 0) {
		//					this->v0d1RowId[v0d1num] = eno;
		//					this->v0d1ColId[v0d1num] = leng_v0d1;
		//					this->v0d1val[v0d1num] = -1 / (this->yn[iy] - this->yn[iy - 1]);
		//					this->v0d1aval[v0d1num] = -1 / ly_avg;

		//					v0d1num++;
		//				}
		//			}
		//			if (iy != this->ny - 1) {   // this node is not on the back plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
		//				if (this->markEdge[eno] == 0) {
		//					this->v0d1RowId[v0d1num] = eno;
		//					this->v0d1ColId[v0d1num] = leng_v0d1;
		//					this->v0d1val[v0d1num] = 1 / (this->yn[iy + 1] - this->yn[iy]);
		//					this->v0d1aval[v0d1num] = 1 / ly_avg;

		//					v0d1num++;
		//				}
		//			}
		//			leng_v0d1++;
		//		}
		//	}
		//}


		for (indj = 0; indj < leng_v0d1; ++indj) {   // when using the frequency domain schema should not normalize the vectors
			this->v0dn[indj] = sqrt(this->v0dn[indj]);
			this->v0dan[indj] = sqrt(this->v0dan[indj]);
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

	void merge_v0db(double block1_x, double block1_y, double block2_x, double block2_y, double block3_x, double block3_y, double sideLen) {

		// generate V0db = [V0d, V0b] non-redundant, using tree search to find non-reduendant for each cluster
		int* visited = (int*)calloc(N_node, sizeof(int));
		myint node;
		int iz, ix, iy, mark;
		myint eno, nno, start, startd;
		leng_v0db = 0;
		v0dbnum = 0;
		v0dnum = 0;
		leng_v0d = 0;

		unordered_map<int, int> cond_remnode;   // conductor and the edges it can remove. First put dielectric edges. If no didletric edges put shared edges
		queue<int> q;
		vector<int> condvisited(this->numCdt, 0);   // 1 this conductor is in the queue, 0 hasn't been visited
		for (int ind = 0; ind < this->numCdt; ++ind) {
			if (condvisited[ind] == 0 && cond_remnode.find(ind) == cond_remnode.end()) {   // this conductor hasn't been considered
				q.push(ind);
				condvisited[ind] = 1;
				while (!q.empty()) {
					int cdt = q.front();
					q.pop();
					mark = 0;
					for (int indi = 0; indi < this->cdtNumNode[cdt]; ++indi) {
						node = this->conductor[cdt].node[indi];
						iz = node / this->N_node_s;
						ix = (node % this->N_node_s) / (this->N_cell_y + 1);
						iy = (node % this->N_node_s) % (this->N_cell_y + 1);

						if (iz != 0) {    // this node is not on the bottom plane
							eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
							nno = node - this->N_node_s;   // the lower node #
							if (this->markEdge[eno] == 0 && this->markNode[nno] == 0 && mark == 0) {
								cond_remnode[cdt] = node;
								mark = 1;
							}
							else if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end() && condvisited[this->markNode[nno] - 1] == 0) {
								q.push(this->markNode[nno] - 1);
								condvisited[this->markNode[nno] - 1] = 1;
							}
						}
						if (iz != this->nz - 1) {   // this node is not on the top plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
							nno = node + this->N_node_s;
							if (this->markEdge[eno] == 0 && this->markNode[nno] == 0 && mark == 0) {
								cond_remnode[cdt] = node;
								mark = 1;
							}
							else if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end() && condvisited[this->markNode[nno] - 1] == 0) {
								q.push(this->markNode[nno] - 1);
								condvisited[this->markNode[nno] - 1] = 1;
							}
						}
						if (ix != 0) {    // this node is not on the left plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
							nno = node - this->N_cell_y - 1;
							if (this->markEdge[eno] == 0 && this->markNode[nno] == 0 && mark == 0) {
								cond_remnode[cdt] = node;
								mark = 1;
							}
							else if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end() && condvisited[this->markNode[nno] - 1] == 0) {
								q.push(this->markNode[nno] - 1);
								condvisited[this->markNode[nno] - 1] = 1;
							}
						}
						if (ix != this->nx - 1) {    // this node is not on the right plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
							nno = node + this->N_cell_y + 1;
							if (this->markEdge[eno] == 0 && this->markNode[nno] == 0 && mark == 0) {
								cond_remnode[cdt] = node;
								mark = 1;
							}
							else if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end() && condvisited[this->markNode[nno] - 1] == 0) {
								q.push(this->markNode[nno] - 1);
								condvisited[this->markNode[nno] - 1] = 1;
							}
						}
						if (iy != 0) {    // this node is not on the front plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
							nno = node - 1;
							if (this->markEdge[eno] == 0 && this->markNode[nno] == 0 && mark == 0) {
								cond_remnode[cdt] = node;
								mark = 1;
							}
							else if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end() && condvisited[this->markNode[nno] - 1] == 0) {
								q.push(this->markNode[nno] - 1);
								condvisited[this->markNode[nno] - 1] = 1;
							}
						}
						if (iy != this->ny - 1) {   // this node is not on the back plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
							nno = node + 1;
							if (this->markEdge[eno] == 0 && this->markNode[nno] == 0 && mark == 0) {
								cond_remnode[cdt] = node;
								mark = 1;
							}
							else if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end() && condvisited[this->markNode[nno] - 1] == 0) {
								q.push(this->markNode[nno] - 1);
								condvisited[this->markNode[nno] - 1] = 1;
							}
						}
					}
					if (cond_remnode.find(cdt) == cond_remnode.end()) {   // if all the edges are shared
						for (int indi = 0; indi < this->cdtNumNode[cdt]; ++indi) {
							node = this->conductor[cdt].node[indi];
							iz = node / this->N_node_s;
							ix = (node % this->N_node_s) / (this->N_cell_y + 1);
							iy = (node % this->N_node_s) % (this->N_cell_y + 1);

							if (iz != 0) {    // this node is not on the bottom plane
								eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
								nno = node - this->N_node_s;   // the lower node #
								if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end()) {   // connect to a new conductor
									cond_remnode[cdt] = node;
									break;
								}
							}
							if (iz != this->nz - 1) {   // this node is not on the top plane
								eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
								nno = node + this->N_node_s;
								if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end()) {
									cond_remnode[cdt] = node;
									break;
								}
							}
							if (ix != 0) {    // this node is not on the left plane
								eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
								nno = node - this->N_cell_y - 1;
								if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end()) {
									cond_remnode[cdt] = node;
									break;
								}
							}
							if (ix != this->nx - 1) {    // this node is not on the right plane
								eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
								nno = node + this->N_cell_y + 1;
								if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end()) {
									cond_remnode[cdt] = node;
									break;
								}
							}
							if (iy != 0) {    // this node is not on the front plane
								eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
								nno = node - 1;
								if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end()) {
									cond_remnode[cdt] = node;
									break;
								}
							}
							if (iy != this->ny - 1) {   // this node is not on the back plane
								eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
								nno = node + 1;
								if (this->markNode[nno] != 0 && this->markNode[nno] != this->markNode[node] && cond_remnode.find(this->markNode[nno] - 1) == cond_remnode.end()) {
									cond_remnode[cdt] = node;
									break;
								}
							}
						}
					}
				}
			}
		}
		//for (int ind = 0; ind < this->numCdt; ++ind) {
		//	cout << "Conductor " << ind << " " << cond_remnode[ind] << endl;
		//}


		//for (int ind = 0; ind < N_node; ++ind) {
		//	if (visited[ind] == 0) {
		//		queue<int> q;
		//		q.push(ind);
		//		visited[ind] = 1;
		//		while (!q.empty()) {
		//			node = q.front();
		//			iz = node / this->N_node_s;
		//			ix = (node % this->N_node_s) / (this->N_cell_y + 1);
		//			iy = (node % this->N_node_s) % (this->N_cell_y + 1);
		//			start = v0dbnum;
		//			startd = v0dnum;
		//			if (iz != 0) {    // this node is not on the bottom plane
		//				eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
		//				nno = node - this->N_node_s;   // the lower node #
		//				if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
		//					v0dbnum++;
		//					if (markNode[node] == 0) {
		//						v0dnum++;
		//					}
		//					if (visited[nno] == 0) {
		//						q.push(nno);
		//						visited[nno] = 1;
		//					}
		//				}
		//			}
		//			if (iz != this->nz - 1) {   // this node is not on the top plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
		//				nno = node + this->N_node_s;
		//				if (this->markEdge[eno] == 0) {
		//					v0dbnum++;
		//					if (markNode[node] == 0) {
		//						v0dnum++;
		//					}
		//					if (visited[nno] == 0) {
		//						q.push(nno);
		//						visited[nno] = 1;
		//					}
		//				}
		//			}
		//			if (ix != 0) {    // this node is not on the left plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
		//				nno = node - this->N_cell_y - 1;
		//				if (this->markEdge[eno] == 0) {
		//					v0dbnum++;
		//					if (markNode[node] == 0) {
		//						v0dnum++;
		//					}
		//					if (visited[nno] == 0) {
		//						q.push(nno);
		//						visited[nno] = 1;
		//					}
		//				}
		//			}
		//			if (ix != this->nx - 1) {    // this node is not on the right plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
		//				nno = node + this->N_cell_y + 1;
		//				if (this->markEdge[eno] == 0) {
		//					v0dbnum++;
		//					if (markNode[node] == 0) {
		//						v0dnum++;
		//					}
		//					if (visited[nno] == 0) {
		//						q.push(nno);
		//						visited[nno] = 1;
		//					}
		//				}
		//			}
		//			if (iy != 0) {    // this node is not on the front plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
		//				nno = node - 1;
		//				if (this->markEdge[eno] == 0) {
		//					v0dbnum++;
		//					if (markNode[node] == 0) {
		//						v0dnum++;
		//					}
		//					if (visited[nno] == 0) {
		//						q.push(nno);
		//						visited[nno] = 1;
		//					}
		//				}
		//			}
		//			if (iy != this->ny - 1) {   // this node is not on the back plane
		//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
		//				nno = node + 1;
		//				if (this->markEdge[eno] == 0) {
		//					v0dbnum++;
		//					if (markNode[node] == 0) {
		//						v0dnum++;
		//					}
		//					if (visited[nno] == 0) {
		//						q.push(nno);
		//						visited[nno] = 1;
		//					}
		//				}
		//			}
		//			leng_v0db++;
		//			if (markNode[node] == 0) {
		//				leng_v0d++;
		//			}
		//			q.pop();
		//		}
		//		// delete the last node's vector
		//		v0dbnum = start;
		//		leng_v0db--;
		//		if (markNode[node] == 0) {
		//			leng_v0d--;
		//			v0dnum = startd;
		//		}
		//	}
		//}

		unordered_map<int, unordered_set<int>> m;   //conductor #, edge #
		int count = 1;
		for (int ind = 0; ind < N_node; ++ind) {
			if (visited[ind] == 0) {
				queue<int> q;
				q.push(ind);
				visited[ind] = 1;
				while (!q.empty()) {
					node = q.front();
					iz = node / this->N_node_s;
					ix = (node % this->N_node_s) / (this->N_cell_y + 1);
					iy = (node % this->N_node_s) % (this->N_cell_y + 1);
					start = v0dbnum;
					startd = v0dnum;
					/* first go through its surrounding to see whether this is the last node of the cluster */
					mark = 0;
					if (iz != 0) {    // this node is not on the bottom plane
						eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
						nno = node - this->N_node_s;   // the lower node #
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;   // this is not the last node
							}
						}
					}
					if (iz != this->nz - 1) {   // this node is not on the top plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
						nno = node + this->N_node_s;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (ix != 0) {    // this node is not on the left plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
						nno = node - this->N_cell_y - 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (ix != this->nx - 1) {    // this node is not on the right plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
						nno = node + this->N_cell_y + 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (iy != 0) {    // this node is not on the front plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
						nno = node - 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (iy != this->ny - 1) {   // this node is not on the back plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
						nno = node + 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}


					if (mark == 1 || q.size() > 1) {    // this node is not the last node
						if (iz != 0) {    // this node is not on the bottom plane
							eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
							nno = node - this->N_node_s;   // the lower node #
							if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
								v0dbnum++;
								if (markNode[node] == 0) {
									mapd[node] = count;
									v0dnum++;
								}
								else {
									m[markNode[node] - 1].insert(eno);
									if (markNode[nno]) {
										m[markNode[nno] - 1].insert(eno);
									}
								}
								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (iz != this->nz - 1) {   // this node is not on the top plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
							nno = node + this->N_node_s;
							if (this->markEdge[eno] == 0) {
								v0dbnum++;
								if (markNode[node] == 0) {
									mapd[node] = count;
									v0dnum++;
								}
								else {
									m[markNode[node] - 1].insert(eno);
									if (markNode[nno]) {
										m[markNode[nno] - 1].insert(eno);
									}
								}
								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (ix != 0) {    // this node is not on the left plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
							nno = node - this->N_cell_y - 1;
							if (this->markEdge[eno] == 0) {
								v0dbnum++;
								if (markNode[node] == 0) {
									mapd[node] = count;
									v0dnum++;
								}
								else {
									m[markNode[node] - 1].insert(eno);
									if (markNode[nno]) {
										m[markNode[nno] - 1].insert(eno);
									}
								}
								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (ix != this->nx - 1) {    // this node is not on the right plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
							nno = node + this->N_cell_y + 1;
							if (this->markEdge[eno] == 0) {
								v0dbnum++;
								if (markNode[node] == 0) {
									mapd[node] = count;
									v0dnum++;
								}
								else {
									m[markNode[node] - 1].insert(eno);
									if (markNode[nno]) {
										m[markNode[nno] - 1].insert(eno);
									}
								}
								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (iy != 0) {    // this node is not on the front plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
							nno = node - 1;
							if (this->markEdge[eno] == 0) {
								v0dbnum++;
								if (markNode[node] == 0) {
									mapd[node] = count;
									v0dnum++;
								}
								else {
									m[markNode[node] - 1].insert(eno);
									if (markNode[nno]) {
										m[markNode[nno] - 1].insert(eno);
									}
								}
								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (iy != this->ny - 1) {   // this node is not on the back plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
							nno = node + 1;
							if (this->markEdge[eno] == 0) {
								v0dbnum++;
								if (markNode[node] == 0) {
									mapd[node] = count;
									v0dnum++;
								}
								else {
									m[markNode[node] - 1].insert(eno);
									if (markNode[nno]) {
										m[markNode[nno] - 1].insert(eno);
									}
								}
								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						leng_v0db++;
						if (markNode[node] == 0) {
							leng_v0d++;
							count++;
						}
					}
					
					q.pop();
				}
			}
		}

		int* markCdt = (int*)calloc(this->numCdt, sizeof(int));    // 1 V0b cannot form the V0g for this cdt
		for (int indi = 0; indi < this->numCdt; ++indi) {
			for (int indj = 0; indj < this->cdtNumNode[indi]; ++indj) {
				iz = this->conductor[indi].node[indj] / this->N_node_s;
				ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
				iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
				if (iz != 0) {    // this node is not on the bottom plane
					eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
					if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {   // cannot find this dielectric surroudning edge, then no V0g for this conductor
						markCdt[indi] = 1;
						break;
					}
				}
				if (iz != this->nz - 1) {   // this node is not on the top plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
					if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
						markCdt[indi] = 1;
						break;
					}
				}
				if (ix != 0) {    // this node is not on the left plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
					if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
						markCdt[indi] = 1;
						break;
					}
				}
				if (ix != this->nx - 1) {    // this node is not on the right plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
					if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
						markCdt[indi] = 1;
						break;
					}
				}
				if (iy != 0) {    // this node is not on the front plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
					if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
						markCdt[indi] = 1;
						break;
					}
				}
				if (iy != this->ny - 1) {   // this node is not on the back plane
					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
					if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
						markCdt[indi] = 1;
						break;
					}
				}
			}
		}

		/* V0g generation */
		for (int indi = 0; indi < this->numCdt; ++indi) {
			if (markCdt[indi] == 0) {   // V0g from V0b
				mark = 0;    // if mark = 0 it means that no V0d2 for this conductor, leng_v0d doesn't increase by 1
				for (int indj = 0; indj < this->cdtNumNode[indi]; indj++) {
					mapd[this->conductor[indi].node[indj]] = count;
					iz = this->conductor[indi].node[indj] / this->N_node_s;
					ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
					iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
					
					if (iz != 0) {    // this node is not on the bottom plane
						eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
						if (this->markEdge[eno] == 0 && this->markNode[(iz - 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this edge is in the dielectric
							v0dnum++;
							v0dbnum++;
							mark = 1;
						}
					}
					if (iz != this->nz - 1) {   // this node is not on the top plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
						if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0dnum++;
							v0dbnum++;
							mark = 1;
						}
					}
					if (ix != 0) {    // this node is not on the left plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0dnum++;
							v0dbnum++;
							mark = 1;
						}
					}
					if (ix != this->nx - 1) {    // this node is not on the right plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0dnum++;
							v0dbnum++;
							mark = 1;
						}
					}
					if (iy != 0) {    // this node is not on the front plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0dnum++;
							v0dbnum++;
							mark = 1;
						}
					}
					if (iy != this->ny - 1) {   // this node is not on the back plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							v0dnum++;
							v0dbnum++;
							mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
						}
					}
				}
				if (mark == 1) {
					leng_v0d++;
					leng_v0db++;
					count++;
				}
			}
		}

		this->v0dbRowId = (myint*)malloc(v0dbnum * sizeof(myint));
		this->v0dbColId = (myint*)malloc(v0dbnum * sizeof(myint));
		this->v0dbval = (double*)malloc(v0dbnum * sizeof(double));
		this->v0dbaval = (double*)malloc(v0dbnum * sizeof(double));
		this->v0dn = (double*)calloc(leng_v0db, sizeof(double));
		this->v0dan = (double*)calloc(leng_v0db, sizeof(double));
		myint nnz = v0dbnum, nnzd = v0dnum;
		double lx_whole_avg = 0;
		double ly_whole_avg = 0;
		double lz_whole_avg = 0;
		lx_whole_avg = (this->xn[this->nx - 1] - this->xn[0]) / (this->nx - 1);
		ly_whole_avg = (this->yn[this->ny - 1] - this->yn[0]) / (this->ny - 1);
		lz_whole_avg = (this->zn[this->nz - 1] - this->zn[0]) / (this->nz - 1);
		v0dbnum = nnzd;    // the starting nnz of V0b
		v0dnum = 0;
		leng_v0db = leng_v0d;
		myint dleng = leng_v0d;
		leng_v0d = 0;
		
		double lx_avg, ly_avg, lz_avg;
		free(visited); visited = (int*)calloc(N_node, sizeof(int));
		
		unordered_map<int, unordered_set<int>> m1;   // conductor #, v0b col #
		for (int ind = 0; ind < N_node; ++ind) {
			if (visited[ind] == 0) {
				queue<int> q;
				q.push(ind);
				visited[ind] = 1;
				int bstart = leng_v0db - dleng;
				int markl = 0;    // whether has dielectric nodes in this group, 1 has dielectric nodes
				while (!q.empty()) {
					node = q.front();
					iz = node / this->N_node_s;
					ix = (node % this->N_node_s) / (this->N_cell_y + 1);
					iy = (node % this->N_node_s) % (this->N_cell_y + 1);
					avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
					start = v0dbnum;
					startd = v0dnum;
					if (markNode[node] == 0) {
						markl = 1;
					}
					
					/* first go through its surrounding to see whether this is the last node of the cluster */
					mark = 0;
					if (iz != 0) {    // this node is not on the bottom plane
						eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
						nno = node - this->N_node_s;   // the lower node #
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;   // this is not the last node
							}
						}
					}
					if (iz != this->nz - 1) {   // this node is not on the top plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
						nno = node + this->N_node_s;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (ix != 0) {    // this node is not on the left plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
						nno = node - this->N_cell_y - 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (ix != this->nx - 1) {    // this node is not on the right plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
						nno = node + this->N_cell_y + 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (iy != 0) {    // this node is not on the front plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
						nno = node - 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					if (iy != this->ny - 1) {   // this node is not on the back plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
						nno = node + 1;
						if (this->markEdge[eno] == 0) {
							if (visited[nno] == 0) {
								mark = 1;
							}
						}
					}
					

					if (mark == 1 || q.size() > 1) {    // this node is not the last node
						if (iz != 0) {    // this node is not on the bottom plane
							eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
							nno = node - this->N_node_s;   // the lower node #
							if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
								if (markNode[node] == 0) {   // node on the dieletric
									this->v0dbRowId[v0dnum] = eno;
									this->v0dbColId[v0dnum] = leng_v0d;
									this->v0dbval[v0dnum] = -1 / (zn[iz] - zn[iz - 1]);
									this->v0dbaval[v0dnum] = -1 / lz_avg;
									this->v0dn[leng_v0d] += pow(1 / (zn[iz] - zn[iz - 1]), 2);
									this->v0dan[leng_v0d] += pow(1 / lz_avg, 2);
									v0dnum++;
								}
								else {    // node on the conductor boundary
									if (markCdt[markNode[node] - 1] == 1 || (markCdt[markNode[node] - 1] == 0 && node != cond_remnode[markNode[node] - 1])) {
										this->v0dbRowId[v0dbnum] = eno;
										this->v0dbColId[v0dbnum] = leng_v0db;
										this->v0dbval[v0dbnum] = -1 / (zn[iz] - zn[iz - 1]);
										this->v0dbaval[v0dbnum] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
										this->v0dn[leng_v0db] += pow(1 / (zn[iz] - zn[iz - 1]), 2);
										this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
										m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
										v0dbnum++;
									}
								}

								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (iz != this->nz - 1) {   // this node is not on the top plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
							nno = node + this->N_node_s;
							if (this->markEdge[eno] == 0) {
								if (markNode[node] == 0) {    // node is on the dielectric
									this->v0dbRowId[v0dnum] = eno;
									this->v0dbColId[v0dnum] = leng_v0d;
									this->v0dbval[v0dnum] = 1 / (zn[iz + 1] - zn[iz]);
									this->v0dbaval[v0dnum] = 1 / lz_avg;
									this->v0dn[leng_v0d] += pow(1 / (zn[iz + 1] - zn[iz]), 2);
									this->v0dan[leng_v0d] += pow(1 / lz_avg, 2);
									v0dnum++;
								}
								else {   // node is on the conductor boundary
									if (markCdt[markNode[node] - 1] == 1 || (markCdt[markNode[node] - 1] == 0 && node != cond_remnode[markNode[node] - 1])) {
										this->v0dbRowId[v0dbnum] = eno;
										this->v0dbColId[v0dbnum] = leng_v0db;
										this->v0dbval[v0dbnum] = 1 / (zn[iz + 1] - zn[iz]);
										this->v0dbaval[v0dbnum] = lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
										this->v0dn[leng_v0db] += pow(1 / (zn[iz + 1] - zn[iz]), 2);
										this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
										m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
										v0dbnum++;
									}
								}

								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (ix != 0) {    // this node is not on the left plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
							nno = node - this->N_cell_y - 1;
							if (this->markEdge[eno] == 0) {
								if (markNode[node] == 0) {    // node is on the dielectric
									this->v0dbRowId[v0dnum] = eno;
									this->v0dbColId[v0dnum] = leng_v0d;
									this->v0dbval[v0dnum] = -1 / (xn[ix] - xn[ix - 1]);
									this->v0dbaval[v0dnum] = -1 / lx_avg;
									this->v0dn[leng_v0d] += pow(1 / (xn[ix] - xn[ix - 1]), 2);
									this->v0dan[leng_v0d] += pow(1 / lx_avg, 2);
									v0dnum++;
								}
								else {    // node is on the conductor boundary
									if (markCdt[markNode[node] - 1] == 1 || (markCdt[markNode[node] - 1] == 0 && node != cond_remnode[markNode[node] - 1])) {
										this->v0dbRowId[v0dbnum] = eno;
										this->v0dbColId[v0dbnum] = leng_v0db;
										this->v0dbval[v0dbnum] = -1 / (xn[ix] - xn[ix - 1]);
										this->v0dbaval[v0dbnum] = -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
										this->v0dn[leng_v0db] += pow(1 / (xn[ix] - xn[ix - 1]), 2);
										this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
										m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
										v0dbnum++;
									}
								}

								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (ix != this->nx - 1) {    // this node is not on the right plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
							nno = node + this->N_cell_y + 1;
							if (this->markEdge[eno] == 0) {
								if (markNode[node] == 0) {    // this node is on the dielectric
									this->v0dbRowId[v0dnum] = eno;
									this->v0dbColId[v0dnum] = leng_v0d;
									this->v0dbval[v0dnum] = 1 / (xn[ix + 1] - xn[ix]);
									this->v0dbaval[v0dnum] = 1 / lx_avg;
									this->v0dn[leng_v0d] += pow(1 / (xn[ix + 1] - xn[ix]), 2);
									this->v0dan[leng_v0d] += pow(1 / lx_avg, 2);

									v0dnum++;
								}
								else {
									if (markCdt[markNode[node] - 1] == 1 || (markCdt[markNode[node] - 1] == 0 && node != cond_remnode[markNode[node] - 1])) {
										this->v0dbRowId[v0dbnum] = eno;
										this->v0dbColId[v0dbnum] = leng_v0db;
										this->v0dbval[v0dbnum] = 1 / (xn[ix + 1] - xn[ix]);
										this->v0dbaval[v0dbnum] = ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
										this->v0dn[leng_v0db] += pow(1 / (xn[ix + 1] - xn[ix]), 2);
										this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
										m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
										v0dbnum++;
									}
								}
								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (iy != 0) {    // this node is not on the front plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
							nno = node - 1;
							if (this->markEdge[eno] == 0) {
								if (markNode[node] == 0) {
									this->v0dbRowId[v0dnum] = eno;
									this->v0dbColId[v0dnum] = leng_v0d;
									this->v0dbval[v0dnum] = -1 / (yn[iy] - yn[iy - 1]);
									this->v0dbaval[v0dnum] = -1 / ly_avg;
									this->v0dn[leng_v0d] += pow(1 / (yn[iy] - yn[iy - 1]), 2);
									this->v0dan[leng_v0d] += pow(1 / ly_avg, 2);

									v0dnum++;
								}
								else {
									if (markCdt[markNode[node] - 1] == 1 || (markCdt[markNode[node] - 1] == 0 && node != cond_remnode[markNode[node] - 1])) {
										this->v0dbRowId[v0dbnum] = eno;
										this->v0dbColId[v0dbnum] = leng_v0db;
										this->v0dbval[v0dbnum] = -1 / (yn[iy] - yn[iy - 1]);
										this->v0dbaval[v0dbnum] = -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
										this->v0dn[leng_v0db] += pow(1 / (yn[iy] - yn[iy - 1]), 2);
										this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
										m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
										v0dbnum++;
									}
								}

								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}
						if (iy != this->ny - 1) {   // this node is not on the back plane
							eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
							nno = node + 1;
							if (this->markEdge[eno] == 0) {
								if (markNode[node] == 0) {
									this->v0dbRowId[v0dnum] = eno;
									this->v0dbColId[v0dnum] = leng_v0d;
									this->v0dbval[v0dnum] = 1 / (yn[iy + 1] - yn[iy]);
									this->v0dbaval[v0dnum] = 1 / ly_avg;
									this->v0dn[leng_v0d] += pow(1 / (yn[iy + 1] - yn[iy]), 2);
									this->v0dan[leng_v0d] += pow(1 / ly_avg, 2);

									v0dnum++;
								}
								else {
									if (markCdt[markNode[node] - 1] == 1 || (markCdt[markNode[node] - 1] == 0 && node != cond_remnode[markNode[node] - 1])) {
										this->v0dbRowId[v0dbnum] = eno;
										this->v0dbColId[v0dbnum] = leng_v0db;
										this->v0dbval[v0dbnum] = 1 / (yn[iy + 1] - yn[iy]);
										this->v0dbaval[v0dbnum] = lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
										this->v0dn[leng_v0db] += pow(1 / (yn[iy + 1] - yn[iy]), 2);
										this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
										m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
										v0dbnum++;
									}
								}

								if (visited[nno] == 0) {
									q.push(nno);
									visited[nno] = 1;
								}
							}
						}

						if (markNode[node] == 0) {
							leng_v0d++;
							//mapd[node] = leng_v0d;
						}
						else if (markCdt[markNode[node] - 1] == 1 || (markCdt[markNode[node] - 1] == 0 && node != cond_remnode[markNode[node] - 1])) {
							leng_v0db++;
						}
						
					}
					else {   // last node
						
						//if the last node's all surrounding nodes are conductor nodes
						if (markl == 0) {
							for (int ii = bstart; ii < leng_v0db - dleng; ++ii) {
								m1[markNode[node] - 1].insert(-(ii + 1));   // all previous node add together
							}
						}
					}
					q.pop();
				}
			}
		}

		cout << "leng_v0dd is " << leng_v0d << endl;
		/* v0g generation */
		double scalar = 1;   // debug use : check with scalar is better
		/* Note: should make the mesh around the conductor with the same size in order to make sure V0da for each conductor is generated correctly */
		for (int indi = 0; indi < this->numCdt; indi++) {
			if (markCdt[indi] == 0) {   // V0g from V0b
				mark = 0;    // if mark = 0 it means that no V0d2 for this conductor, leng_v0d doesn't increase by 1

				start = v0dnum;
				myint cnode = this->cdtNumNode[indi];
				for (int indj = 0; indj < this->cdtNumNode[indi]; indj++) {

					iz = this->conductor[indi].node[indj] / this->N_node_s;
					ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
					iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
					avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
					if (iz != 0) {    // this node is not on the bottom plane
						eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge

						if (this->markEdge[eno] == 0 && this->markNode[(iz - 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this edge is in the dielectric
							this->v0dbRowId[v0dnum] = eno;
							this->v0dbColId[v0dnum] = leng_v0d;
							this->v0dbval[v0dnum] = -1 / (this->zn[iz] - this->zn[iz - 1]);
							this->v0dbaval[v0dnum] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / lz_avg;
							v0dnum++;
							mark = 1;
						}
					}
					if (iz != this->nz - 1) {   // this node is not on the top plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
						if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							this->v0dbRowId[v0dnum] = eno;
							this->v0dbColId[v0dnum] = leng_v0d;
							this->v0dbval[v0dnum] = 1 / (this->zn[iz + 1] - this->zn[iz]);
							this->v0dbaval[v0dnum] = lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / lz_avg;
							v0dnum++;
							mark = 1;
						}
					}
					if (ix != 0) {    // this node is not on the left plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							this->v0dbRowId[v0dnum] = eno;
							this->v0dbColId[v0dnum] = leng_v0d;
							this->v0dbval[v0dnum] = -1 / (this->xn[ix] - this->xn[ix - 1]);
							this->v0dbaval[v0dnum] = -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / lx_avg;
							v0dnum++;
							mark = 1;
						}
					}
					if (ix != this->nx - 1) {    // this node is not on the right plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							this->v0dbRowId[v0dnum] = eno;
							this->v0dbColId[v0dnum] = leng_v0d;
							this->v0dbval[v0dnum] = 1 / (this->xn[ix + 1] - this->xn[ix]);
							this->v0dbaval[v0dnum] = ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / lx_avg;
							v0dnum++;
							mark = 1;
						}
					}
					if (iy != 0) {    // this node is not on the front plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							this->v0dbRowId[v0dnum] = eno;
							this->v0dbColId[v0dnum] = leng_v0d;
							this->v0dbval[v0dnum] = -1 / (this->yn[iy] - this->yn[iy - 1]);
							this->v0dbaval[v0dnum] = -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / ly_avg
							v0dnum++;
							mark = 1;
						}
					}
					if (iy != this->ny - 1) {   // this node is not on the back plane
						eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
						if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
							this->v0dbRowId[v0dnum] = eno;
							this->v0dbColId[v0dnum] = leng_v0d;
							this->v0dbval[v0dnum] = 1 / (this->yn[iy + 1] - this->yn[iy]);
							this->v0dbaval[v0dnum] = lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / ly_avg;
							v0dnum++;
							mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
						}
					}
				}
				if (mark == 1) {
					for (int indj = start; indj < v0dnum; indj++) {
						this->v0dn[leng_v0d] += pow(this->v0dbval[indj], 2);
						this->v0dan[leng_v0d] += pow(this->v0dbaval[indj], 2);
					}
					leng_v0d++;
				}
			}
		}
		
		for (int ind = 0; ind < leng_v0db; ++ind) {
			this->v0dn[ind] = sqrt(this->v0dn[ind]);
			this->v0dan[ind] = sqrt(this->v0dan[ind]);
		}
		//ofstream out;
		//out.open("map.txt", ofstream::trunc | ofstream::out);
		//for (int ind = 0; ind < N_node; ++ind) {
		//	out << mapd[ind] << endl;
		//}
		//out.close();

		free(visited); visited = NULL;
	}

	//void merge_v0db(double block1_x, double block1_y, double block2_x, double block2_y, double block3_x, double block3_y, double sideLen, int* v0bsize) {

	//	int* visited;
	//	clock_t t1;
	//	double t, ta;
	//	double ratio;
	//	double startx, starty;    // the start coordinates of each block
	//	queue<int> st;    // dfs stack
	//	int indsize;
	//	myint indi = 0;
	//	int indx, indy;
	//	int mark;
	//	int indnum;
	//	int* markLayerNode = (int*)calloc(this->N_node_s, sizeof(int));
	//	/* Mark layer nodes from port sides */
	//	for (int indPort = 0; indPort < this->numPorts; indPort++) {

	//		for (int indPortSide = 0; indPortSide < this->portCoor[indPort].portCnd.size(); indPortSide++) {
	//			myint indCdt = this->portCoor[indPort].portCnd[indPortSide] - 1; // Conductor index for this port side

	//			for (int indCdtNode = 0; indCdtNode < this->cdtNumNode[indCdt]; indCdtNode++) {
	//				markLayerNode[this->conductor[indCdt].node[indCdtNode] % (this->N_node_s)] = 1;
	//			}
	//		}
	//	}

	//	leng_v0db = 0;
	//	v0dbnum = 0;

	//	/* First assign a larger number of storage, don't need to calculate the entries twice */
	//	//myint *v0dbRowId = (myint*)malloc(2 * this->outedge * sizeof(myint));
	//	//myint *v0dbColId = (myint*)malloc(2 * this->outedge * sizeof(myint));
	//	//double *v0dbval = (double*)malloc(2 * this->outedge * sizeof(double));
	//	///*myint *v0d1aRowId = (myint*)malloc(this->N_edge * sizeof(myint));
	//	//myint *v0d1aColId = (myint*)malloc(this->N_edge * sizeof(myint));*/
	//	//double *v0dbaval = (double*)malloc(2 * this->outedge * sizeof(double));

	//	/* V0d1 generation */
	//	int count = 1;    /* count which box it is */
	//	clock_t t2 = clock();
	//	unordered_map<myint, double> va, v;
	//	vector<set<myint>> node_group;
	//	set<myint> base;
	//	myint eno, nno;
	//	double lx_avg, ly_avg, lz_avg;
	//	t = 0.;
	//	ta = 0.;
	//	myint node1, node2;
	//	int nodegs;   // node group #
	//	myint indj;
	//	myint inz, inx, iny;
	//	myint iz, ix, iy;
	//	for (int iz = 0; iz < this->nz; iz++) {    // merge on each layer
	//		visited = (int*)calloc(this->nx * this->ny, sizeof(int));
	//		for (int ix = 0; ix < this->nx; ix++) {
	//			for (int iy = 0; iy < this->ny; iy++) {
	//				nno = iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy;
	//				if ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(nno) == this->lbdn.end() && this->ubdn.find(nno) == this->ubdn.end())) {   // if this node is not boundary node
	//					if (visited[ix * (this->N_cell_y + 1) + iy] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] == 0) {
	//						if (markLayerNode[ix * (this->N_cell_y + 1) + iy] == 0 && !this->markProSide[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this point is not visited and it is outside the conductor, not in the projection of the excited conductor
	//							startx = this->xn[ix];
	//							starty = this->yn[iy];
	//							node_group.push_back(base);
	//							nodegs = node_group.size() - 1;
	//							node_group[nodegs].insert(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);
	//							mapd[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
	//							st.push(ix * (this->N_cell_y + 1) + iy);
	//							visited[ix * (this->N_cell_y + 1) + iy] = 1;

	//							while (!st.empty()) {
	//								indx = (st.front()) / (this->N_cell_y + 1);
	//								indy = st.front() % (this->N_cell_y + 1);
	//								if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
	//									if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block1_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block1_y) {   // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() + this->N_cell_y + 1] == 0
	//											&& visited[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
	//											&& markLayerNode[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
	//											&& !this->markProSide[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited, not in the markLayerNode (Projection of the excited conductors), not in the projection side, not among the boundary nodes
	//											st.push((indx + 1) * (this->N_cell_y + 1) + indy);
	//											visited[(indx + 1) * (this->N_cell_y + 1) + indy] = 1;
	//											mapd[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);
	//										}
	//									}
	//								}
	//								if (indx != 0) {    // it must have a left x edge, thus left x node
	//									if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block1_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block1_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() - this->N_cell_y - 1] == 0
	//											&& visited[(indx - 1) * (this->N_cell_y + 1) + indy] == 0
	//											&& markLayerNode[(indx - 1) * (this->N_cell_y + 1) + indy] == 0 && !this->markProSide[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited, not in the markLayerNode (Projection of the excited conductors), not in the projection side, not among the boundary nodes

	//											st.push((indx - 1) * (this->N_cell_y + 1) + indy);
	//											visited[(indx - 1) * (this->N_cell_y + 1) + indy] = 1;
	//											mapd[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
	//										}
	//									}
	//								}
	//								if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
	//									if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block1_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block1_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() + 1] == 0
	//											&& visited[indx * (this->N_cell_y + 1) + indy + 1] == 0
	//											&& markLayerNode[indx * (this->N_cell_y + 1) + indy + 1] == 0
	//											&& !this->markProSide[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx) * (this->N_cell_y + 1) + indy + 1);
	//											visited[(indx) * (this->N_cell_y + 1) + indy + 1] = 1;
	//											mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
	//										}
	//									}
	//								}
	//								if (indy != 0) {    // it must have a closer y edge, thus closer y node
	//									if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block1_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block1_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
	//											&& visited[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
	//											&& markLayerNode[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
	//											&& !this->markProSide[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx) * (this->N_cell_y + 1) + indy - 1);
	//											visited[(indx) * (this->N_cell_y + 1) + indy - 1] = 1;
	//											mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
	//										}
	//									}
	//								}
	//								st.pop();
	//							}

	//							for (auto ndi : node_group[nodegs]) {
	//								indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
	//								indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
	//								avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
	//								if (iz != 0) {    // this node is not on the bottom plane
	//									eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (iz != this->nz - 1) {   // this node is not on the top plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indx != 0) {    // this node is not on the left plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}

	//								}
	//								if (indx != this->nx - 1) {    // this node is not on the right plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indy != 0) {    // this node is not on the front plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indy != this->ny - 1) {   // this node is not on the back plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//							}
	//							leng_v0db++;

	//							count++;
	//						}

	//						else if (markLayerNode[ix * (this->N_cell_y + 1) + iy] == 1 && !this->markProSide[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {//&& this->exciteCdtLayer[iz] == 1) {    // this point is not visited and it is outside the conductor, in the projection of the excited conductor
	//							startx = this->xn[ix];
	//							starty = this->yn[iy];

	//							mapd[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
	//							st.push(ix * (this->N_cell_y + 1) + iy);
	//							visited[ix * (this->N_cell_y + 1) + iy] = 1;
	//							node_group.push_back(base);
	//							nodegs = node_group.size() - 1;
	//							node_group[nodegs].insert(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);
	//							while (!st.empty()) {
	//								mark = 0;
	//								indx = (st.front()) / (this->N_cell_y + 1);
	//								indy = st.front() % (this->N_cell_y + 1);

	//								if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
	//									if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() + this->N_cell_y + 1] == 0
	//											&& visited[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
	//											&& markLayerNode[(indx + 1) * (this->N_cell_y + 1) + indy] == 1
	//											&& !this->markProSide[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited and not the boundary nodes

	//											st.push((indx + 1) * (this->N_cell_y + 1) + indy);
	//											visited[(indx + 1) * (this->N_cell_y + 1) + indy] = 1;
	//											mapd[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);
	//										}
	//									}
	//								}

	//								if (indx != 0) {    // it must have a left x edge, thus left x node
	//									if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() - this->N_cell_y - 1] == 0
	//											&& visited[(indx - 1) * (this->N_cell_y + 1) + indy] == 0
	//											&& markLayerNode[(indx - 1) * (this->N_cell_y + 1) + indy] == 1
	//											&& !this->markProSide[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx - 1) * (this->N_cell_y + 1) + indy);
	//											visited[(indx - 1) * (this->N_cell_y + 1) + indy] = 1;
	//											mapd[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
	//										}
	//									}
	//								}
	//								if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
	//									if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block2_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() + 1] == 0
	//											&& visited[indx * (this->N_cell_y + 1) + indy + 1] == 0
	//											&& markLayerNode[indx * (this->N_cell_y + 1) + indy + 1] == 1
	//											&& !this->markProSide[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx) * (this->N_cell_y + 1) + indy + 1);
	//											visited[(indx) * (this->N_cell_y + 1) + indy + 1] = 1;
	//											mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
	//										}
	//									}
	//								}
	//								if (indy != 0) {    // it must have a closer y edge, thus closer y node
	//									if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block2_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
	//											&& visited[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
	//											&& markLayerNode[(indx) * (this->N_cell_y + 1) + indy - 1] == 1
	//											&& !this->markProSide[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx) * (this->N_cell_y + 1) + indy - 1);
	//											visited[(indx) * (this->N_cell_y + 1) + indy - 1] = 1;
	//											mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
	//										}
	//									}
	//								}
	//								st.pop();
	//							}

	//							for (auto ndi : node_group[nodegs]) {
	//								indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
	//								indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
	//								avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
	//								if (iz != 0) {    // this node is not on the bottom plane
	//									eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (iz != this->nz - 1) {   // this node is not on the top plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indx != 0) {    // this node is not on the left plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indx != this->nx - 1) {    // this node is not on the right plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indy != 0) {    // this node is not on the front plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indy != this->ny - 1) {   // this node is not on the back plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//							}
	//							leng_v0db++;

	//							count++;

	//						}

	//						else {
	//							startx = this->xn[ix];
	//							starty = this->yn[iy];

	//							mapd[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy] = count;
	//							st.push(ix * (this->N_cell_y + 1) + iy);
	//							visited[ix * (this->N_cell_y + 1) + iy] = 1;
	//							node_group.push_back(base);
	//							nodegs = node_group.size() - 1;
	//							node_group[nodegs].insert(iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy);

	//							while (!st.empty()) {
	//								mark = 0;
	//								indx = (st.front()) / (this->N_cell_y + 1);
	//								indy = st.front() % (this->N_cell_y + 1);

	//								if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
	//									if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block3_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block3_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() + this->N_cell_y + 1] == 0
	//											&& visited[(indx + 1) * (this->N_cell_y + 1) + indy] == 0
	//											&& this->markProSide[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in the sideLen && this node is not among the boundary nodes

	//											st.push((indx + 1) * (this->N_cell_y + 1) + indy);
	//											visited[(indx + 1) * (this->N_cell_y + 1) + indy] = 1;
	//											mapd[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);
	//										}
	//									}
	//								}
	//								if (indx != 0) {    // it must have a left x edge, thus left x node
	//									if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block3_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block3_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() - this->N_cell_y - 1] == 0
	//											&& visited[(indx - 1) * (this->N_cell_y + 1) + indy] == 0
	//											&& this->markProSide[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx - 1) * (this->N_cell_y + 1) + indy);
	//											visited[(indx - 1) * (this->N_cell_y + 1) + indy] = 1;
	//											mapd[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
	//										}
	//									}
	//								}
	//								if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
	//									if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block3_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block3_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() + 1] == 0
	//											&& visited[indx * (this->N_cell_y + 1) + indy + 1] == 0
	//											&& this->markProSide[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx) * (this->N_cell_y + 1) + indy + 1);
	//											visited[(indx) * (this->N_cell_y + 1) + indy + 1] = 1;
	//											mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
	//										}
	//									}
	//								}
	//								if (indy != 0) {    // it must have a closer y edge, thus closer y node
	//									if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block3_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block3_y) {    // this node is within the block area
	//										if (this->markNode[iz * this->N_node_s + st.front() - 1] == 0
	//											&& visited[(indx) * (this->N_cell_y + 1) + indy - 1] == 0
	//											&& this->markProSide[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1]
	//											&& ((iz != 0 && iz != this->nz - 1) || (this->lbdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->lbdn.end() && this->ubdn.find(iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1) == this->ubdn.end()))) {    // this node is in dielectric and this node is not visited

	//											st.push((indx) * (this->N_cell_y + 1) + indy - 1);
	//											visited[(indx) * (this->N_cell_y + 1) + indy - 1] = 1;
	//											mapd[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = count;
	//											node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
	//										}
	//									}
	//								}
	//								st.pop();
	//							}
	//							for (auto ndi : node_group[nodegs]) {
	//								indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
	//								indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
	//								avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
	//								if (iz != 0) {    // this node is not on the bottom plane
	//									eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (iz != this->nz - 1) {   // this node is not on the top plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indx != 0) {    // this node is not on the left plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indx != this->nx - 1) {    // this node is not on the right plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indy != 0) {    // this node is not on the front plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//								if (indy != this->ny - 1) {   // this node is not on the back plane
	//									eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
	//									compute_edgelink(eno, node1, node2);
	//									if (node1 != ndi) {
	//										if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//									else if (node2 != ndi) {
	//										if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//											v0dbnum++;
	//										}
	//									}
	//								}
	//							}
	//							leng_v0db++;
	//							count++;
	//						}
	//					}
	//				}
	//			}
	//		}
	//		free(visited); visited = NULL;
	//	}

	//	leng_v0dd = leng_v0db;
	//	/* V0g generation */
	//	int markref = 0;   // mark is one reference conductor is excluded in V0d2
	//	for (indi = 0; indi < this->numCdt; indi++) {
	//		//cout << this->conductor[indi].markPort << " ";
	//		if (this->conductor[indi].markPort <= -1 && markref == 0) {   // if this conductor is the reference conductor, no V0d2 corresponding to it
	//			markref = 1;
	//			continue;
	//		}
	//		else {
	//			mark = 0;    // if mark = 0 it means that no V0d2 for this conductor, leng_v0d doesn't increase by 1
	//						 //v.clear();
	//						 //va.clear();
	//			for (indj = 0; indj < this->cdtNumNode[indi]; indj++) {
	//				mapd[this->conductor[indi].node[indj]] = count;
	//				iz = this->conductor[indi].node[indj] / this->N_node_s;
	//				ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
	//				iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
	//				avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
	//				if (iz != 0) {    // this node is not on the bottom plane
	//					eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
	//					if (this->markEdge[eno] == 0 && this->markNode[(iz - 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this edge is in the dielectric
	//						v0dbnum++;
	//						mark = 1;
	//					}
	//				}
	//				if (iz != this->nz - 1) {   // this node is not on the top plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//					if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//						v0dbnum++;
	//						mark = 1;
	//					}
	//				}
	//				if (ix != 0) {    // this node is not on the left plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//						v0dbnum++;
	//						mark = 1;
	//					}
	//				}
	//				if (ix != this->nx - 1) {    // this node is not on the right plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//						v0dbnum++;
	//						mark = 1;
	//					}
	//				}
	//				if (iy != 0) {    // this node is not on the front plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//						v0dbnum++;
	//						mark = 1;
	//					}
	//				}
	//				if (iy != this->ny - 1) {   // this node is not on the back plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//					if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//						v0dbnum++;
	//						mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
	//					}
	//				}
	//			}
	//			if (mark == 1) {
	//				leng_v0db++;
	//				count++;
	//			}
	//		}
	//	}

	//	/* Sparse matrix construction for V0d1 */
	//	this->v0dbRowId = (myint*)malloc(v0dbnum * sizeof(myint));
	//	this->v0dbColId = (myint*)malloc(v0dbnum * sizeof(myint));
	//	this->v0dbval = (double*)malloc(v0dbnum * sizeof(double));
	//	this->v0dbaval = (double*)malloc(v0dbnum * sizeof(double));
	//	this->v0dn = (double*)calloc(leng_v0db, sizeof(double));
	//	this->v0dan = (double*)calloc(leng_v0db, sizeof(double));

	//	double lx_whole_avg = 0;
	//	double ly_whole_avg = 0;
	//	double lz_whole_avg = 0;
	//	lx_whole_avg = (this->xn[this->nx - 1] - this->xn[0]) / (this->nx - 1);
	//	ly_whole_avg = (this->yn[this->ny - 1] - this->yn[0]) / (this->ny - 1);
	//	lz_whole_avg = (this->zn[this->nz - 1] - this->zn[0]) / (this->nz - 1);
	//	leng_v0db = 0;
	//	v0dbnum = 0;
	//	for (nodegs = 0; nodegs < node_group.size(); nodegs++) {
	//		for (auto ndi : node_group[nodegs]) {
	//			iz = ndi / (this->N_node_s);
	//			indx = (ndi % this->N_node_s) / (this->N_cell_y + 1);
	//			indy = (ndi % this->N_node_s) % (this->N_cell_y + 1);
	//			avg_length(iz, indy, indx, lx_avg, ly_avg, lz_avg);
	//			if (iz != 0) {    // this node is not on the bottom plane
	//				eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the lower edge
	//				compute_edgelink(eno, node1, node2);
	//				if (node1 != ndi) {
	//					if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = -1 / (this->zn[iz] - this->zn[iz - 1]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->zn[iz] - this->zn[iz - 1]), 2);
	//						this->v0dbaval[v0dbnum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lz_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//				else if (node2 != ndi) {
	//					if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = -1 / (this->zn[iz] - this->zn[iz - 1]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->zn[iz] - this->zn[iz - 1]), 2);
	//						this->v0dbaval[v0dbnum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lz_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//			}
	//			if (iz != this->nz - 1) {   // this node is not on the top plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + indx * (this->N_cell_y + 1) + indy;    // the upper edge
	//				compute_edgelink(eno, node1, node2);
	//				if (node1 != ndi) {
	//					if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = 1 / (this->zn[iz + 1] - this->zn[iz]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->zn[iz + 1] - this->zn[iz]), 2);
	//						this->v0dbaval[v0dbnum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lz_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//				else if (node2 != ndi) {
	//					if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = 1 / (this->zn[iz + 1] - this->zn[iz]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->zn[iz + 1] - this->zn[iz]), 2);
	//						this->v0dbaval[v0dbnum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lz_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//			}
	//			if (indx != 0) {    // this node is not on the left plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (indx - 1) * (this->N_cell_y + 1) + indy;    // the left edge
	//				compute_edgelink(eno, node1, node2);
	//				if (node1 != ndi) {
	//					if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = -1 / (this->xn[indx] - this->xn[indx - 1]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->xn[indx] - this->xn[indx - 1]), 2);
	//						this->v0dbaval[v0dbnum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lx_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//				else if (node2 != ndi) {
	//					if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = -1 / (this->xn[indx] - this->xn[indx - 1]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->xn[indx] - this->xn[indx - 1]), 2);
	//						this->v0dbaval[v0dbnum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lx_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//			}
	//			if (indx != this->nx - 1) {    // this node is not on the right plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + indx * (this->N_cell_y + 1) + indy;    // the right edge
	//				compute_edgelink(eno, node1, node2);
	//				if (node1 != ndi) {
	//					if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = 1 / (this->xn[indx + 1] - this->xn[indx]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->xn[indx + 1] - this->xn[indx]), 2);
	//						this->v0dbaval[v0dbnum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lx_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//				else if (node2 != ndi) {
	//					if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = 1 / (this->xn[indx + 1] - this->xn[indx]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->xn[indx + 1] - this->xn[indx]), 2);
	//						this->v0dbaval[v0dbnum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / lx_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//			}
	//			if (indy != 0) {    // this node is not on the front plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy - 1;    // the front edge
	//				compute_edgelink(eno, node1, node2);
	//				if (node1 != ndi) {
	//					if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = -1 / (this->yn[indy] - this->yn[indy - 1]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->yn[indy] - this->yn[indy - 1]), 2);
	//						this->v0dbaval[v0dbnum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / ly_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//				else if (node2 != ndi) {
	//					if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = -1 / (this->yn[indy] - this->yn[indy - 1]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->yn[indy] - this->yn[indy - 1]), 2);
	//						this->v0dbaval[v0dbnum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / ly_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//			}
	//			if (indy != this->ny - 1) {   // this node is not on the back plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + indx * this->N_cell_y + indy;    // the back edge
	//				compute_edgelink(eno, node1, node2);
	//				if (node1 != ndi) {
	//					if (node_group[nodegs].find(node1) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = 1 / (this->yn[indy + 1] - this->yn[indy]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->yn[indy + 1] - this->yn[indy]), 2);
	//						this->v0dbaval[v0dbnum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / ly_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//				else if (node2 != ndi) {
	//					if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
	//						this->v0dbRowId[v0dbnum] = eno;
	//						this->v0dbColId[v0dbnum] = leng_v0db;
	//						this->v0dbval[v0dbnum] = 1 / (this->yn[indy + 1] - this->yn[indy]);
	//						this->v0dn[leng_v0db] += pow(1 / (this->yn[indy + 1] - this->yn[indy]), 2);
	//						this->v0dbaval[v0dbnum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//						this->v0dan[leng_v0db] += pow(1 / ly_avg, 2);
	//						v0dbnum++;
	//					}
	//				}
	//			}
	//		}
	//		leng_v0db++;
	//	}

	//	/* v0g generation */
	//	int start;
	//	markref = 0;
	//	double scalar = 1;   // debug use : check with scalar is better
	//						 /* Note: should make the mesh around the conductor with the same size in order to make sure V0da for each conductor is generated correctly */
	//	for (indi = 0; indi < this->numCdt; indi++) {
	//		if (this->conductor[indi].markPort <= -1 && markref == 0) {
	//			markref = 1;
	//			continue;
	//		}

	//		start = v0dbnum;
	//		mark = 0;
	//		myint cnode = this->cdtNumNode[indi];
	//		for (indj = 0; indj < this->cdtNumNode[indi]; indj++) {

	//			iz = this->conductor[indi].node[indj] / this->N_node_s;
	//			ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
	//			iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
	//			avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
	//			if (iz != 0) {    // this node is not on the bottom plane
	//				eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge

	//				if (this->markEdge[eno] == 0 && this->markNode[(iz - 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {    // this edge is in the dielectric
	//					this->v0dbRowId[v0dbnum] = eno;
	//					this->v0dbColId[v0dbnum] = leng_v0db;
	//					this->v0dbval[v0dbnum] = -1 / (this->zn[iz] - this->zn[iz - 1]);
	//					this->v0dbaval[v0dbnum] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / lz_avg;
	//					v0dbnum++;
	//					mark = 1;
	//				}
	//			}
	//			if (iz != this->nz - 1) {   // this node is not on the top plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//				if (this->markEdge[eno] == 0 && this->markNode[(iz + 1) * this->N_node_s + ix * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//					this->v0dbRowId[v0dbnum] = eno;
	//					this->v0dbColId[v0dbnum] = leng_v0db;
	//					this->v0dbval[v0dbnum] = 1 / (this->zn[iz + 1] - this->zn[iz]);
	//					this->v0dbaval[v0dbnum] = lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / lz_avg;
	//					v0dbnum++;
	//					mark = 1;
	//				}
	//			}
	//			if (ix != 0) {    // this node is not on the left plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//				if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix - 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//					this->v0dbRowId[v0dbnum] = eno;
	//					this->v0dbColId[v0dbnum] = leng_v0db;
	//					this->v0dbval[v0dbnum] = -1 / (this->xn[ix] - this->xn[ix - 1]);
	//					this->v0dbaval[v0dbnum] = -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / lx_avg;
	//					v0dbnum++;
	//					mark = 1;
	//				}
	//			}
	//			if (ix != this->nx - 1) {    // this node is not on the right plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//				if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + (ix + 1) * (this->N_cell_y + 1) + iy] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//					this->v0dbRowId[v0dbnum] = eno;
	//					this->v0dbColId[v0dbnum] = leng_v0db;
	//					this->v0dbval[v0dbnum] = 1 / (this->xn[ix + 1] - this->xn[ix]);
	//					this->v0dbaval[v0dbnum] = ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / lx_avg;
	//					v0dbnum++;
	//					mark = 1;
	//				}
	//			}
	//			if (iy != 0) {    // this node is not on the front plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//				if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy - 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//					this->v0dbRowId[v0dbnum] = eno;
	//					this->v0dbColId[v0dbnum] = leng_v0db;
	//					this->v0dbval[v0dbnum] = -1 / (this->yn[iy] - this->yn[iy - 1]);
	//					this->v0dbaval[v0dbnum] = -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// -1 / ly_avg;
	//					v0dbnum++;
	//					mark = 1;
	//				}
	//			}
	//			if (iy != this->ny - 1) {   // this node is not on the back plane
	//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//				if (this->markEdge[eno] == 0 && this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy + 1] != this->markNode[iz * this->N_node_s + ix * (this->N_cell_y + 1) + iy]) {
	//					this->v0dbRowId[v0dbnum] = eno;
	//					this->v0dbColId[v0dbnum] = leng_v0db;
	//					this->v0dbval[v0dbnum] = 1 / (this->yn[iy + 1] - this->yn[iy]);
	//					this->v0dbaval[v0dbnum] = lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg) / scalar;// 1 / ly_avg;
	//					v0dbnum++;
	//					mark = 1;    // mark = 1 means that V0d1 has entries for this conductor, leng_v0d will increase by 1
	//				}
	//			}
	//		}
	//		if (mark == 1) {
	//			for (indj = start; indj < v0dbnum; indj++) {
	//				this->v0dn[leng_v0db] += pow(this->v0dbval[indj], 2);
	//				this->v0dan[leng_v0db] += pow(this->v0dbaval[indj], 2);
	//			}
	//			leng_v0db++;
	//		}
	//	}
	//}

	//void merge_v0b(myint& v0bnum, myint& leng_v0b) {
	//	/* generation of V0b */
	//	int* markcond = (int*)calloc(this->numCdt, sizeof(int));   // mark if this conductor is already considered
	//	unordered_set<int> remnode;   // the node that should be removed in V0b
	//	int ix, iy, iz, mark;
	//	myint node1, node2;
	//	v0bnum = 0;
	//	leng_v0b = 0;

	//	/* V0b generation */
	//	int* visited = (int*)calloc(this->N_node, sizeof(int));
	//	myint node;
	//	myint eno, nno, start, startd;
	//	leng_v0b = 0;
	//	v0bnum = 0;
	//	

	//	for (int ind = 0; ind < N_node; ++ind) {
	//		if (visited[ind] == 0 && this->markNode[ind] == 0) {   // start the tree search from a dielectric node
	//			queue<int> q;
	//			int markc = 0;   // 0 need to remove the last node, 1 don't need to remove the last node
	//			q.push(ind);
	//			visited[ind] = 1;
	//			while (!q.empty()) {
	//				node = q.front();
	//				iz = node / this->N_node_s;
	//				ix = (node % this->N_node_s) / (this->N_cell_y + 1);
	//				iy = (node % this->N_node_s) % (this->N_cell_y + 1);
	//				if (markNode[node] != 0)
	//					start = v0bnum;

	//				if (iz != 0) {    // this node is not on the bottom plane
	//					eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
	//					nno = node - this->N_node_s;   // the lower node #
	//					if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
	//						if (markNode[node] != 0) {
	//							v0bnum++;
	//						}
	//						if (visited[nno] == 0) {
	//							q.push(nno);
	//							visited[nno] = 1;
	//						}
	//					}
	//				}
	//				if (iz != this->nz - 1) {   // this node is not on the top plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//					nno = node + this->N_node_s;
	//					if (this->markEdge[eno] == 0) {
	//						if (markNode[node] != 0) {
	//							v0bnum++;
	//						}
	//						if (visited[nno] == 0) {
	//							q.push(nno);
	//							visited[nno] = 1;
	//						}
	//					}
	//				}
	//				if (ix != 0) {    // this node is not on the left plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//					nno = node - this->N_cell_y - 1;
	//					if (this->markEdge[eno] == 0) {
	//						if (markNode[node] != 0) {
	//							v0bnum++;
	//						}
	//						if (visited[nno] == 0) {
	//							q.push(nno);
	//							visited[nno] = 1;
	//						}
	//					}
	//				}
	//				if (ix != this->nx - 1) {    // this node is not on the right plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//					nno = node + this->N_cell_y + 1;
	//					if (this->markEdge[eno] == 0) {
	//						if (markNode[node] != 0) {
	//							v0bnum++;
	//						}
	//						if (visited[nno] == 0) {
	//							q.push(nno);
	//							visited[nno] = 1;
	//						}
	//					}
	//				}
	//				if (iy != 0) {    // this node is not on the front plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//					nno = node - 1;
	//					if (this->markEdge[eno] == 0) {
	//						if (markNode[node] != 0) {
	//							v0bnum++;
	//						}
	//						if (visited[nno] == 0) {
	//							q.push(nno);
	//							visited[nno] = 1;
	//						}
	//					}
	//				}
	//				if (iy != this->ny - 1) {   // this node is not on the back plane
	//					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//					nno = node + 1;
	//					if (this->markEdge[eno] == 0) {
	//						if (markNode[node] != 0) {
	//							v0bnum++;
	//						}
	//						if (visited[nno] == 0) {
	//							q.push(nno);
	//							visited[nno] = 1;
	//						}
	//					}
	//				}
	//				if (markNode[node] != 0) {
	//					leng_v0b++;
	//				}
	//				q.pop();
	//				if (this->markNode[node] != 0 && markcond[this->markNode[node] - 1] == 0) {
	//					markcond[this->markNode[node] - 1] = 1;
	//					markc = 1;
	//					leng_v0b--;
	//					v0bnum = start;
	//				}
	//			}
	//			if (markc == 0) {
	//				// delete the last node's vector
	//				v0bnum = start;
	//				leng_v0b--;
	//			}
	//		}
	//	}
	//	int* visitedcond = (int*)calloc(this->numCdt, sizeof(int));
	//	for (int i = 0; i < this->numCdt; ++i) {
	//		if (markcond[i] == 1 && visitedcond[i] == 0) {   // start from an outer conductor
	//			queue<int> qu;
	//			qu.push(i);
	//			visitedcond[i] = 1;
	//			while (!qu.empty()) {
	//				int cond = qu.front();
	//				for (int j = 0; j < this->cdtNumNode[cond]; ++j) {
	//					int node = this->conductor[cond].node[j];
	//					if (visited[node] == 0) {
	//						queue<int> q;
	//						q.push(node);
	//						visited[node] = 1;
	//						int markc = 0;   // 0 need to remove the last node, 1 don't need to remove the last node
	//						while (!q.empty()) {
	//							int nod = q.front();
	//							iz = nod / this->N_node_s;
	//							ix = (nod % this->N_node_s) / (this->N_cell_y + 1);
	//							iy = (nod % this->N_node_s) % (this->N_cell_y + 1);
	//							if (this->markNode[nod] != 0)
	//								start = v0bnum;
	//							
	//							if (iz != 0) {    // this node is not on the bottom plane
	//								eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
	//								nno = nod - this->N_node_s;   // the lower node #
	//								if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
	//									if (markNode[nod] != 0) {
	//										v0bnum++;
	//									}
	//									if (visited[nno] == 0) {
	//										q.push(nno);
	//										visited[nno] = 1;
	//									}
	//								}
	//							}
	//							if (iz != this->nz - 1) {   // this node is not on the top plane
	//								eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//								nno = nod + this->N_node_s;
	//								if (this->markEdge[eno] == 0) {
	//									if (markNode[nod] != 0) {
	//										v0bnum++;
	//									}
	//									if (visited[nno] == 0) {
	//										q.push(nno);
	//										visited[nno] = 1;
	//									}
	//								}
	//							}
	//							if (ix != 0) {    // this node is not on the left plane
	//								eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//								nno = nod - this->N_cell_y - 1;
	//								if (this->markEdge[eno] == 0) {
	//									if (markNode[nod] != 0) {
	//										v0bnum++;
	//									}
	//									if (visited[nno] == 0) {
	//										q.push(nno);
	//										visited[nno] = 1;
	//									}
	//								}
	//							}
	//							if (ix != this->nx - 1) {    // this node is not on the right plane
	//								eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//								nno = nod + this->N_cell_y + 1;
	//								if (this->markEdge[eno] == 0) {
	//									if (markNode[nod] != 0) {
	//										v0bnum++;
	//									}
	//									if (visited[nno] == 0) {
	//										q.push(nno);
	//										visited[nno] = 1;
	//									}
	//								}
	//							}
	//							if (iy != 0) {    // this node is not on the front plane
	//								eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//								nno = nod - 1;
	//								if (this->markEdge[eno] == 0) {
	//									if (markNode[nod] != 0) {
	//										v0bnum++;
	//									}
	//									if (visited[nno] == 0) {
	//										q.push(nno);
	//										visited[nno] = 1;
	//									}
	//								}
	//							}
	//							if (iy != this->ny - 1) {   // this node is not on the back plane
	//								eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//								nno = nod + 1;
	//								if (this->markEdge[eno] == 0) {
	//									if (markNode[nod] != 0) {
	//										v0bnum++;
	//									}
	//									if (visited[nno] == 0) {
	//										q.push(nno);
	//										visited[nno] = 1;
	//									}
	//								}
	//							}
	//							if (this->markNode[nod] != 0) {
	//								leng_v0b++;
	//							}
	//							q.pop();
	//							if (this->markNode[nod] != 0 && markcond[this->markNode[nod] - 1] == 0) {
	//								markcond[this->markNode[nod] - 1] = 1;
	//								markc = 1;
	//								leng_v0b--;
	//								v0bnum = start;
	//							}
	//							if (visitedcond[this->markNode[nod] - 1] == 0) {
	//								qu.push(this->markNode[nod] - 1);
	//								visitedcond[this->markNode[nod] - 1] = 1;
	//							}
	//						}
	//						if (markc == 0) {
	//							// delete the last node's vector
	//							v0bnum = start;
	//							leng_v0b--;
	//						}
	//					}
	//				}
	//				qu.pop();
	//			}
	//		}
	//	}




	//	//// generate V0db = [V0d, V0b] non-redundant, using tree search to find non-reduendant for each cluster
	//	//int* visited = (int*)calloc(N_node, sizeof(int));
	//	//myint node;
	//	//int iz, ix, iy, mark;
	//	//myint eno, nno, start, startd;
	//	//leng_v0db = 0;
	//	//v0dbnum = 0;
	//	//v0dnum = 0;
	//	//leng_v0d = 0;

	//	//for (int ind = 0; ind < N_node; ++ind) {
	//	//	if (visited[ind] == 0) {
	//	//		queue<int> q;
	//	//		q.push(ind);
	//	//		visited[ind] = 1;
	//	//		while (!q.empty()) {
	//	//			node = q.front();
	//	//			iz = node / this->N_node_s;
	//	//			ix = (node % this->N_node_s) / (this->N_cell_y + 1);
	//	//			iy = (node % this->N_node_s) % (this->N_cell_y + 1);
	//	//			start = v0dbnum;
	//	//			startd = v0dnum;
	//	//			if (iz != 0) {    // this node is not on the bottom plane
	//	//				eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
	//	//				nno = node - this->N_node_s;   // the lower node #
	//	//				if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
	//	//					v0dbnum++;
	//	//					if (markNode[node] == 0) {
	//	//						v0dnum++;
	//	//					}
	//	//					if (visited[nno] == 0) {
	//	//						q.push(nno);
	//	//						visited[nno] = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (iz != this->nz - 1) {   // this node is not on the top plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//	//				nno = node + this->N_node_s;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					v0dbnum++;
	//	//					if (markNode[node] == 0) {
	//	//						v0dnum++;
	//	//					}
	//	//					if (visited[nno] == 0) {
	//	//						q.push(nno);
	//	//						visited[nno] = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (ix != 0) {    // this node is not on the left plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//	//				nno = node - this->N_cell_y - 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					v0dbnum++;
	//	//					if (markNode[node] == 0) {
	//	//						v0dnum++;
	//	//					}
	//	//					if (visited[nno] == 0) {
	//	//						q.push(nno);
	//	//						visited[nno] = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (ix != this->nx - 1) {    // this node is not on the right plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//	//				nno = node + this->N_cell_y + 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					v0dbnum++;
	//	//					if (markNode[node] == 0) {
	//	//						v0dnum++;
	//	//					}
	//	//					if (visited[nno] == 0) {
	//	//						q.push(nno);
	//	//						visited[nno] = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (iy != 0) {    // this node is not on the front plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//	//				nno = node - 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					v0dbnum++;
	//	//					if (markNode[node] == 0) {
	//	//						v0dnum++;
	//	//					}
	//	//					if (visited[nno] == 0) {
	//	//						q.push(nno);
	//	//						visited[nno] = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (iy != this->ny - 1) {   // this node is not on the back plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//	//				nno = node + 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					v0dbnum++;
	//	//					if (markNode[node] == 0) {
	//	//						v0dnum++;
	//	//					}
	//	//					if (visited[nno] == 0) {
	//	//						q.push(nno);
	//	//						visited[nno] = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			leng_v0db++;
	//	//			if (markNode[node] == 0) {
	//	//				leng_v0d++;
	//	//			}
	//	//			q.pop();
	//	//		}
	//	//		// delete the last node's vector
	//	//		v0dbnum = start;
	//	//		leng_v0db--;
	//	//		if (markNode[node] == 0) {
	//	//			leng_v0d--;
	//	//			v0dnum = startd;
	//	//		}
	//	//	}
	//	//}

	//	//this->v0dbRowId = (myint*)malloc(v0dbnum * sizeof(myint));
	//	//this->v0dbColId = (myint*)malloc(v0dbnum * sizeof(myint));
	//	//this->v0dbval = (double*)malloc(v0dbnum * sizeof(double));
	//	//this->v0dbaval = (double*)malloc(v0dbnum * sizeof(double));
	//	//this->v0dn = (double*)calloc(leng_v0db, sizeof(double));
	//	//this->v0dan = (double*)calloc(leng_v0db, sizeof(double));
	//	//myint nnz = v0dbnum, nnzd = v0dnum;
	//	//double lx_whole_avg = 0;
	//	//double ly_whole_avg = 0;
	//	//double lz_whole_avg = 0;
	//	//lx_whole_avg = (this->xn[this->nx - 1] - this->xn[0]) / (this->nx - 1);
	//	//ly_whole_avg = (this->yn[this->ny - 1] - this->yn[0]) / (this->ny - 1);
	//	//lz_whole_avg = (this->zn[this->nz - 1] - this->zn[0]) / (this->nz - 1);
	//	//v0dbnum = nnzd;    // the starting nnz of V0b
	//	//v0dnum = 0;
	//	//leng_v0db = leng_v0d;
	//	//myint dleng = leng_v0d;
	//	//leng_v0d = 0;

	//	//double lx_avg, ly_avg, lz_avg;
	//	//free(visited); visited = (int*)calloc(N_node, sizeof(int));

	//	//unordered_map<int, unordered_map<int, int>> m;   //conductor #, edge #, V0b col # (replication)
	//	//unordered_map<int, unordered_set<int>> m1;   // conductor #, v0b col #

	//	//for (int ind = 0; ind < N_node; ++ind) {
	//	//	if (visited[ind] == 0) {
	//	//		queue<int> q;
	//	//		q.push(ind);
	//	//		visited[ind] = 1;
	//	//		int bstart = leng_v0db - dleng;
	//	//		int markl = 0;    // whether has dielectric nodes in this group, 1 has dielectric nodes
	//	//		while (!q.empty()) {
	//	//			node = q.front();
	//	//			iz = node / this->N_node_s;
	//	//			ix = (node % this->N_node_s) / (this->N_cell_y + 1);
	//	//			iy = (node % this->N_node_s) % (this->N_cell_y + 1);
	//	//			avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
	//	//			start = v0dbnum;
	//	//			startd = v0dnum;
	//	//			if (markNode[node] == 0) {
	//	//				markl = 1;
	//	//			}
	//	//			/* first go through its surrounding to see whether this is the last node of the cluster */
	//	//			mark = 0;
	//	//			if (iz != 0) {    // this node is not on the bottom plane
	//	//				eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
	//	//				nno = node - this->N_node_s;   // the lower node #
	//	//				if (this->markEdge[eno] == 0) {
	//	//					if (visited[nno] == 0) {
	//	//						mark = 1;   // this is not the last node
	//	//					}
	//	//				}
	//	//			}
	//	//			if (iz != this->nz - 1) {   // this node is not on the top plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//	//				nno = node + this->N_node_s;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					if (visited[nno] == 0) {
	//	//						mark = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (ix != 0) {    // this node is not on the left plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//	//				nno = node - this->N_cell_y - 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					if (visited[nno] == 0) {
	//	//						mark = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (ix != this->nx - 1) {    // this node is not on the right plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//	//				nno = node + this->N_cell_y + 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					if (visited[nno] == 0) {
	//	//						mark = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (iy != 0) {    // this node is not on the front plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//	//				nno = node - 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					if (visited[nno] == 0) {
	//	//						mark = 1;
	//	//					}
	//	//				}
	//	//			}
	//	//			if (iy != this->ny - 1) {   // this node is not on the back plane
	//	//				eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//	//				nno = node + 1;
	//	//				if (this->markEdge[eno] == 0) {
	//	//					if (visited[nno] == 0) {
	//	//						mark = 1;
	//	//					}
	//	//				}
	//	//			}


	//	//			if (mark == 1 || q.size() > 1) {    // this node is not the last node
	//	//				if (iz != 0) {    // this node is not on the bottom plane
	//	//					eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
	//	//					nno = node - this->N_node_s;   // the lower node #
	//	//					if (this->markEdge[eno] == 0) {    // this edge is in the dielectric and the other node is outside the conductor
	//	//						if (markNode[node] == 0) {   // node on the dieletric
	//	//							this->v0dbRowId[v0dnum] = eno;
	//	//							this->v0dbColId[v0dnum] = leng_v0d;
	//	//							this->v0dbval[v0dnum] = -1 / (zn[iz] - zn[iz - 1]);
	//	//							this->v0dbaval[v0dnum] = -1 / lz_avg;
	//	//							this->v0dn[leng_v0d] += pow(1 / (zn[iz] - zn[iz - 1]), 2);
	//	//							this->v0dan[leng_v0d] += pow(1 / lz_avg, 2);
	//	//							v0dnum++;
	//	//						}
	//	//						else {    // node on the conductor boundary
	//	//							this->v0dbRowId[v0dbnum] = eno;
	//	//							this->v0dbColId[v0dbnum] = leng_v0db;
	//	//							this->v0dbval[v0dbnum] = -1 / (zn[iz] - zn[iz - 1]);
	//	//							this->v0dbaval[v0dbnum] = -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//	//							this->v0dn[leng_v0db] += pow(1 / (zn[iz] - zn[iz - 1]), 2);
	//	//							this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
	//	//							m[markNode[node] - 1][eno] = leng_v0db - dleng;
	//	//							m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
	//	//							if (markNode[nno]) {
	//	//								m[markNode[nno] - 1][eno] = leng_v0db - dleng;
	//	//							}
	//	//							v0dbnum++;
	//	//						}

	//	//						if (visited[nno] == 0) {
	//	//							q.push(nno);
	//	//							visited[nno] = 1;
	//	//						}
	//	//					}
	//	//				}
	//	//				if (iz != this->nz - 1) {   // this node is not on the top plane
	//	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//	//					nno = node + this->N_node_s;
	//	//					if (this->markEdge[eno] == 0) {
	//	//						if (markNode[node] == 0) {    // node is on the dielectric
	//	//							this->v0dbRowId[v0dnum] = eno;
	//	//							this->v0dbColId[v0dnum] = leng_v0d;
	//	//							this->v0dbval[v0dnum] = 1 / (zn[iz + 1] - zn[iz]);
	//	//							this->v0dbaval[v0dnum] = 1 / lz_avg;
	//	//							this->v0dn[leng_v0d] += pow(1 / (zn[iz + 1] - zn[iz]), 2);
	//	//							this->v0dan[leng_v0d] += pow(1 / lz_avg, 2);
	//	//							v0dnum++;
	//	//						}
	//	//						else {   // node is on the conductor boundary
	//	//							this->v0dbRowId[v0dbnum] = eno;
	//	//							this->v0dbColId[v0dbnum] = leng_v0db;
	//	//							this->v0dbval[v0dbnum] = 1 / (zn[iz + 1] - zn[iz]);
	//	//							this->v0dbaval[v0dbnum] = lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//	//							this->v0dn[leng_v0db] += pow(1 / (zn[iz + 1] - zn[iz]), 2);
	//	//							this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
	//	//							m[markNode[node] - 1][eno] = leng_v0db - dleng;
	//	//							m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
	//	//							if (markNode[nno]) {
	//	//								m[markNode[nno] - 1][eno] = leng_v0db - dleng;
	//	//							}
	//	//							v0dbnum++;
	//	//						}

	//	//						if (visited[nno] == 0) {
	//	//							q.push(nno);
	//	//							visited[nno] = 1;
	//	//						}
	//	//					}
	//	//				}
	//	//				if (ix != 0) {    // this node is not on the left plane
	//	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//	//					nno = node - this->N_cell_y - 1;
	//	//					if (this->markEdge[eno] == 0) {
	//	//						if (markNode[node] == 0) {    // node is on the dielectric
	//	//							this->v0dbRowId[v0dnum] = eno;
	//	//							this->v0dbColId[v0dnum] = leng_v0d;
	//	//							this->v0dbval[v0dnum] = -1 / (xn[ix] - xn[ix - 1]);
	//	//							this->v0dbaval[v0dnum] = -1 / lx_avg;
	//	//							this->v0dn[leng_v0d] += pow(1 / (xn[ix] - xn[ix - 1]), 2);
	//	//							this->v0dan[leng_v0d] += pow(1 / lx_avg, 2);
	//	//							v0dnum++;
	//	//						}
	//	//						else {    // node is on the conductor boundary
	//	//							this->v0dbRowId[v0dbnum] = eno;
	//	//							this->v0dbColId[v0dbnum] = leng_v0db;
	//	//							this->v0dbval[v0dbnum] = -1 / (xn[ix] - xn[ix - 1]);
	//	//							this->v0dbaval[v0dbnum] = -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//	//							this->v0dn[leng_v0db] += pow(1 / (xn[ix] - xn[ix - 1]), 2);
	//	//							this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
	//	//							m[markNode[node] - 1][eno] = leng_v0db - dleng;
	//	//							m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
	//	//							if (markNode[nno]) {
	//	//								m[markNode[nno] - 1][eno] = leng_v0db - dleng;
	//	//							}
	//	//							v0dbnum++;
	//	//						}

	//	//						if (visited[nno] == 0) {
	//	//							q.push(nno);
	//	//							visited[nno] = 1;
	//	//						}
	//	//					}
	//	//				}
	//	//				if (ix != this->nx - 1) {    // this node is not on the right plane
	//	//					eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//	//					nno = node + this->N_cell_y + 1;
	//	//					if (this->markEdge[eno] == 0) {
	//	//						if (markNode[node] == 0) {    // this node is on the dielectric
	//	//							this->v0dbRowId[v0dnum] = eno;
	//	//							this->v0dbColId[v0dnum] = leng_v0d;
	//	//							this->v0dbval[v0dnum] = 1 / (xn[ix + 1] - xn[ix]);
	//	//							this->v0dbaval[v0dnum] = 1 / lx_avg;
	//	//							this->v0dn[leng_v0d] += pow(1 / (xn[ix + 1] - xn[ix]), 2);
	//	//							this->v0dan[leng_v0d] += pow(1 / lx_avg, 2);

	//	//							v0dnum++;
	//	//						}
	//	//						else {
	//	//							this->v0dbRowId[v0dbnum] = eno;
	//	//							this->v0dbColId[v0dbnum] = leng_v0db;
	//	//							this->v0dbval[v0dbnum] = 1 / (xn[ix + 1] - xn[ix]);
	//	//							this->v0dbaval[v0dbnum] = ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//	//							this->v0dn[leng_v0db] += pow(1 / (xn[ix + 1] - xn[ix]), 2);
	//	//							this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
	//	//							m[markNode[node] - 1][eno] = leng_v0db - dleng;
	//	//							m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
	//	//							if (markNode[nno]) {
	//	//								m[markNode[nno] - 1][eno] = leng_v0db - dleng;
	//	//							}
	//	//							v0dbnum++;
	//	//						}
	//	//						if (visited[nno] == 0) {
	//	//							q.push(nno);
	//	//							visited[nno] = 1;
	//	//						}
	//	//					}
	//	//				}
	//	//				if (iy != 0) {    // this node is not on the front plane
	//	//					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//	//					nno = node - 1;
	//	//					if (this->markEdge[eno] == 0) {
	//	//						if (markNode[node] == 0) {
	//	//							this->v0dbRowId[v0dnum] = eno;
	//	//							this->v0dbColId[v0dnum] = leng_v0d;
	//	//							this->v0dbval[v0dnum] = -1 / (yn[iy] - yn[iy - 1]);
	//	//							this->v0dbaval[v0dnum] = -1 / ly_avg;
	//	//							this->v0dn[leng_v0d] += pow(1 / (yn[iy] - yn[iy - 1]), 2);
	//	//							this->v0dan[leng_v0d] += pow(1 / ly_avg, 2);

	//	//							v0dnum++;
	//	//						}
	//	//						else {
	//	//							this->v0dbRowId[v0dbnum] = eno;
	//	//							this->v0dbColId[v0dbnum] = leng_v0db;
	//	//							this->v0dbval[v0dbnum] = -1 / (yn[iy] - yn[iy - 1]);
	//	//							this->v0dbaval[v0dbnum] = -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//	//							this->v0dn[leng_v0db] += pow(1 / (yn[iy] - yn[iy - 1]), 2);
	//	//							this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
	//	//							m[markNode[node] - 1][eno] = leng_v0db - dleng;
	//	//							m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
	//	//							if (markNode[nno]) {
	//	//								m[markNode[nno] - 1][eno] = leng_v0db - dleng;
	//	//							}
	//	//							v0dbnum++;
	//	//						}

	//	//						if (visited[nno] == 0) {
	//	//							q.push(nno);
	//	//							visited[nno] = 1;
	//	//						}
	//	//					}
	//	//				}
	//	//				if (iy != this->ny - 1) {   // this node is not on the back plane
	//	//					eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//	//					nno = node + 1;
	//	//					if (this->markEdge[eno] == 0) {
	//	//						if (markNode[node] == 0) {
	//	//							this->v0dbRowId[v0dnum] = eno;
	//	//							this->v0dbColId[v0dnum] = leng_v0d;
	//	//							this->v0dbval[v0dnum] = 1 / (yn[iy + 1] - yn[iy]);
	//	//							this->v0dbaval[v0dnum] = 1 / ly_avg;
	//	//							this->v0dn[leng_v0d] += pow(1 / (yn[iy + 1] - yn[iy]), 2);
	//	//							this->v0dan[leng_v0d] += pow(1 / ly_avg, 2);

	//	//							v0dnum++;
	//	//						}
	//	//						else {
	//	//							this->v0dbRowId[v0dbnum] = eno;
	//	//							this->v0dbColId[v0dbnum] = leng_v0db;
	//	//							this->v0dbval[v0dbnum] = 1 / (yn[iy + 1] - yn[iy]);
	//	//							this->v0dbaval[v0dbnum] = lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
	//	//							this->v0dn[leng_v0db] += pow(1 / (yn[iy + 1] - yn[iy]), 2);
	//	//							this->v0dan[leng_v0db] += pow(this->v0dbaval[v0dbnum], 2);
	//	//							m[markNode[node] - 1][eno] = leng_v0db - dleng;
	//	//							m1[markNode[node] - 1].insert(leng_v0db - dleng + 1);
	//	//							if (markNode[nno]) {
	//	//								m[markNode[nno] - 1][eno] = leng_v0db - dleng;
	//	//							}
	//	//							v0dbnum++;
	//	//						}

	//	//						if (visited[nno] == 0) {
	//	//							q.push(nno);
	//	//							visited[nno] = 1;
	//	//						}
	//	//					}
	//	//				}

	//	//				if (markNode[node] == 0) {
	//	//					leng_v0d++;
	//	//					mapd[node] = leng_v0d;
	//	//				}
	//	//				else {
	//	//					leng_v0db++;
	//	//				}
	//	//			}
	//	//			else {   // last node

	//	//					 //if the last node's all surrounding nodes are conductor nodes
	//	//				if (markl == 0) {
	//	//					for (int ii = bstart; ii < leng_v0db - dleng; ++ii) {
	//	//						m1[markNode[node] - 1].insert(-(ii + 1));   // all previous node add together
	//	//					}
	//	//				}
	//	//			}
	//	//			q.pop();
	//	//		}
	//	//	}
	//	//}
	//	//for (int ind = 0; ind < leng_v0db; ++ind) {
	//	//	this->v0dn[ind] = sqrt(this->v0dn[ind]);
	//	//	this->v0dan[ind] = sqrt(this->v0dan[ind]);
	//	//}

	//	//int* markCdt = (int*)calloc(this->numCdt, sizeof(int));    // 1 V0b cannot form the V0g for this cdt
	//	//for (int indi = 0; indi < this->numCdt; ++indi) {
	//	//	for (int indj = 0; indj < this->cdtNumNode[indi]; ++indj) {
	//	//		iz = this->conductor[indi].node[indj] / this->N_node_s;
	//	//		ix = (this->conductor[indi].node[indj] % this->N_node_s) / (this->N_cell_y + 1);
	//	//		iy = (this->conductor[indi].node[indj] % this->N_node_s) % (this->N_cell_y + 1);
	//	//		avg_length(iz, iy, ix, lx_avg, ly_avg, lz_avg);
	//	//		if (iz != 0) {    // this node is not on the bottom plane
	//	//			eno = (iz - 1) * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the lower edge
	//	//			if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {   // cannot find this dielectric surroudning edge, then no V0g for this conductor
	//	//				markCdt[indi] = 1;
	//	//				break;
	//	//			}
	//	//		}
	//	//		if (iz != this->nz - 1) {   // this node is not on the top plane
	//	//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_edge_s + ix * (this->N_cell_y + 1) + iy;    // the upper edge
	//	//			if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
	//	//				markCdt[indi] = 1;
	//	//				break;
	//	//			}
	//	//		}
	//	//		if (ix != 0) {    // this node is not on the left plane
	//	//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (ix - 1) * (this->N_cell_y + 1) + iy;    // the left edge
	//	//			if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
	//	//				markCdt[indi] = 1;
	//	//				break;
	//	//			}
	//	//		}
	//	//		if (ix != this->nx - 1) {    // this node is not on the right plane
	//	//			eno = iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + ix * (this->N_cell_y + 1) + iy;    // the right edge
	//	//			if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
	//	//				markCdt[indi] = 1;
	//	//				break;
	//	//			}
	//	//		}
	//	//		if (iy != 0) {    // this node is not on the front plane
	//	//			eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy - 1;    // the front edge
	//	//			if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
	//	//				markCdt[indi] = 1;
	//	//				break;
	//	//			}
	//	//		}
	//	//		if (iy != this->ny - 1) {   // this node is not on the back plane
	//	//			eno = iz * (this->N_edge_s + this->N_edge_v) + ix * this->N_cell_y + iy;    // the back edge
	//	//			if (this->markEdge[eno] == 0 && m[indi].find(eno) == m[indi].end()) {
	//	//				markCdt[indi] = 1;
	//	//				break;
	//	//			}
	//	//		}
	//	//	}
	//	//}

	//	///* form Xg and Xga, there are the same */
	//	//leng_Xg = 0;
	//	//for (int indi = 0; indi < this->numCdt; ++indi) {
	//	//	if (markCdt[indi] == 0) {   // V0g can be formed
	//	//		for (auto mi : m1[indi]) {
	//	//			leng_Xg++;
	//	//		}
	//	//	}
	//	//}
	//	//XgRowId = (myint*)calloc(leng_Xg, sizeof(myint));
	//	//XgColId = (myint*)calloc(leng_Xg, sizeof(myint));
	//	//Xgval = (double*)calloc(leng_Xg, sizeof(double));
	//	//leng_Xg = 0;
	//	//int count_cdt = 0;
	//	//for (int indi = 0; indi < this->numCdt; ++indi) {
	//	//	if (markCdt[indi] == 0) {   // V0g can be formed
	//	//		for (auto mi : m1[indi]) {
	//	//			XgRowId[leng_Xg] = abs(mi) - 1;
	//	//			XgColId[leng_Xg] = count_cdt;
	//	//			Xgval[leng_Xg] = mi / abs(mi);
	//	//			leng_Xg++;
	//	//		}
	//	//		count_cdt++;
	//	//	}
	//	//}

	//	////ofstream out;
	//	////out.open("map.txt", ofstream::trunc | ofstream::out);
	//	////for (int ind = 0; ind < N_node; ++ind) {
	//	////	out << mapd[ind] << endl;
	//	////}
	//	////out.close();

	//	//free(visited); visited = NULL;
	//}

	void generateMd(myint* map, myint v0d1num, myint v0d1anum, myint leng_v0d1) {
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
					Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]));// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
			}
			else if (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
				inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
				iny = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + (inx + 1) * (this->N_cell_y + 1) + iny;
				if (map[node1] != this->v0d1ColId[indi] + 1 && map[node1] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]));// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
			}
			else {    // this edge is along y axis
				inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
				iny = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny + 1;
				if (map[node1] != this->v0d1ColId[indi] + 1 && map[node1] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]));// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
			}
		}



		leng_Md = 0;
		for (indi = 0; indi < leng_v0d1; indi++) {
			leng_Md += Ad1[indi].size();
		}

		this->MdRowId = (myint*)calloc(leng_Md, sizeof(myint));
		this->MdColId = (myint*)calloc(leng_Md, sizeof(myint));
		this->Mdval = (double*)calloc(leng_Md, sizeof(double));
		myint indj = 0;
		ofstream out;
		//out.open("Ad.txt", std::ofstream::trunc | std::ofstream::out);
		for (indi = 0; indi < leng_v0d1; indi++) {
			vector<pair<myint, double>> v(Ad1[indi].begin(), Ad1[indi].end());
			sort(v.begin(), v.end());
			for (auto adi : v) {
				//if (abs(adi.second) > 1e-8) {
				this->MdRowId[indj] = indi;
				this->MdColId[indj] = adi.first;
				this->Mdval[indj] = adi.second;
				//out << this->AdRowId[indj] << " " << this->AdColId[indj] << " ";
				//out << setprecision(15) << this->Adval[indj] << endl;
				indj++;
				//}
			}
			v.clear();
		}
		Ad1.clear();
		//out.close();

	}

	void generateAdeps_v0d1(myint* map, myint v0d1num, myint v0d1anum, myint leng_v0d1) {
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
					Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]));// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
			}
			else if (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
				inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
				iny = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + (inx + 1) * (this->N_cell_y + 1) + iny;
				if (map[node1] != this->v0d1ColId[indi] + 1 && map[node1] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]));// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
			}
			else {    // this edge is along y axis
				inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
				iny = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny + 1;
				if (map[node1] != this->v0d1ColId[indi] + 1 && map[node1] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node1] - 1] += this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else if (map[node2] != this->v0d1ColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0d1ColId[indi]][map[node2] - 1] += this->v0d1aval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]);// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0d1ColId[indi]][this->v0d1ColId[indi]] += abs(this->v0d1aval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]));// / ((this->v0dan[this->v0d1ColId[indi]]) * (this->v0dn[this->v0d1ColId[indi]]));
				}
			}
		}



		leng_Ad = 0;
		for (indi = 0; indi < leng_v0d1; indi++) {
			leng_Ad += Ad1[indi].size();
		}

		this->AdRowId = (myint*)calloc(leng_Ad, sizeof(myint));
		this->AdColId = (myint*)calloc(leng_Ad, sizeof(myint));
		this->Adval = (double*)calloc(leng_Ad, sizeof(double));
		myint indj = 0;
		ofstream out;
		//out.open("Ad.txt", std::ofstream::trunc | std::ofstream::out);
		for (indi = 0; indi < leng_v0d1; indi++) {
			vector<pair<myint, double>> v(Ad1[indi].begin(), Ad1[indi].end());
			sort(v.begin(), v.end());
			for (auto adi : v) {
				//if (abs(adi.second) > 1e-8) {
				this->AdRowId[indj] = indi;
				this->AdColId[indj] = adi.first;
				this->Adval[indj] = adi.second;
				//out << this->AdRowId[indj] << " " << this->AdColId[indj] << " ";
				//out << setprecision(15) << this->Adval[indj] << endl;
				indj++;
				//}
			}
			v.clear();
		}
		Ad1.clear();
		//out.close();

	}


	/* Ad = normalized(V0da)' * D_eps * normalized(V0d)
	*/
	void generateAdeps_v0db(myint* map, myint v0d1num, myint v0d1anum, myint leng_v0d1) {
		unordered_map<myint, unordered_map<myint, double>> Ad1;
		myint indi, inz, inx, iny, node1, node2;

		for (indi = 0; indi < v0d1num; indi++) {
			if (this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // this edge is along z axis
				inz = this->v0dbRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) / (this->N_cell_y + 1);
				iny = ((this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) % (this->N_cell_y + 1);
				node1 = inz * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
				node2 = (inz + 1) * this->N_node_s + (this->N_cell_y + 1) * inx + iny;
				if (map[node1] != this->v0dbColId[indi] + 1 && map[node1] != 0) {
					Ad1[this->v0dbColId[indi]][map[node1] - 1] += this->v0dbaval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += this->v0dbaval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
				else if (map[node2] != this->v0dbColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0dbColId[indi]][map[node2] - 1] += this->v0dbaval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += this->v0dbaval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += abs(this->v0dbaval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0dbRowId[indi]));// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
			}
			else if (this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
				inz = this->v0dbRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
				iny = ((this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + (inx + 1) * (this->N_cell_y + 1) + iny;
				if (map[node1] != this->v0dbColId[indi] + 1 && map[node1] != 0) {
					Ad1[this->v0dbColId[indi]][map[node1] - 1] += this->v0dbaval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += this->v0dbaval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
				else if (map[node2] != this->v0dbColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0dbColId[indi]][map[node2] - 1] += this->v0dbaval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += this->v0dbaval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += abs(this->v0dbaval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0dbRowId[indi]));// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
			}
			else {    // this edge is along y axis
				inz = this->v0dbRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = (this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
				iny = (this->v0dbRowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny + 1;
				if (map[node1] != this->v0dbColId[indi] + 1 && map[node1] != 0) {
					Ad1[this->v0dbColId[indi]][map[node1] - 1] += this->v0dbaval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[map[node1] - 1]));
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += this->v0dbaval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
				else if (map[node2] != this->v0dbColId[indi] + 1 && map[node2] != 0) {
					Ad1[this->v0dbColId[indi]][map[node2] - 1] += this->v0dbaval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[map[node2] - 1]));
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += this->v0dbaval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0dbRowId[indi]);// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ad1[this->v0dbColId[indi]][this->v0dbColId[indi]] += abs(this->v0dbaval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0dbRowId[indi]));// / ((this->v0dan[this->v0dbColId[indi]]) * (this->v0dn[this->v0dbColId[indi]]));
				}
			}
		}

		leng_Ad = 0;
		for (indi = 0; indi < leng_v0d1; indi++) {
			leng_Ad += Ad1[indi].size();
		}

		this->AdRowId = (myint*)calloc(leng_Ad, sizeof(myint));
		this->AdColId = (myint*)calloc(leng_Ad, sizeof(myint));
		this->Adval = (double*)calloc(leng_Ad, sizeof(double));
		myint indj = 0;
		ofstream out;
		//out.open("Ad.txt", std::ofstream::trunc | std::ofstream::out);
		for (indi = 0; indi < leng_v0d1; indi++) {
			vector<pair<myint, double>> v(Ad1[indi].begin(), Ad1[indi].end());
			sort(v.begin(), v.end());
			for (auto adi : v) {
				//if (abs(adi.second) > 1e-8) {
				this->AdRowId[indj] = indi;
				this->AdColId[indj] = adi.first;
				this->Adval[indj] = adi.second;
				//out << this->AdRowId[indj] + 1 << " " << this->AdColId[indj] + 1 << " ";
				//out << setprecision(15) << this->Adval[indj] << endl;
				indj++;
				//}
			}
			v.clear();
		}
		Ad1.clear();
		//out.close();

	}

	/* Generate V0c */
	void merge_v0c(double block_x, double block_y, double block2_x, double block2_y) {

		int* visited;
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
		int mapc_count = 1;
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
				if (this->conductor[ic].markPort == -2)
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
						mapc[this->conductor[ic].node[jc]] = mapc_count;
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


										st.push((indx + 1) * (this->N_cell_y + 1) + indy);
										visited[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = 1;
										mapc[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);

									}
								}
							}
							if (indx != 0) {    // it must have a left x edge, thus left x node
								if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block_y) {    // this node is within the block area
									if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (this->N_cell_x + 1) * (indx - 1) + indy] == markcond
										&& visited[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] == 0
										&& (this->conductor[ic].markPort == 0 && iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy != this->conductor[ic].node[this->cdtNumNode[ic] - 1]) || (this->conductor[ic].markPort == -1 && iz != 0 && iz != this->nz - 1)) {    // this node is in conductor and this node is not visited


										st.push((indx - 1) * (this->N_cell_y + 1) + indy);
										visited[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = 1;
										mapc[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
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


										st.push((indx) * (this->N_cell_y + 1) + indy + 1);
										visited[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = 1;
										mapc[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
										//mark = 1;

										//continue;
									}
								}
							}
							if (indy != 0) {    // it must have a closer y edge, thus closer y node
								if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block_y) {    // this node is within the block area

									if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * indx + indy - 1] == markcond
										&& visited[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] == 0
										&& (this->conductor[ic].markPort == 0 && iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1 != this->conductor[ic].node[this->cdtNumNode[ic] - 1]) || (this->conductor[ic].markPort == -1 && iz != 0 && iz != this->nz - 1)) {    // this node is in conductor and this node is not visited

										st.push((indx) * (this->N_cell_y + 1) + indy - 1);
										visited[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = 1;
										mapc[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
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
						leng_v0c++;
						mapc_count++;

					}
				}

				if (leng_v0c > this->acu_cnno.back())
					this->acu_cnno.push_back(leng_v0c);

				//free(visited); visited = NULL;
			}
			else {   // markPort > 0 then the conductor is not reference conductor, no need to judge the boundary nodes


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
						mapc[this->conductor[ic].node[jc]] = mapc_count;
						node_group.push_back(base);
						nodegs = node_group.size() - 1;
						node_group[nodegs].insert(this->conductor[ic].node[jc]);

						while (!st.empty()) {

							mark = 0;
							indx = (st.front()) / (this->N_cell_y + 1);
							indy = st.front() % (this->N_cell_y + 1);

							if (indx != this->nx - 1) {    // it must have a right x edge, thus right x node
								if ((this->xn[indx + 1] - startx) >= 0 && (this->xn[indx + 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
									if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (this->N_cell_x + 1) * indx + indy] == markcond && visited[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] == 0 && (iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy != this->conductor[ic].node[this->cdtNumNode[ic] - 1])) {    // this node is in conductor and this node is not visited


										st.push((indx + 1) * (this->N_cell_y + 1) + indy);
										visited[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = 1;
										mapc[iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx + 1) * (this->N_cell_y + 1) + indy);
										//mark = 1;

										//continue;
									}
								}
							}
							if (indx != 0) {    // it must have a left x edge, thus left x node
								if ((this->xn[indx - 1] - startx) >= 0 && (this->xn[indx - 1] - startx) <= block2_x && (this->yn[indy] - starty) >= 0 && (this->yn[indy] - starty) <= block2_y) {    // this node is within the block area
									if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * (this->N_cell_x + 1) + (this->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] == 0 && (iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy != this->conductor[ic].node[this->cdtNumNode[ic] - 1])) {    // this node is in conductor and this node is not visited


										st.push((indx - 1) * (this->N_cell_y + 1) + indy);
										visited[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = 1;
										mapc[iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx - 1) * (this->N_cell_y + 1) + indy);
										//mark = 1;

										//continue;
									}
								}
							}
							if (indy != this->ny - 1) {    // it must have a farther y edge, thus farther y node
								if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy + 1] - starty) >= 0 && (this->yn[indy + 1] - starty) <= block2_y) {    // this node is within the block area
									if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * indx + indy] == markcond && visited[iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1] == 0 && (iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy + 1 != this->conductor[ic].node[this->cdtNumNode[ic] - 1])) {    // this node is in conductor and this node is not visited


										st.push((indx) * (this->N_cell_y + 1) + indy + 1);
										visited[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = 1;
										mapc[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy + 1);
										//mark = 1;

										//continue;
									}
								}
							}
							if (indy != 0) {    // it must have a closer y edge, thus closer y node
								if ((this->xn[indx] - startx) >= 0 && (this->xn[indx] - startx) <= block2_x && (this->yn[indy - 1] - starty) >= 0 && (this->yn[indy - 1] - starty) <= block2_y) {    // this node is within the block area
									if (this->markEdge[iz * (this->N_edge_s + this->N_edge_v) + this->N_cell_y * indx + indy - 1] == markcond && visited[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] == 0 && (iz * this->N_node_s + indx * (this->N_cell_y + 1) + indy - 1 != this->conductor[ic].node[this->cdtNumNode[ic] - 1])) {    // this node is in conductor and this node is not visited

										st.push((indx) * (this->N_cell_y + 1) + indy - 1);
										visited[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = 1;
										mapc[iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1] = mapc_count;
										node_group[nodegs].insert(iz * this->N_node_s + (indx) * (this->N_cell_y + 1) + indy - 1);
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
						leng_v0c++;
						mapc_count++;

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

		this->v0cn = (double*)calloc(leng_v0c, sizeof(double));
		this->v0can = (double*)calloc(leng_v0c, sizeof(double));

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
							this->v0cn[leng_v0c] += pow(1 / (this->zn[iz] - this->zn[iz - 1]), 2);
							v0cnum++;
							this->v0caval[v0canum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lz_avg, 2);
							v0canum++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0cRowId[v0cnum] = eno;
							this->v0cColId[v0cnum] = leng_v0c;
							this->v0cval[v0cnum] = -1 / (this->zn[iz] - this->zn[iz - 1]);
							this->v0cn[leng_v0c] += pow(1 / (this->zn[iz] - this->zn[iz - 1]), 2);
							v0cnum++;
							this->v0caval[v0canum] = -1 / lz_avg; // -lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lz_avg, 2);
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
							this->v0cn[leng_v0c] += pow(1 / (this->zn[iz + 1] - this->zn[iz]), 2);
							v0cnum++;
							this->v0caval[v0canum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lz_avg, 2);
							v0canum++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0cRowId[v0cnum] = eno;
							this->v0cColId[v0cnum] = leng_v0c;
							this->v0cval[v0cnum] = 1 / (this->zn[iz + 1] - this->zn[iz]);
							this->v0cn[leng_v0c] += pow(1 / (this->zn[iz + 1] - this->zn[iz]), 2);
							v0cnum++;
							this->v0caval[v0canum] = 1 / lz_avg; // lx_avg * ly_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lz_avg, 2);
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
							this->v0cn[leng_v0c] += pow(1 / (this->xn[indx] - this->xn[indx - 1]), 2);
							v0cnum++;
							this->v0caval[v0canum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lx_avg, 2);
							v0canum++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0cRowId[v0cnum] = eno;
							this->v0cColId[v0cnum] = leng_v0c;
							this->v0cval[v0cnum] = -1 / (this->xn[indx] - this->xn[indx - 1]);
							this->v0cn[leng_v0c] += pow(1 / (this->xn[indx] - this->xn[indx - 1]), 2);
							v0cnum++;
							this->v0caval[v0canum] = -1 / lx_avg; // -ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lx_avg, 2);
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
							this->v0cn[leng_v0c] += pow(1 / (this->xn[indx + 1] - this->xn[indx]), 2);
							v0cnum++;
							this->v0caval[v0canum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lx_avg, 2);
							v0canum++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0cRowId[v0cnum] = eno;
							this->v0cColId[v0cnum] = leng_v0c;
							this->v0cval[v0cnum] = 1 / (this->xn[indx + 1] - this->xn[indx]);
							this->v0cn[leng_v0c] += pow(1 / (this->xn[indx + 1] - this->xn[indx]), 2);
							v0cnum++;
							this->v0caval[v0canum] = 1 / lx_avg; // ly_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / lx_avg, 2);
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
							this->v0cn[leng_v0c] += pow(1 / (this->yn[indy] - this->yn[indy - 1]), 2);
							v0cnum++;
							this->v0caval[v0canum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / ly_avg, 2);
							v0canum++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0cRowId[v0cnum] = eno;
							this->v0cColId[v0cnum] = leng_v0c;
							this->v0cval[v0cnum] = -1 / (this->yn[indy] - this->yn[indy - 1]);
							this->v0cn[leng_v0c] += pow(1 / (this->yn[indy] - this->yn[indy - 1]), 2);
							v0cnum++;
							this->v0caval[v0canum] = -1 / ly_avg; // -lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / ly_avg, 2);
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
							this->v0cn[leng_v0c] += pow(1 / (this->yn[indy + 1] - this->yn[indy]), 2);
							v0cnum++;
							this->v0caval[v0canum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / ly_avg, 2);
							v0canum++;
						}
					}
					else if (node2 != ndi) {
						if (node_group[nodegs].find(node2) == node_group[nodegs].end()) {
							this->v0cRowId[v0cnum] = eno;
							this->v0cColId[v0cnum] = leng_v0c;
							this->v0cval[v0cnum] = 1 / (this->yn[indy + 1] - this->yn[indy]);
							this->v0cn[leng_v0c] += pow(1 / (this->yn[indy + 1] - this->yn[indy]), 2);
							v0cnum++;
							this->v0caval[v0canum] = 1 / ly_avg; // lx_avg * lz_avg / (lx_whole_avg * ly_whole_avg * lz_whole_avg);
							this->v0can[leng_v0c] += pow(1 / ly_avg, 2);
							v0canum++;
						}
					}
				}
			}

			leng_v0c++;
			leng_v0ca++;
		}
		for (int ind = 0; ind < leng_v0c; ++ind) {   // when in the frequency schema should not normalize the vectors
			this->v0cn[ind] = sqrt(this->v0cn[ind]);
			this->v0can[ind] = sqrt(this->v0can[ind]);
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

	void generateAc(myint* map, myint v0cnum, myint v0canum, myint leng_v0c) {
		/* Ac = normalized(V0ca)'*D_sig*normalized(V0c) */
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
					Ac[this->v0cColId[indi]][map[node1] - 1] += this->v0caval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[map[node1] - 1]));
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
				else if (map[node2] != this->v0cColId[indi] + 1 && map[node2] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][map[node2] - 1] += this->v0caval[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[map[node2] - 1]));
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
				else if (this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += abs(this->v0caval[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * SIGMA);// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
			}
			else if (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
				inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
				iny = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + (inx + 1) * (this->N_cell_y + 1) + iny;
				if (map[node1] != this->v0cColId[indi] + 1 && map[node1] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][map[node1] - 1] += this->v0caval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[map[node1] - 1]));
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
				else if (map[node2] != this->v0cColId[indi] + 1 && map[node2] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][map[node2] - 1] += this->v0caval[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[map[node2] - 1]));
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
				else if (this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += abs(this->v0caval[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * SIGMA);// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
			}
			else {    // this edge is along y axis
				inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
				iny = (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
				node1 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny;
				node2 = inz * this->N_node_s + inx * (this->N_cell_y + 1) + iny + 1;
				if (map[node1] != this->v0cColId[indi] + 1 && map[node1] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][map[node1] - 1] += this->v0caval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[map[node1] - 1]));
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
				else if (map[node2] != this->v0cColId[indi] + 1 && map[node2] != 0 && this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][map[node2] - 1] += this->v0caval[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[map[node2] - 1]));
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += this->v0caval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * SIGMA;// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
				else if (this->markEdge[this->v0cRowId[indi]] != 0) {
					Ac[this->v0cColId[indi]][this->v0cColId[indi]] += abs(this->v0caval[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * SIGMA);// / ((this->v0can[this->v0cColId[indi]]) * (this->v0cn[this->v0cColId[indi]]));
				}
			}
		}
		leng_Ac = 0;
		for (indi = 0; indi < leng_v0c; indi++) {
			leng_Ac += Ac[indi].size();
		}
		this->AcRowId = (myint*)calloc(leng_Ac, sizeof(myint));
		this->AcColId = (myint*)calloc(leng_Ac, sizeof(myint));
		this->Acval = (double*)calloc(leng_Ac, sizeof(double));
		indj = 0;
		int k = 1;
		//ofstream out;
		//out.open("Ac.txt", std::ofstream::out | std::ofstream::trunc);
		for (indi = 0; indi < leng_v0c; indi++) {
			vector<pair<myint, double>> v(Ac[indi].begin(), Ac[indi].end());
			sort(v.begin(), v.end());
			for (auto aci : v) {
				//if (abs(aci.second) > 1e4) {    // normlized(v0ca)'*D_sig*normalized(v0c)
				this->AcRowId[indj] = indi;
				this->AcColId[indj] = aci.first;
				this->Acval[indj] = aci.second;
				//out << this->AcRowId[indj] << " " << this->AcColId[indj] << " " << this->Acval[indj] << endl;
				if (this->AcRowId[indj] >= this->acu_cnno[k]) {
					this->cindex.push_back(indj - 1);
					k++;
				}
				indj++;
				//}
			}
			v.clear();
		}
		//out.close();
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

		for (i = 0; i < 2 * (this->N_edge - this->bden); i++) {
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
		for (i = 0; i < k; i++) {
			for (j = 0; j < k; j++) {
				if (i <= j + 1) {
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
			if (abs(wi[i]) > eps) {
				select[i] = 1;    // calculate the non-zero eigenvalues' eigenvectors
								  //cout << wr[i] << " " << wi[i] << endl;
				mm++;
			}
			else {
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
				if (select[j] == 1) {   // always generate the complex conjugate pairs eigenvectors
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
		for (i = 0; i < mm; i++) {
			for (j = 0; j < (this->N_edge - this->bden); j++) {
				this->Vh[i * (this->N_edge - this->bden) + j].real = 0;
				this->Vh[i * (this->N_edge - this->bden) + j].imag = 0;
				for (int in = 0; in < k; in++) {
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
		double* alphar, *alphai, *beta, *vl, *vr, *lscale, *rscale;
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
				val[i] += getEps(this->mapEdgeR[rowId[i]]);
				if (this->markEdge[this->mapEdgeR[rowId[i]]]) {
					val[i] += dt * SIGMA;
				}
			}
			i++;
		}
	}

	/* Calculate the reference */
	void reference1(int freqNo, int sourcePort, complex<double>* xr) {
		/* Use pardiso to solve the original equation
		freqNo : the No. of the frequency
		sourcePort : the number of the sourcePort
		xr : solution */
		double freq = this->freqNo2freq(freqNo);
		myint size = this->N_edge - this->bden;
		myint* RowId1 = (myint*)malloc((size + 1) * sizeof(myint));
		int count = 0;
		int indi = 0;
		int k = 0;
		complex<double>* valc;
		valc = (complex<double>*)calloc(this->leng_S, sizeof(complex<double>));
		complex<double>* J;
		J = (complex<double>*)calloc((this->N_edge - this->bden), sizeof(complex<double>));
		int indz, indy, temp;

		for (int sourcePortSide = 0; sourcePortSide < this->portCoor[sourcePort].multiplicity; sourcePortSide++) {
			for (int inde = 0; inde < this->portCoor[sourcePort].portEdge[sourcePortSide].size(); inde++) {
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
		while (indi < nnz) {
			start = this->SRowId[indi];
			while (indi < nnz && this->SRowId[indi] == start) {
				valc[indi] += this->Sval[indi]; // val[indi] is real
				if (this->SRowId[indi] == this->SColId[indi]) {
					if (this->markEdge[this->mapEdgeR[this->SRowId[indi]]] != 0) {
						complex<double> addedPart(-(2. * M_PI * freq) * this->getEps(this->SRowId[indi]), SIGMA);
						valc[indi] += (2. * M_PI * freq) * addedPart;
					}
					else {
						complex<double> addedPart(-(2. * M_PI * freq) * this->getEps(this->SRowId[indi]), 0);
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
		void* pt[64];

		/* Pardiso control parameters */
		myint iparm[64];
		myint maxfct, mnum, phase, error, msglvl, solver;
		double dparm[64];
		int v0csin;
		myint perm;

		/* Auxiliary variables */
		char* var;

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
		complex<double>* ddum;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, valc, RowId1, this->SColId, &perm, &nrhs, iparm, &msglvl, J, xr, &error);
		if (error != 0) {
			printf("\nERROR during numerical factorization: %d", error);
			exit(2);
		}

		phase = -1;     // Release internal memory
		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, &ddum, RowId1, this->SColId, &perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

		free(RowId1); RowId1 = NULL;
		free(valc); valc = NULL;
		free(J); J = NULL;

	}

	/* Calculate two node's length and the averaged length around one node */
	void nodeLength(myint node1, myint node2, int no, double& l, double& la) {
		/* node1 : the first node
		node2 : the second node
		no : which node is used to calculate the averaged length
		l : length between node1 and node 2
		la : the averaged length around node1 */
		int indx1, indy1, indz1, indx2, indy2, indz2;
		double lxa, lya, lza;
		compute_node_index(node1, indx1, indy1, indz1);
		compute_node_index(node2, indx2, indy2, indz2);
		if (no == 1) {
			avg_length(indz1, indy1, indx1, lxa, lya, lza);
		}
		else if (no == 2) {
			avg_length(indz2, indy2, indx2, lxa, lya, lza);
		}
		if (indx1 != indx2) {   // the edge is along x direction
			l = xn[indx2] - xn[indx1];
			la = lxa;
		}
		else if (indy1 != indy2) {    // the edge is along y direction
			l = yn[indy2] - yn[indy1];
			la = lya;
		}
		else if (indz1 != indz2) {    // the edge is along z direction
			l = zn[indz2] - zn[indz1];
			la = lza;
		}
	}

	/* Calculate the averaged length */
	void avg_length(int iz, int iy, int ix, double& lx, double& ly, double& lz) {    // given a node, we can know its averaged lengths along x, y, z directions
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
	void compute_edgelink(myint eno, myint& node1, myint& node2) {
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
	void compute_node_index(myint node, int& indx, int& indy, int& indz) {
		/* node : the node number before adding the boundary condition
		indx : the node's x index
		indy : the node's y index
		indz : the node's z index */
		indz = node / this->N_node_s;
		indx = (node % this->N_node_s) / (this->N_cell_y + 1);
		indy = (node % this->N_node_s) % (this->N_cell_y + 1);


	}

	/* Construct Z parameters with V0 and Vh */
	void Construct_Z_V0_Vh(complex<double>* x, int freqNo, int sourcePort) {
		/* x: field distribution, after boundary edges removed
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
			cout << "Z is " << this->x[freqNo * (this->numPorts * this->numPorts) + indPort + this->numPorts * sourcePort] << endl;
		}
	}

	/* Construct Z parameters with V0 */
	void Construct_Z_V0(complex<double>* x, int sourcePort) {
		/* x: field distribution from V0 solution
		sourcePort: port no. */
		myint inz, inx, iny;
		double leng;
		double freq;
		complex<double> Zresult;

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

		if (this->nfreq > 1) {

			for (int id = 1; id < this->nfreq; id++) {

				freq = freqNo2freq(id);

				// Report the results beyond the first and append to storage object
				for (int indi = 0; indi < this->numPorts; indi++) {
					for (int indj = 0; indj < this->numPorts; indj++) {
						Zresult = this->x[indj + indi * this->numPorts].real() + (1i) * this->x[indj + indi * this->numPorts].imag() * this->freqStart * this->freqUnit / freq;
						this->x[id * (this->numPorts * this->numPorts) + indi + indj * this->numPorts] = Zresult;
					}
				}
			}
		}
	}

	double getVoltage(double* x, int sourcePort) {
		double vol = 0;
		int inz, inx, iny;
		double leng;

		int indPortSide = 0; // Only deal with first port side to get voltage
		for (int indEdge = 0; indEdge < this->portCoor[sourcePort].portEdge[indPortSide].size(); indEdge++) {
			myint thisEdge = this->portCoor[sourcePort].portEdge[indPortSide][indEdge];
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

			vol -= x[this->mapEdge[thisEdge]] * leng * (this->portCoor[sourcePort].portDirection[indPortSide] * 1.0); // Accumulating responses due to each response edge line integral (V)
		}

		return vol;
	}


	/* print Z for just V0 */
	void print_z() {
		int indi, indj;
		double freq;
		complex<double> Zresult;

		for (int id = 0; id < nfreq; id++) {
			freq = freqNo2freq(id);
			cout << "Z-parameters at frequency " << freq << endl;
			for (int idi = 0; idi < this->numPorts; idi++) {
				for (int idj = 0; idj < this->numPorts; idj++) {
					cout << this->x[id * (this->numPorts * this->numPorts) + idi + idj * this->numPorts] << " ";
				}
				cout << endl;
			}
		}
	}

	/* freqNo to freq */
	double freqNo2freq(int freqNo) {
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
		/* Generate the matrix [V0a'*D*V0, V0a'*D; D*V0, D+L] rowwise */
		myint i;
		int indx, indy, indz, mark;

		unordered_map<myint, unordered_map<myint, double>> Ll;
		myint indi, inz, inx, iny, node1, node2;

		/* Compute [V0da' * D_eps * V0d, V0da' * D_eps * V0c, V0da' * D_eps] */
		for (indi = 0; indi < v0d1num; indi++) {
			if (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // this edge is along z axis
																								  //cout << "z\n";
				inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) / (this->N_cell_y + 1);
				iny = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) % (this->N_cell_y + 1);
				compute_edgelink(this->v0d1RowId[indi], node1, node2);
				//cout << mapd[node1] << " " << mapd[node2] << " " << this->v0d1ColIdo[indi] + 1 << endl;
				/* Generate V0da' * D_eps * V0d */
				if (mapd[node1] != this->v0d1ColIdo[indi] + 1 && mapd[node1] != 0) {
					//cout << this->v0d1avalo[indi] << endl;
					//cout << this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] << endl;
					//cout << (this->v0dan[this->v0d1ColIdo[indi]]) << " " << (this->v0dn[mapd[node1] - 1]) << endl;
					Ll[this->v0d1ColIdo[indi]][mapd[node1] - 1] += this->v0d1avalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[mapd[node1] - 1]));   // do the normalization
					//cout << this->v0d1avalo[indi] << endl;
					//cout << this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] << endl;
					//cout << (this->v0dan[this->v0d1ColIdo[indi]]) << " " << (this->v0dn[this->v0d1ColIdo[indi]]) << endl;
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += this->v0d1avalo[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));
				}
				else if (mapd[node2] != this->v0d1ColIdo[indi] + 1 && mapd[node2] != 0) {
					//cout << this->v0d1avalo[indi] << endl;
					//cout << this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] << endl;
					//cout << (this->v0dan[this->v0d1ColIdo[indi]]) << " " << (this->v0dn[mapd[node2] - 1]) << endl;
					Ll[this->v0d1ColIdo[indi]][mapd[node2] - 1] += this->v0d1avalo[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[mapd[node2] - 1]));
					//cout << this->v0d1avalo[indi] << endl;
					//cout << this->stackEpsn[(this->v0d1RowId[indi] + this->N_edge_v) / (this->N_edge_s + this->N_edge_v)] << endl;
					//cout << (this->v0dan[this->v0d1ColIdo[indi]]) << " " << (this->v0dn[this->v0d1ColIdo[indi]]) << endl;
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += this->v0d1avalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {

					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += abs(this->v0d1avalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi])) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));

				}
				/* Generate V0da' * D_eps * V0c */
				if (mapc[node1] != mapc[node2]) {   // if mapc[node1] == mapc[node2] this means this edge doesn't exsit in V0c
					if (mapc[node1] != 0) {
						Ll[this->v0d1ColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0d1avalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
					}
					if (mapc[node2] != 0) {
						Ll[this->v0d1ColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += -this->v0d1avalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
					}
				}

				/* Compute V0da' * D_eps */
				Ll[this->v0d1ColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0d1RowId[indi]]] += this->v0d1avalo[indi] * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]));
			}
			else if (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
																																//cout << "x\n";
				inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
				iny = ((this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
				compute_edgelink(this->v0d1RowId[indi], node1, node2);
				/* Generate V0da' * D_eps * V0d */
				if (mapd[node1] != this->v0d1ColIdo[indi] + 1 && mapd[node1] != 0) {
					Ll[this->v0d1ColIdo[indi]][mapd[node1] - 1] += this->v0d1avalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[mapd[node1] - 1]));
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += this->v0d1avalo[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));
				}
				else if (mapd[node2] != this->v0d1ColIdo[indi] + 1 && mapd[node2] != 0) {
					Ll[this->v0d1ColIdo[indi]][mapd[node2] - 1] += this->v0d1avalo[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[mapd[node2] - 1]));
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += this->v0d1avalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += abs(this->v0d1avalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi])) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));
				}
				/* Generate V0da' * D_eps * V0c */
				if (mapc[node1] != mapc[node2]) {   // if mapc[node1] == mapc[node2] this means this edge doesn't exsit in V0c
					if (mapc[node1] != 0) {
						Ll[this->v0d1ColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0d1avalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
					}
					if (mapc[node2] != 0) {
						Ll[this->v0d1ColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += -this->v0d1avalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
					}
				}

				/* Compute V0da' * D_eps */
				Ll[this->v0d1ColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0d1RowId[indi]]] += this->v0d1avalo[indi] * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]));
			}
			else {    // this edge is along y axis
					  //cout << "y\n";
				inz = this->v0d1RowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
				iny = (this->v0d1RowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
				compute_edgelink(this->v0d1RowId[indi], node1, node2);
				/* Generate V0da' * D_eps * V0d */
				if (mapd[node1] != this->v0d1ColIdo[indi] + 1 && mapd[node1] != 0) {
					Ll[this->v0d1ColIdo[indi]][mapd[node1] - 1] += this->v0d1avalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[mapd[node1] - 1]));
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += this->v0d1avalo[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));
				}
				else if (mapd[node2] != this->v0d1ColIdo[indi] + 1 && mapd[node2] != 0) {
					Ll[this->v0d1ColIdo[indi]][mapd[node2] - 1] += this->v0d1avalo[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[mapd[node2] - 1]));
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += this->v0d1avalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[v0d1ColIdo[indi]]));
				}
				else {//if (map[this->edgelink[this->v0d1aRowId[indi] * 2]] == 0 || map[this->edgelink[this->v0d1aRowId[indi] * 2] + 1] == 0) {
					Ll[this->v0d1ColIdo[indi]][this->v0d1ColIdo[indi]] += abs(this->v0d1avalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi])) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0dn[this->v0d1ColIdo[indi]]));
				}
				/* Generate V0da' * D_eps * V0c */
				if (mapc[node1] != mapc[node2]) {   // if mapc[node1] == mapc[node2] this means this edge doesn't exsit in V0c
					if (mapc[node1] != 0) {
						Ll[this->v0d1ColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0d1avalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
					}
					if (mapc[node2] != 0) {
						Ll[this->v0d1ColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += -this->v0d1avalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
					}
				}
				/* Compute V0da' * D_eps */
				Ll[this->v0d1ColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0d1RowId[indi]]] += this->v0d1avalo[indi] * this->getEps(this->v0d1RowId[indi]) / ((this->v0dan[this->v0d1ColIdo[indi]]));
			}
		}

		/* Compute [V0ca' * D_eps * V0d, V0ca' * (D_eps + dt * D_sig) * V0c, V0ca' * (D_eps + dt * D_sig)] */
		for (indi = 0; indi < v0canum; ++indi) {    // the upper and lower planes are PEC
			if (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v) >= this->N_edge_s) {    // this edge is along z axis
				inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) / (this->N_cell_y + 1);
				iny = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - this->N_edge_s) % (this->N_cell_y + 1);
				compute_edgelink(this->v0cRowId[indi], node1, node2);
				/* Compute V0ca' * D_eps * V0d */
				if (mapd[node1] != mapd[node2]) {   // if mapd[node1] == mapd[node2] this means this edge doesn't exsit in V0d
					if (mapd[node1] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][mapd[node1] - 1] += this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0cRowId[indi]) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0dn[mapd[node1] - 1]));
					}
					if (mapd[node2] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][mapd[node2] - 1] += -this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * this->getEps(this->v0cRowId[indi]) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0dn[mapd[node2] - 1]));
					}
				}

				/* Compute V0ca' * (D_eps + dt * D_sig) * V0c */
				if (mapc[node1] != this->v0cColIdo[indi] + 1 && mapc[node1] != 0) {
					if (this->markEdge[this->v0cRowId[indi]]) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
				}
				else if (mapc[node2] != this->v0cColIdo[indi] + 1 && mapc[node2] != 0) {
					if (this->markEdge[this->v0cRowId[indi]]) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += this->v0cavalo[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += this->v0cavalo[indi] * (-1) / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
				}
				else {
					if (this->markEdge[this->v0cRowId[indi]]) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += abs(this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]])));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += abs(this->v0cavalo[indi] * 1 / (this->zn[inz + 1] - this->zn[inz]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]])));
					}
				}

				/* Compute V0ca' * (D_eps + dt * D_sig) */
				if (this->markEdge[this->v0cRowId[indi]])
					Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]] += this->v0cavalo[indi] * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]));
				else
					Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]] += this->v0cavalo[indi] * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]));
			}
			else if (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v) >= (this->N_cell_y) * (this->N_cell_x + 1)) {    // this edge is along x axis
				inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) / (this->N_cell_y + 1);
				iny = ((this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) - (this->N_cell_y) * (this->N_cell_x + 1)) % (this->N_cell_y + 1);
				compute_edgelink(this->v0cRowId[indi], node1, node2);

				/* Compute V0ca' * D_eps * V0d */
				if (mapd[node1] != mapd[node2]) {   // if mapd[node1] == mapd[node2] this means this edge doesn't exsit in V0d
					if (mapd[node1] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][mapd[node1] - 1] += this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0cRowId[indi]) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0dn[mapd[node1] - 1]));
					}
					if (mapd[node2] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][mapd[node2] - 1] += -this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * this->getEps(this->v0cRowId[indi]) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0dn[mapd[node2] - 1]));
					}
				}

				/* Compute V0ca' * (D_eps + dt * D_sig) * V0c */
				if (mapc[node1] != this->v0cColIdo[indi] + 1 && mapc[node1] != 0) {
					if (this->markEdge[this->v0cRowId[indi]]) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
				}
				else if (mapc[node2] != this->v0cColIdo[indi] + 1 && mapc[node2] != 0) {
					if (this->markEdge[this->v0cRowId[indi]]) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += this->v0cavalo[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += this->v0cavalo[indi] * (-1) / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
				}
				else {
					if (this->markEdge[this->v0cRowId[indi]] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += abs(this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]])));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += abs(this->v0cavalo[indi] * 1 / (this->xn[inx + 1] - this->xn[inx]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]])));
					}
				}
				/* Compute V0ca' * (D_eps + dt * D_sig) */
				if (this->markEdge[this->v0cRowId[indi]])
					Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]] += this->v0cavalo[indi] * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]));
				else
					Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]] += this->v0cavalo[indi] * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]));
			}
			else {    // this edge is along y axis
				inz = this->v0cRowId[indi] / (this->N_edge_s + this->N_edge_v);
				inx = (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) / this->N_cell_y;
				iny = (this->v0cRowId[indi] % (this->N_edge_s + this->N_edge_v)) % this->N_cell_y;
				compute_edgelink(this->v0cRowId[indi], node1, node2);
				/* Compute V0ca' *D_eps * V0d */
				if (mapd[node1] != mapd[node2]) {   // if mapd[node1] == mapd[node2] this means this edge doesn't exsit in V0d
					if (mapd[node1] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][mapd[node1] - 1] += this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0cRowId[indi]) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0dn[mapd[node1] - 1]));
					}
					if (mapd[node2] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][mapd[node2] - 1] += -this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * this->getEps(this->v0cRowId[indi]) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0dn[mapd[node2] - 1]));
					}
				}

				/* Compute V0ca' * (D_eps + dt * D_sig) * V0c */
				if (mapc[node1] != this->v0cColIdo[indi] + 1 && mapc[node1] != 0) {
					if (this->markEdge[this->v0cRowId[indi]]) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node1] - 1] += this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node1] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
				}
				else if (mapc[node2] != this->v0cColIdo[indi] + 1 && mapc[node2] != 0) {
					if (this->markEdge[this->v0cRowId[indi]] != 0) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += this->v0cavalo[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + mapc[node2] - 1] += this->v0cavalo[indi] * (-1) / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[mapc[node2] - 1]));
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]]));
					}

				}
				else {
					if (this->markEdge[this->v0cRowId[indi]]) {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += abs(this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]])));
					}
					else {
						Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + this->v0cColIdo[indi]] += abs(this->v0cavalo[indi] * 1 / (this->yn[iny + 1] - this->yn[iny]) * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]) * (this->v0cn[this->v0cColIdo[indi]])));
					}

				}
				/* Compute V0ca' * (D_eps + dt * D_sig) */
				if (this->markEdge[this->v0cRowId[indi]])
					Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]] += this->v0cavalo[indi] * (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) / ((this->v0can[this->v0cColIdo[indi]]));
				else
					Ll[leng_v0d1 + this->v0cColIdo[indi]][leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]] += this->v0cavalo[indi] * (this->getEps(this->v0cRowId[indi])) / ((this->v0can[this->v0cColIdo[indi]]));

			}
		}


		/* Compute [D_eps * V0d, (D_eps + dt * D_sig) * V0c, D_eps + dt * D_sig + dt^2 * L] */
		for (indi = 0; indi < v0d1num; indi++) {
			Ll[leng_v0d1 + leng_v0c + this->mapEdge[this->v0d1RowId[indi]]][this->v0d1ColIdo[indi]] += this->getEps(this->v0d1RowId[indi]) * this->v0d1valo[indi] / ((this->v0dn[this->v0d1ColIdo[indi]]));
			//if (this->v0d1ColIdo[indi] >= 4137)
			//	cout << this->v0d1ColIdo[indi] << " ";
		}

		for (indi = 0; indi < v0cnum; ++indi) {
			if (markEdge[this->v0cRowId[indi]])
				Ll[leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]][leng_v0d1 + this->v0cColIdo[indi]] += (this->getEps(this->v0cRowId[indi]) + dt * SIGMA) * this->v0cvalo[indi] / ((this->v0cn[this->v0cColIdo[indi]]));
			else {
				Ll[leng_v0d1 + leng_v0c + this->mapEdge[this->v0cRowId[indi]]][leng_v0d1 + this->v0cColIdo[indi]] += (this->getEps(this->v0cRowId[indi])) * this->v0cvalo[indi] / ((this->v0cn[this->v0cColIdo[indi]]));
			}
			//if (leng_v0d1 + this->v0cColIdo[indi] >= 4137)
			//	cout << leng_v0d1 + this->v0cColIdo[indi] << " ";
		}

		for (indi = 0; indi < Lleng; ++indi) {
			if (LrowId[indi] == LcolId[indi]) {
				if (markEdge[this->mapEdgeR[LrowId[indi]]])
					Ll[leng_v0d1 + leng_v0c + LrowId[indi]][leng_v0d1 + leng_v0c + LcolId[indi]] += this->getEps(this->mapEdgeR[LrowId[indi]]) + dt * SIGMA + dt * dt * Lval[indi];
				else
					Ll[leng_v0d1 + leng_v0c + LrowId[indi]][leng_v0d1 + leng_v0c + LcolId[indi]] += this->getEps(this->mapEdgeR[LrowId[indi]]) + dt * dt * Lval[indi];
			}
			else {
				Ll[leng_v0d1 + leng_v0c + LrowId[indi]][leng_v0d1 + leng_v0c + LcolId[indi]] += dt * dt * Lval[indi];
			}
			//if (leng_v0d1 + leng_v0c + LcolId[indi] >= 4137)
			//	cout << leng_v0d1 + leng_v0c + LcolId[indi] << " ";
		}
		//cout << endl;
		leng_Ll = 0;
		for (indi = 0; indi < leng_v0d1 + leng_v0c + N_edge - bden; indi++) {
			leng_Ll += Ll[indi].size();
		}

		this->LlRowId = (myint*)calloc(leng_Ll, sizeof(myint));
		this->LlColId = (myint*)calloc(leng_Ll, sizeof(myint));
		this->Llval = (double*)calloc(leng_Ll, sizeof(double));
		myint indj = 0;

		for (indi = 0; indi < leng_v0d1 + leng_v0c + N_edge - bden; indi++) {
			vector<pair<myint, double>> v(Ll[indi].begin(), Ll[indi].end());
			sort(v.begin(), v.end());
			for (auto vi : v) {
				//if (abs(vi.second) > 1e-8) {
				this->LlRowId[indj] = indi;
				this->LlColId[indj] = vi.first;
				this->Llval[indj] = vi.second;
				indj++;
				//}
			}
			v.clear();
		}
		Ll.clear();
	}

	void generateLaplacianRight(double* bm, double* b) {
		/* Lr is the result of [V0a'*b;b]
		Note : b is of the size N_edge - bden */
		int indx, indy, indz;
		double lxa, lya, lza;
		myint eno, i;

		/* Generate V0a'*b */
		for (i = 0; i < v0d1num; ++i) {
			bm[v0d1ColIdo[i]] += v0d1avalo[i] * b[mapEdge[v0d1RowId[i]]] / (v0dan[v0d1ColIdo[i]]);
		}
		for (i = 0; i < v0cnum; ++i) {
			bm[v0cColIdo[i] + leng_v0d1] += v0cavalo[i] * b[mapEdge[v0cRowId[i]]] / (v0can[v0cColIdo[i]]);
		}

		/* Generate b */
		for (i = 0; i < N_edge - bden; ++i) {
			bm[leng_v0d1 + leng_v0c + i] = b[i];
		}

	}

	/* Destructor */
	~fdtdMesh() {
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
		free(this->mapc);
		free(this->mapd);
		free(this->LlRowId);
		free(this->LlColId);
		free(this->Llval);
		free(this->v0dn);
		free(this->v0dan);
		free(this->v0cn);
		free(this->v0can);
		free(this->mapio);
		free(this->mapioR);
		free(this->LooRowId);
		free(this->LooRowId1);
		free(this->LooColId);
		free(this->Looval);
		free(this->mapEdge);
		free(this->mapEdgeR);
	}
};




int meshAndMark(fdtdMesh* sys, unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi, unordered_set<double>* portCoorx, unordered_set<double>* portCoory);
int matrixConstruction(fdtdMesh* sys);
int portSet(fdtdMesh* sys, unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi);
int COO2CSR(vector<int>& rowId, vector<int>& ColId, vector<double>& val);
int paraGenerator(fdtdMesh* sys, unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi);
int COO2CSR_malloc(myint* rowId, myint* ColId, double* val, myint totalnum, myint leng, myint* rowId1);
int generateStiff(fdtdMesh* sys);
int setsideLen(int node, double sideLen, int* markLayerNode, int* markProSide, fdtdMesh* sys);
int generateStiff(fdtdMesh* sys);
int mklMatrixMulti_nt(fdtdMesh* sys, myint& leng_A, myint* aRowId, myint* aColId, double* aval, myint arow, myint acol, myint* bRowId, myint* bColId, double* bval);
int find_Vh_central(fdtdMesh* sys, int sourcePort);
int find_Vh_back(fdtdMesh* sys, int sourcePort, sparse_matrix_t& v0ct, sparse_matrix_t& v0cat, sparse_matrix_t& v0dt, sparse_matrix_t& v0dat, myint* A12RowId, myint* A12ColId, double* A12val, myint leng_A12, myint* A21RowId, myint* A21ColId, double* A21val, myint leng_A21, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, myint* SRowId1, myint* SColId, double* Sval, sparse_matrix_t& Ll);   // findVh.cpp
int matrix_multi(char operation, lapack_complex_double* a, myint arow, myint acol, lapack_complex_double* b, myint brow, myint bcol, lapack_complex_double* tmp3);
int reference(fdtdMesh* sys, int freqNo, myint* RowId, myint* ColId, double* val, complex<double>* xsol);
int plotTime(fdtdMesh* sys, int sourcePort, double* u0d, double* u0c);
myint generateLaplacian_count(fdtdMesh* sys);    // count how many nnz in L, and return the number
int generateLaplacian(fdtdMesh* sys, myint* rowId, myint* colId, double* val);    // generate the Laplacian matrix
void get1122Block_count(myint leng11, myint leng22, fdtdMesh* sys, myint& leng_A12, myint& leng_A21, myint &leng_A22);   // findVh.cpp
void get1122Block(myint leng11, myint leng22, fdtdMesh* sys, myint* A12RowId, myint* A12ColId, double* A12val, myint* A21RowId, myint* A21ColId, double* A21val, myint* A22RowId, myint* A22ColId, double* A22val);   // findVh.cpp
int mkl_gmres(fdtdMesh* sys, double* bm, double* x, sparse_matrix_t Ll, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, sparse_matrix_t v0ct, sparse_matrix_t v0cat, sparse_matrix_t v0dt, sparse_matrix_t v0dat);   // findVh.cpp
int combine_x(double* x, fdtdMesh* sys, double* xr);   // findVh.cpp
int applyPrecond(fdtdMesh* sys, double* b1, double* b2, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, sparse_matrix_t v0ct, sparse_matrix_t v0cat, sparse_matrix_t v0dt, sparse_matrix_t v0dat, int choice);   // findVh.cpp
int solveBackMatrix(fdtdMesh* sys, double* bm, double* x, sparse_matrix_t& v0ct, sparse_matrix_t& v0cat, sparse_matrix_t& v0dt, sparse_matrix_t& v0dat, myint* A12RowId, myint* A12ColId, double* A12val, myint leng_A12, myint* A21RowId, myint* A21ColId, double* A21val, myint leng_A21, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22);    // findVh.cpp
int solveA11Matrix(fdtdMesh* sys, double* rhs, sparse_matrix_t& v0ct, sparse_matrix_t& v0cat, sparse_matrix_t& v0dt, sparse_matrix_t& v0dat, double* sol);    // findVh.cpp
int sparseMatrixVecMul(myint* rowId, myint* colId, double* val, myint leng, double* v1, double* v2);   // findVh.cpp
int sparseMatrixVecMul_c(myint* rowId, myint* colId, double* val, myint leng, complex<double>* v1, complex<double>* v2);   // findVh.cpp
int pardisoSolve(myint* rowId, myint* colId, double* val, double* rsc, double* xsol, myint size);   // findVh.cpp
int storeTimeRespValue(fdtdMesh* sys, double** resp, int ind, double* xr);   // findVh.cpp
int mklFFT(fdtdMesh* sys, double* time, complex<double>* freq, int N);   // findVh.cpp
																		 //void sparseMatrixSum(fdtdMesh* sys, sparse_matrix_t& A, myint* browId1, myint* bcolId, double* bval, myint Rows);   // matrixCon.cpp
void sparseMatrixSum(fdtdMesh* sys, sparse_matrix_t& A, myint* browId1, myint* bcolId, double* bval, myint Rows, myint** rRowId, myint** rColId, double** rval, myint& lengr); // matrixCon.cpp
void sparseMatrixSum1(fdtdMesh* sys, myint* arowId1, myint* acolId, double* aval, myint* browId1, myint* bcolId, double* bval, myint Rows);  // findVh.cpp
void matrixOutside(fdtdMesh* sys, myint* ArowStart, myint* ArowEnd, myint* AcolId, double* Aval, myint ARows, myint* AoorowId, myint* AoocolId, double* Aooval, double scale);   // matrixCon.cpp
void matrixOutside_count(fdtdMesh* sys, myint* ArowStart, myint* ArowEnd, myint* AcolId, double* Aval, myint ARows, myint& leng_Aoo);   // matrixCon.cpp
void matrixInsideOutside(myint* mapio, myint outside, myint inside, myint* rowId, myint* colId, double* val, myint leng, myint* ooRowId, myint* ooColId, double* ooval, myint* oiRowId, myint* oiColId, double* oival, myint* ioRowId, myint* ioColId, double* ioval);   // matrixCon.cpp
void matrixInsideOutside_count(myint* mapio, myint outside, myint inside, myint* rowId, myint* colId, myint leng, myint& lengoo, myint& lengoi, myint& lengio);   // matrixCon.cpp
int mkl_gmres_A(double* bm, double* x, myint* ARowId, myint* AColId, double* Aval, myint leng_A, myint N);   // findVh.cpp
int solveFreqIO(fdtdMesh* sys, int freqi, complex<double>* y, myint* MooRowId, myint* MooRowId1, myint* MooColId, double* Mooval, myint lengoo, myint* SioRowId, myint* SioColId, double* Sioval, myint lengio, myint* PRowId1, myint* PColId, double* Pval, myint leng_P);   // findVh.cpp
void solveOO(fdtdMesh* sys, double freq, double* Jo, sparse_matrix_t& v0dt, sparse_matrix_t& v0dat, double* xo);  // findVh.cpp
void solveOO_pardiso(fdtdMesh* sys, myint* MooRowId1, myint* MooColId, double* Mooval, myint leng, myint N, double* Jo, double* xo, int rhs_s);   // findVh.cpp
int checkCSRRepeats(myint* RowIdStart, myint* RowIdEnd, myint* ColId, double* val, myint& leng, myint N);   // matrixCon.cpp
int reference_oo(fdtdMesh *sys, int freqNo, myint *RowId1, myint *ColId, double *val, double* xsol); // generateStiff.cpp
int mkl_gmres_A_P(double* bm, double* x, myint* ARowId, myint* AColId, double* Aval, myint leng_A, myint N, myint* PRowId1, myint* PColId, double* Pval, myint   leng_P);   // findVh.cpp
int solveV0L(fdtdMesh* sys, int freqNo, myint* LrowId, myint* LcolId, double* Lval, myint leng, myint N, sparse_matrix_t& V0dt, sparse_matrix_t& V0dat, sparse_matrix_t& V0ct, sparse_matrix_t& V0cat);   // findVh.cpp
void generateL1(fdtdMesh* sys, double freq, myint* LRowId, myint* LColId, double* Lval, myint leng_L, myint N, myint* L1RowId, myint* L1ColId, double* L1val, myint& leng_L1);   // matrixCon.cpp
#endif