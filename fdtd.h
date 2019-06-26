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
#include <iostream>
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

 //HYPRE and MKL data type control
#define LARGE_SYSTEM (1)
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
#define SIGMA (5.8e+7)
#define FDTD_MAXC (256*6)
#define STACKNUM (20)
//#define SOLVERLENGTH (128)
#define DOUBLEMAX (1.e+30)
#define DOUBLEMIN (-1.e+30)
#define MINDISFRACXY (1.0e-6) // Fraction setting minimum discretization retained in x- or y-directions after node merging in terms of smaller of x-extent or y-extent
#define MINDISFRACZ (0.1) // Fraction setting minimum discretization retained in z-direction after node merging in terms of distance between closest layers
#define MAXDISFRACX (0.1) // Fraction setting largest discretization in x-direction in terms of x-extent
#define MAXDISFRACY (MAXDISFRACX) // Fraction setting largest discretization in y-direction in terms of y-extent
#define MAXDISLAYERZ (1.) // Largest discretization in z-direction represented as fewest nodes placed between closest layers (1. = distance between closest layers, 2. = half distance between closest layers)
#define DT (1.e-16)

// HYPRE control macros
#define HYPRE_METHOD (3) // 1 = AMG, 2 = PCG with AMG Preconditioner, 3 = Flexible GMRES with AMG Preconditioner
#define HYPRE_CONV_TOL (1.e-4) // Convergence relative tolerance for HYPRE
#define HYPRE_MAX_ITER (100) // Maximum iterations for HYPRE

// Debug testing macros (comment out if not necessary)
#define PRINT_NODE_COORD
#define PRINT_DIS_COUNT (1)
#define SKIP_MARK_CELL
//#define PRINT_PORT_SET
//#define PRINT_V0D_BLOCKS
#define SKIP_PARDISO // Remove PARDISO solver code
#define SKIP_VH // Turn on to save a lot of time
#define SKIP_STIFF_REFERENCE 


// Function-like macros
#define NELEMENT(x) (sizeof(x) / sizeof((x)[0]))


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
    int cdtNodeind;
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
    double *x;
    double *y;
    double *z;
    double x1, x2;
    double y1, y2;
    double z1, z2;
    int portDirection;
    int *node;
    int portCnd;

    /* Default Constructor */
    fdtdPort(){
        x = NULL;
        y = NULL;
        z = NULL;
        portDirection = 0;
        node = NULL;
        portCnd = 0;
    }

    /* Print Function */
    void print();

    /* Destructor */
    ~fdtdPort(){
        portDirection = 0;
        portCnd = 0;
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

    double *nodepos;
    double *Epoints;
    myint *edgelink;
    double *Hpoints;
    vector<vector<pair<myint,double>>> nodeEdge;    // for each node which edge is connected with this node
    vector<vector<pair<myint, double>>> nodeEdgea;    // for each node which edge is connected with this node

    /* Upper and lower PEC */
    int *bd_node1;   //lower PEC
    int *bd_node2;   //upper PEC
    int *bd_edge;

    /* Layer stack up parameters */
    int numStack;
    vector<double> stackEps;
    vector<double> stackBegCoor;
    vector<double> stackEndCoor;
    vector<string> stackName;
    double *eps;
    vector<double> stackEpsn;
    double *stackCdtMark;

    /* Conductor parameters */
    vector<fdtdOneCondct> conductorIn;
    myint numCdtRow;    // how many input rows
    myint numCdt;       // number of isolated conductors in design
    int *markEdge;    // mark if this edge is inside a conductor
    int *markCell;
    myint *cdtNumNode;
    double *sig;
    fdtdCdt *conductor;
    int *markNode;    // mark this node if it is inside the conductor
    vector<vector<int>> edgeCell;    // for each cell which edge is around it
    vector<vector<double>> edgeCellArea;    // for each cell the area of the perpendicular rectangle
    vector<int> acu_cnno; // accumulated conductor number of nodes
    vector<int> cindex;
    int *exciteCdtLayer;
    vector<unordered_set<int>> cond2condIn;
    int *markProSide;

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
    myint *SRowId, *SColId;
    double *Sval;
    myint leng_S;

    /* Solution storage */
    complex<double> *y;
    complex<double> *x;    // the solution involving all the sourcePorts

    /* Port coordinates */
    int numPorts;
    fdtdPort *portCoor;
    vector<vector<int>> portEdge;
    vector<double> portArea;

    /* Current sources */
    double *J;

    /* Current V0c,s^T*I matrix */
    complex<double> *v0csJ;
    complex<double> *Y;

    /* Default Constructor */
    fdtdMesh(){
        numCdtRow = (myint)0;
        leng_S = (myint)0;

        // Set all pointers to NULL for safety
        xn = NULL;
        yn = NULL;
        zn = NULL;
        xnu = NULL;
        ynu = NULL;
        znu = NULL;
        nodepos = NULL;
        Epoints = NULL;
        edgelink = NULL;
        Hpoints = NULL;
        bd_node1 = NULL;
        bd_node2 = NULL;
        bd_edge = NULL;
        eps = NULL;
        stackCdtMark = NULL;
        markEdge = NULL;
        markCell = NULL;
        cdtNumNode = NULL;
        sig = NULL;
        conductor = NULL;
        markNode = NULL;
        exciteCdtLayer = NULL;
        markProSide = NULL;
        patch = NULL;
        bound = NULL;
        v0cRowId = NULL;
        v0cColId = NULL;
        v0cColIdo = NULL;
        v0cval = NULL;
        v0cvalo = NULL;
        v0caRowId = NULL;
        v0caColId = NULL;
        v0caColIdo = NULL;
        v0caval = NULL;
        v0cavalo = NULL;
        v0c2y0c2 = NULL;
        v0c2y0c2o = NULL;
        yc = NULL;
        v0cy0c = NULL;
        AcRowId = NULL;
        AcRowId1 = NULL;
        AcColId = NULL;
        Acval = NULL;
        AdRowId = NULL;
        AdRowId1 = NULL;
        AdColId = NULL;
        Adval = NULL;
        crhs = NULL;
        v0d1RowId = NULL;
        v0d1ColId = NULL;
        v0d1ColIdo = NULL;
        v0d1val = NULL;
        v0d1valo = NULL;
        v0d1aRowId = NULL;
        v0d1aColId = NULL;
        v0d1aColIdo = NULL;
        v0d1aval = NULL;
        v0d1avalo = NULL;
        v0d2RowId = NULL;
        v0d2ColId = NULL;
        v0d2ColIdo = NULL;
        v0d2val = NULL;
        v0d2valo = NULL;
        v0d2aRowId = NULL;
        v0d2aColId = NULL;
        v0d2aColIdo = NULL;
        v0d2aval = NULL;
        v0d2avalo = NULL;
        yd = NULL;
        Vh = NULL;
        SRowId = NULL;
        SColId = NULL;
        Sval = NULL;
        y = NULL;
        x = NULL;
        portCoor = NULL;
        J = NULL;
        v0csJ = NULL;
        Y = NULL;
    }

    /* Print Function */
    void print();

    /* Destructor */
    ~fdtdMesh(){
        numCdtRow = (myint)0;
        leng_S = (myint)0;
    }
};




int meshAndMark(fdtdMesh* sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi, unordered_set<double> *portCoorx, unordered_set<double> *portCoory);
int compute_edgelink(fdtdMesh *sys, myint eno, myint &node1, myint &node2);
int parameterConstruction(fdtdMesh* sys, unordered_map<double,int> xi, unordered_map<double,int> yi, unordered_map<double,int> zi);
bool polyIn(double x, double y, fdtdMesh *sys, int inPoly);
int fdtdStringWord(char*, char *word[]);
double fdtdGetValue(char *str);
int compareString(char *a, char *b);
//int nodeAddAvg(vector<int> &rowId, vector<int> &colId, vector<double> &val, int index, int size, fdtdMesh *sys);
void freePara(fdtdMesh *sys);
int matrixConstruction(fdtdMesh *sys);
int portSet(fdtdMesh *sys, unordered_map<double,int> xi, unordered_map<double,int> yi, unordered_map<double,int> zi);
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
int pardisoSolve_c(fdtdMesh *sys, double *rhs, double *solution, int nodestart, int nodeend, int indstart, int indend);
int pardisoSolve_r(fdtdMesh *sys, complex<double> *rhs, int *RowId, int *ColId, complex<double> *val, int nnz, int size, complex<double> *solution);
int COO2CSR_malloc(myint *rowId, myint *ColId, double *val, myint totalnum, myint leng, myint *rowId1);
int generateStiff(fdtdMesh *sys);
int merge_v0d1(fdtdMesh *sys, double block1_x, double block1_y, double block2_x, double block2_y, double block3_x, double block3_y, myint &v0d1num, myint &leng_v0d1, myint &v0d1anum, myint &leng_v0d1a, myint *map, double sideLen);
int merge_v0c(fdtdMesh *sys, double block_x, double block_y, double block2_x, double block2_y, myint &v0cnum, myint &leng_v0c, myint &v0canum, myint &leng_v0ca, myint *map);
int setsideLen(int node, double sideLen, int *markLayerNode, int *markProSide, fdtdMesh *sys);
int generateStiff(fdtdMesh *sys);
int mklMatrixMulti_nt(fdtdMesh *sys, myint &leng_A, myint *aRowId, myint *aColId, double *aval, myint arow, myint acol, myint *bRowId, myint *bColId, double *bval);
int find_Vh(fdtdMesh *sys, lapack_complex_double *u0, lapack_complex_double *u0a, int sourcePort);
int matrix_multi(char operation, lapack_complex_double *a, myint arow, myint acol, lapack_complex_double *b, myint brow, myint bcol, lapack_complex_double *tmp3);
int reference(fdtdMesh *sys, lapack_complex_double *x, myint *RowId, myint *ColId, double *val);
int plotTime(fdtdMesh *sys, int sourcePort, double *u0d, double *u0c);
int avg_length(fdtdMesh *sys, int iz, int iy, int ix, double &lx, double &ly, double &lz);
#endif
