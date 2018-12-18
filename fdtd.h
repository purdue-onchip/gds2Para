#ifndef _FDTD_H
#define _FDTD_H
//#include "stdafx.h"
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <set>
#include <complex>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <utility>
//#include <mkl>


using namespace std;

#define PI (3.1415926)
#define MU (4*PI*1.e-7)
#define CSPED (299792458.)
#define EPSILON0 (8.854e-12)
#define SIGMA (5.7e+7)
#define FDTD_MAXC (256*6)
#define STACKNUM (20)
//#define length (128)
#define DOUBLEMAX (1.e+30)
#define DOUBLEMIN (-1.e-30)


class fdtdOneCondct {
public:
    int cdtName;
    int vcc;
    int plaStackId;
    int portId;
    int numNode;
    int numVert;
    double *x;
    double *y;
    double *z;
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double sigma;
    int *cdtInNode;
    int *markCdtInNode;
};

typedef class {
public:
    int *node;
    int markPort;
    int *portNode;
    int portind;
    int cdtNodeind;
} fdtdCdt;

typedef class {
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
} fdtdBound;

typedef class {
public:
    double *x;
    double *y;
    double *z;
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double x1, x2;
    double y1, y2;
    double z1, z2;
    int portDirection;
} fdtdPort;

typedef class{
public:
    double patchArea;
    int *node;
}fdtdPatch;


class fdtdMesh {
    /*mesh information*/
public: 
    double lengthUnit;
    /* frequency parameter */
    double freqUnit;
    double freqStart;
    double freqEnd;
    int nfreq;
    int freqScale;


    int ix;    //x direction number of periods
    int iy;    //y direction number of periods

    int nx, ny, nz;    //number of nodes along x, y, z

    double *xn, *yn, *zn;    //coordinates of the nodes along x, y, z


    int N_cell_x;
    int N_cell_y;
    int N_cell_z;

    double factor_x;
    double factor_y;
    double factor_z;
    double xlim1, xlim2;
    double ylim1, ylim2;
    double zlim1, zlim2;

    double txmin, txmax;
    double tymin, tymax;
    double tzmin, tzmax;
    
    int N_edge;
    int N_edge_s;
    int N_edge_v;

    int N_node;
    int N_node_s;

    int N_patch;
    int N_patch_s;
    int N_patch_v;

    double *nodepos;
    double *Epoints;
    int *edgelink;
    double *Hpoints;
    vector<vector<pair<int,double> > > nodeEdge;    // for each node which edge is connected with this node

    

    /*upper and lower PEC*/
    int *bd_node1;   //lower PEC
    int *bd_node2;   //upper PEC
    int *bd_edge;

    /*material parameter*/
    int numStack;
    double *stackEps;
    double *stackBegCoor;
    double *stackEndCoor;
    vector<string> stackName;
    double    *eps;

    /*conductor parameter*/
    fdtdOneCondct *conductorIn;
    int numCdtRow;   //how many input rows
    int numCdt;
    int *markEdge;    //mark if this edge is inside a conductor
    int *cdtNumNode;
    double *sig;
    fdtdCdt *conductor;
    int *markNode;    // mark this node if it is inside the conductor
    

    /* patch information */
    fdtdPatch *patch;

    /*boundary*/
    fdtdBound *bound;

    /*V0c row, column, value*/
    vector<int> v0cRowId;
    vector<int> v0cColId;
    vector<int> v0cColIdo;
    vector<double> v0cval;
    vector<double> v0cvalo;
    vector<int> v0c2RowId;
    vector<int> v0c2ColId;
    vector<int> v0c2ColIdo;
    vector<double> v0c2val;
    vector<int> y0c2RowId;
    vector<int> y0c2ColId;
    vector<double> y0c2val;
    double *v0c2y0c2;
    double *v0c2y0c2o;
    double *yc;
    double *v0cy0c;

    /* V0c'*D_sig*V0c row, column, value*/
    vector<int> AcRowId;
    vector<int> AcColId;
    vector<double> Acval;
    vector<int> crhsRowId;
    vector<int> crhsColId;
    vector<double> crhsval;
    double *crhs;


    /* V0d row, column, value*/
    vector<int> v0d1RowId;
    vector<int> v0d1ColId;
    vector<double> v0d1val;
    vector<double> v0d1valo;
    vector<int> v0d2RowId;
    vector<int> v0d2ColId;
    vector<int> v0d2ColIdo;
    vector<double> v0d2val;
    vector<double> v0d2valo;
    double *yd_re, *yd_im;


    /* final solution */
    complex<double> *y;
    complex<double> *x;    // the solution evolving all the sourcePorts


    /*port coordinates*/
    int numPorts;
    fdtdPort *portCoor;
    vector<vector<int> > portEdge;
    vector<double> portArea;


    /*current source*/
    double *J;
    

    /* current V0c,s^T*I matrix */
    complex<double> *v0csJ;
    complex<double> *Y;
    complex<double> *voltage;
    
    fdtdMesh(){
        numCdtRow = 0;
    }
    ~fdtdMesh(){
        numCdtRow = 0;
    }
};




int readInput(fdtdMesh* sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi);
int parameterConstruction(fdtdMesh* sys, unordered_map<double,int> xi, unordered_map<double,int> yi, unordered_map<double,int> zi);
bool polyIn(double x, double y, fdtdMesh *sys, int inPoly);
int fdtdStringWord(char*, char *word[]);
double fdtdGetValue(char *str);
void freePara(fdtdMesh *sys);
int matrixConstruction(fdtdMesh *sys);
int portSet(fdtdMesh *sys, unordered_map<double,int> xi, unordered_map<double,int> yi, unordered_map<double,int> zi);
int matrixMulti(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> bRowId, vector<int> bColId, vector<double> bval, vector<int>& cRowId, vector<int>& cColId, vector<double>& cval);
// the first is read row by row, and the second one is read column by column
int COO2CSR(vector<int>& rowId, vector<int>& ColId, vector<double>& val);
int mvMulti(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int>& bRowId, vector<int>& bColId, vector<double>& bval, double *index_val, int size);
int nodeAdd(vector<int> &rowId, vector<int> &colId, vector<double> &val, vector<int> index, fdtdMesh *sys);
int interativeSolver(int N, int nrhs, double *rhs, int *ia, int *ja, double *a, double *solution, fdtdMesh *sys);
int output(fdtdMesh *sys);
int paraGenerator(fdtdMesh *sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi);
int yParaGenerator(fdtdMesh *sys);
#endif
