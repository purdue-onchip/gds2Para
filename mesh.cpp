//#include "stdafx.h"
#include "fdtd.h"

static bool comp(pair<double, int> a, pair<double, int> b)
{
    return a.second <= b.second;
};

int meshAndMark(fdtdMesh *sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi, unordered_set<double> *portCoorx, unordered_set<double> *portCoory)
{
    int lyr;
    int i, j, k, m;
    int status;
    int count;
    int mark;
    double xmin, xmax, ymin, ymax, xwid, ywid;

    /*cout << " Print the conductor information: " << endl;
    for (i = 0; i < sys->numCdtRow; i++){
        cout << sys->conductorIn[i].xmin << " " << sys->conductorIn[i].xmax << " " << sys->conductorIn[i].ymin << " " << sys->conductorIn[i].ymax << " " << sys->conductorIn[i].zmin << " " << sys->conductorIn[i].zmax << endl;
    }*/

    /* Generate the mesh nodes based on conductorIn information */
    int numNode = 0;
    double *xOrigOld, *yOrigOld, *zOrigOld;
    double disMin = 1.e-9;
    double disMaxx, disMaxy;   // the max discretization in x, y, z directions
    for (i = 0; i < sys->numCdtRow; i++){
        numNode += sys->conductorIn[i].numVert;
    }
    //cout << numNode << endl;
    xOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    yOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    zOrigOld = (double*)calloc(2 * sys->numStack + 2 * sys->numPorts, sizeof(double));

    j = 0;
    for (i = 0; i < sys->numCdtRow; i++){
        for (k = 0; k < sys->conductorIn[i].numVert; k++){
            xOrigOld[j] = sys->conductorIn[i].x[k];
            yOrigOld[j] = sys->conductorIn[i].y[k];
            j++;
        }
    }

    for (i = 0; i < sys->numPorts; i++){
        xOrigOld[j] = sys->portCoor[i].x1;
        yOrigOld[j] = sys->portCoor[i].y1;
        j++;
        xOrigOld[j] = sys->portCoor[i].x2;
        yOrigOld[j] = sys->portCoor[i].y2;
        j++;
    }
    j = 0;
    for (i = 0; i < sys->numStack; i++){
        zOrigOld[j] = sys->stackBegCoor[i];
        j++;
        zOrigOld[j] = sys->stackEndCoor[i];
        j++;
    }
    for (i = 0; i < sys->numPorts; i++){
        zOrigOld[j] = sys->portCoor[i].z1;
        j++;
        zOrigOld[j] = sys->portCoor[i].z2;
        j++;
    }

    /*******************************************************************************************/
    sort(xOrigOld, xOrigOld + numNode + 2 * sys->numPorts);
    sys->nx = 1;
    xmin = xOrigOld[0];
    xmax = xOrigOld[numNode + 2 * sys->numPorts - 1];
    disMaxx = (xmax - xmin) / 5;

    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(xOrigOld[i] - xOrigOld[i - 1]) > disMin){
            sys->nx++;
        }
    }
    double *xn = (double*)calloc(numNode + 6 * sys->numPorts + 10, sizeof(double));
    xn[0] = xOrigOld[0];
    j = 0;
    sys->nx = 1;
    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(xOrigOld[i] - xOrigOld[i - 1]) > disMin && abs(xOrigOld[i] - xOrigOld[i - 1]) <= disMaxx){
            j++;
            xn[j] = xOrigOld[i];
            sys->nx++;
        }
        else if (abs(xOrigOld[i] - xOrigOld[i - 1]) > disMin && abs(xOrigOld[i] - xOrigOld[i - 1]) > disMaxx){
            while (abs(xOrigOld[i] - xn[j]) > disMaxx){
                j++;
                xn[j] = xn[j - 1] + disMaxx;
                sys->nx++;
            }
            if (abs(xOrigOld[i] - xn[j]) > disMin){
                sys->nx++;
            }
            j++;
            xn[j] = xOrigOld[i];
        }
        else{
            j++;
            xn[j] = xOrigOld[i];
        }
    }
    int countx = j;
    //sort(xn, xn + countx + 1);


    sys->xnu = (double*)calloc(sys->nx, sizeof(double));

    j = 0;
    sys->xnu[0] = xn[0];
    xi[sys->xnu[0]] = j;
    double first, second;
    for (i = 1; i <= countx; i++){    // set the discretization length around port to be equal
        if (abs(xn[i] - xn[i - 1]) > disMin){
            j++;
            sys->xnu[j] = xn[i];
            xi[sys->xnu[j]] = j;
        }
        else{
            xi[xn[i]] = j;
        }
    }
    sys->nx = j + 1;

    /*vector<pair<double, int> > v(xi.begin(), xi.end());
    sort(v.begin(), v.end(), comp);
    for (i = 0; i < v.size(); i++){
    cout << v[i].first << " " << v[i].second << endl;
    }
    cout << sys->nx << endl;*/



    /***************************************************************************/
    sort(yOrigOld, yOrigOld + numNode + 2 * sys->numPorts);

    sys->ny = 1;
    ymin = yOrigOld[0];
    ymax = yOrigOld[numNode + 2 * sys->numPorts - 1];
    disMaxy = (ymax - ymin) / 10;


    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(yOrigOld[i] - yOrigOld[i - 1]) > disMin){
            sys->ny++;
        }
    }
    double *yn = (double*)calloc(numNode + 6 * sys->numPorts + 10, sizeof(double));
    yn[0] = yOrigOld[0];
    j = 0;
    sys->ny = 1;
    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(yOrigOld[i] - yOrigOld[i - 1]) > disMin && abs(yOrigOld[i] - yOrigOld[i - 1]) <= disMaxy){
            j++;
            yn[j] = yOrigOld[i];
            sys->ny++;
        }
        else if (abs(yOrigOld[i] - yOrigOld[i - 1]) > disMin && abs(yOrigOld[i] - yOrigOld[i - 1]) > disMaxy){
            while (abs(yOrigOld[i] - yn[j]) > disMaxy){
                j++;
                yn[j] = yn[j - 1] + disMaxy;
                sys->ny++;
            }
            if (abs(yOrigOld[i] - yn[j]) > disMin){
                sys->ny++;
            }
            j++;
            yn[j] = yOrigOld[i];
        }
        else{
            j++;
            yn[j] = yOrigOld[i];
        }
    }

    int county = j;
    //sort(yn, yn + county + 1);

    sys->ynu = (double*)calloc(sys->ny, sizeof(double));


    j = 0;
    sys->ynu[0] = yn[0];
    yi[sys->ynu[0]] = j;

    for (i = 1; i <= county; i++){    // set the discretization length around port to be equal
        if (abs(yn[i] - yn[i - 1]) > disMin){
            j++;
            sys->ynu[j] = yn[i];
            yi[sys->ynu[j]] = j;
        }
        else{
            yi[yn[i]] = j;
        }
    }
    sys->ny = j + 1;
    /*vector<pair<double, int> > v(yi.begin(), yi.end());
    sort(v.begin(), v.end(), comp);
    for (i = 0; i < v.size(); i++){
    cout << v[i].first << " " << v[i].second << endl;
    }
    cout << sys->ny << endl;*/



    /********************************************************************************/
    sort(zOrigOld, zOrigOld + 2 * sys->numStack + 2 * sys->numPorts);
    sys->nz = 1;
    for (i = 1; i < 2 * sys->numStack + 2 * sys->numPorts; i++){
        if (abs(zOrigOld[i] - zOrigOld[i - 1]) > disMin){
            sys->nz++;
        }
    }

    double *zn = (double*)calloc(2 * sys->numStack + 6 * sys->numPorts, sizeof(double));
    for (i = 0; i < 2 * sys->numStack + 2 * sys->numPorts; i++){
        zn[i] = zOrigOld[i];
    }
    int countz = 2 * sys->numStack + 2 * sys->numPorts - 1;
    //sort(zn, zn + countz + 1);

    sys->znu = (double*)calloc(sys->nz, sizeof(double));
    sys->znu[0] = zOrigOld[0];
    zi[sys->znu[0]] = 0;
    j = 0;
    for (i = 1; i < 2 * sys->numStack + 2 * sys->numPorts; i++){
        if (abs(zOrigOld[i] - zOrigOld[i - 1]) > disMin){
            j++;
            sys->znu[j] = zOrigOld[i];
            zi[sys->znu[j]] = j;
        }
        else{
            zi[zOrigOld[i]] = j;    // considering the double precision, put all the show up y in the map
        }
    }


    /*************************************************************************************/
    /* make the edge around the port be the same */
    vector<double> porteqx;    // store the eq leng caused x coordinates around each port
    vector<double> porteqy;    // store the eq leng caused y coordinates around each port
    vector<double> porteqz;    // store the eq leng caused z coordinates around each port
    /*double mini;
    for (i = 0; i < sys->numPorts; i++){
        if (sys->portCoor[i].x1 == sys->portCoor[i].x2){
            mini = DOUBLEMAX;
            if (sys->portCoor[i].x1 - sys->xnu[xi[sys->portCoor[i].x1] - 1] < mini){
                mini = sys->portCoor[i].x1 - sys->xnu[xi[sys->portCoor[i].x1] - 1];
            }
            if (sys->xnu[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1 < mini){
                mini = sys->xnu[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1;
            }
            if (sys->portCoor[i].y1 - sys->ynu[yi[sys->portCoor[i].y1] - 1] < mini){
                mini = sys->portCoor[i].y1 - sys->ynu[yi[sys->portCoor[i].y1] - 1];
            }
            if (sys->portCoor[i].z1 - sys->znu[zi[sys->portCoor[i].z1] - 1] < mini){
                mini = sys->portCoor[i].z1 - sys->znu[zi[sys->portCoor[i].z1] - 1];
            }
            if (sys->ynu[yi[sys->portCoor[i].y2] + 1] - sys->portCoor[i].y2 < mini){
                mini = sys->ynu[yi[sys->portCoor[i].y2] + 1] - sys->portCoor[i].y2;
            }
            if (sys->znu[zi[sys->portCoor[i].z2] + 1] - sys->portCoor[i].z2 < mini){
                mini = sys->znu[zi[sys->portCoor[i].z2] + 1] - sys->portCoor[i].z2;
            }
            porteqx.push_back(sys->portCoor[i].x1 - mini);
            porteqx.push_back(sys->portCoor[i].x1 + mini);
            porteqy.push_back(sys->portCoor[i].y1 - mini);
            porteqy.push_back(sys->portCoor[i].y2 + mini);
            porteqz.push_back(sys->portCoor[i].z1 - mini);
            porteqz.push_back(sys->portCoor[i].z2 + mini);

            porteqy.push_back(sys->portCoor[i].y1 + mini);
            porteqy.push_back(sys->portCoor[i].y2 - mini);
            porteqz.push_back(sys->portCoor[i].z1 + mini);
            porteqz.push_back(sys->portCoor[i].z2 - mini);
        }
        else if (sys->portCoor[i].y1 == sys->portCoor[i].y2){
            mini = DOUBLEMAX;
            if (sys->portCoor[i].y1 - sys->ynu[yi[sys->portCoor[i].y1] - 1] < mini){
                mini = sys->portCoor[i].y1 - sys->ynu[yi[sys->portCoor[i].y1] - 1];
            }
            if (sys->ynu[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1 < mini){
                mini = sys->ynu[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1;
            }
            if (sys->portCoor[i].x1 - sys->xnu[xi[sys->portCoor[i].x1] - 1] < mini){
                mini = sys->portCoor[i].x1 - sys->xnu[xi[sys->portCoor[i].x1] - 1];
            }
            if (sys->portCoor[i].z1 - sys->znu[zi[sys->portCoor[i].z1] - 1] < mini){
                mini = sys->portCoor[i].z1 - sys->znu[zi[sys->portCoor[i].z1] - 1];
            }
            if (sys->xnu[xi[sys->portCoor[i].x2] + 1] - sys->portCoor[i].x2 < mini){
                mini = sys->xnu[xi[sys->portCoor[i].x2] + 1] - sys->portCoor[i].x2;
            }
            if (sys->znu[zi[sys->portCoor[i].z2] + 1] - sys->portCoor[i].z2 < mini){
                mini = sys->znu[zi[sys->portCoor[i].z2] + 1] - sys->portCoor[i].z2;
            }
            porteqx.push_back(sys->portCoor[i].x1 - mini);
            porteqx.push_back(sys->portCoor[i].x2 + mini);
            porteqy.push_back(sys->portCoor[i].y1 - mini);
            porteqy.push_back(sys->portCoor[i].y1 + mini);
            porteqz.push_back(sys->portCoor[i].z1 - mini);
            porteqz.push_back(sys->portCoor[i].z2 + mini);

            porteqx.push_back(sys->portCoor[i].x1 + mini);
            porteqx.push_back(sys->portCoor[i].x2 - mini);
            porteqz.push_back(sys->portCoor[i].z1 + mini);
            porteqz.push_back(sys->portCoor[i].z2 - mini);
        }
        else if (sys->portCoor[i].z1 == sys->portCoor[i].z2){
            mini = DOUBLEMAX;
            if (sys->portCoor[i].z1 - sys->znu[zi[sys->portCoor[i].z1] - 1] < mini){
                mini = sys->portCoor[i].z1 - sys->znu[zi[sys->portCoor[i].z1] - 1];
            }
            if (sys->znu[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1 < mini){
                mini = sys->znu[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1;
            }
            if (sys->portCoor[i].x1 - sys->xnu[xi[sys->portCoor[i].x1] - 1] < mini){
                mini = sys->portCoor[i].x1 - sys->xnu[xi[sys->portCoor[i].x1] - 1];
            }
            if (sys->portCoor[i].y1 - sys->ynu[yi[sys->portCoor[i].y1] - 1] < mini){
                mini = sys->portCoor[i].y1 - sys->ynu[yi[sys->portCoor[i].y1] - 1];
            }
            if (sys->xnu[xi[sys->portCoor[i].x2] + 1] - sys->portCoor[i].x2 < mini){
                mini = sys->xnu[xi[sys->portCoor[i].x2] + 1] - sys->portCoor[i].x2;
            }
            if (sys->ynu[yi[sys->portCoor[i].y2] + 1] - sys->portCoor[i].y2 < mini){
                mini = sys->ynu[yi[sys->portCoor[i].y2] + 1] - sys->portCoor[i].y2;
            }
            porteqx.push_back(sys->portCoor[i].x1 - mini);
            porteqx.push_back(sys->portCoor[i].x2 + mini);
            porteqy.push_back(sys->portCoor[i].y1 - mini);
            porteqy.push_back(sys->portCoor[i].y2 + mini);
            porteqz.push_back(sys->portCoor[i].z1 - mini);
            porteqz.push_back(sys->portCoor[i].z1 + mini);

            porteqx.push_back(sys->portCoor[i].x1 + mini);
            porteqx.push_back(sys->portCoor[i].x2 - mini);
            porteqy.push_back(sys->portCoor[i].y1 + mini);
            porteqy.push_back(sys->portCoor[i].y2 - mini);
        }

    }
    for (i = 0; i < porteqx.size(); i++){
        xn[countx + 1 + i] = porteqx[i];
    }
    for (i = 0; i < porteqy.size(); i++){
        yn[county + 1 + i] = porteqy[i];
    }
    for (i = 0; i < porteqz.size(); i++){
        zn[countz + 1 + i] = porteqz[i];
    }*/



    sort(xn, xn + countx + 1 + porteqx.size());
    xi.clear();
    sys->xn = (double*)calloc(sys->nx + 4 * sys->numPorts, sizeof(double));
    j = 0;
    sys->xn[0] = xn[0];
    xi[sys->xn[0]] = j;
    for (i = 1; i <= countx + porteqx.size(); i++){    // set the discretization length around port to be equal
        if (abs(xn[i] - xn[i - 1]) > disMin){
            j++;
            sys->xn[j] = xn[i];
            xi[sys->xn[j]] = j;
        }
        else{
            xi[xn[i]] = j;
        }
    }
    sys->nx = j + 1;


    sort(yn, yn + county + 1 + porteqy.size());
    yi.clear();
    sys->yn = (double*)calloc(sys->ny + 4 * sys->numPorts, sizeof(double));
    j = 0;
    sys->yn[0] = yn[0];
    yi[sys->yn[0]] = j;
    for (i = 1; i <= county + porteqy.size(); i++){    // set the discretization length around port to be equal
        if (abs(yn[i] - yn[i - 1]) > disMin){
            j++;
            sys->yn[j] = yn[i];
            yi[sys->yn[j]] = j;
        }
        else{
            yi[yn[i]] = j;
        }
    }
    sys->ny = j + 1;

    sort(zn, zn + countz + 1 + porteqz.size());
    cout << endl;
    zi.clear();
    sys->zn = (double*)calloc(sys->nz + 4 * sys->numPorts, sizeof(double));
    j = 0;
    sys->zn[0] = zn[0];
    zi[sys->zn[0]] = j;
    for (i = 1; i <= countz + porteqz.size(); i++){    // set the discretization length around port to be equal
        if (abs(zn[i] - zn[i - 1]) > disMin){
            j++;
            sys->zn[j] = zn[i];
            zi[sys->zn[j]] = j;
        }
        else{
            zi[zn[i]] = j;
        }
    }
    sys->nz = j + 1;

    sys->stackEpsn = (double*)calloc(sys->nz - 1, sizeof(double));
    i = 0;
    j = 0;
    while (i < sys->nz - 1){
        if ((sys->zn[i] + sys->zn[i + 1]) / 2 >= sys->stackBegCoor[j] && (sys->zn[i] + sys->zn[i + 1]) / 2 <= sys->stackEndCoor[j]){
            sys->stackEpsn[i] = sys->stackEps[j];
            i++;
        }
        else{
            j++;
        }
    }
    /*for (i = 0; i < sys->nz - 1; i++){
        cout << sys->stackEpsn[i] << endl;
    }*/

    /*for (i = 0; i < sys->nx; i++){
        cout << sys->xn[i] << " ";
    }
    cout << "\n" << endl;
    for (i = 0; i < sys->ny; i++){
        cout << sys->yn[i] << " ";
    }
    cout << "\n" << endl;
    for (i = 0; i < sys->nz; i++){
        cout << sys->zn[i] << " ";
    }
    cout << "\n" << endl;*/
    /*vector<pair<double, int> > vx(xi.begin(), xi.end());
    sort(vx.begin(), vx.end(), comp);
    for (int i = 0; i < vx.size(); i++){
    cout << vx[i].first << " " << vx[i].second << ", ";
    }
    cout << endl;
    cout << sys->nx << endl;

    vector<pair<double, int> > vy(yi.begin(), yi.end());
    sort(vy.begin(), vy.end(), comp);
    for (int i = 0; i < vy.size(); i++){
    cout << vy[i].first << " " << vy[i].second << ", ";
    }
    cout << endl;
    cout << sys->ny << endl;

    vector<pair<double, int> > vz(zi.begin(), zi.end());
    sort(vz.begin(), vz.end(), comp);
    for (int i = 0; i < vz.size(); i++){
    cout << vz[i].first << " " <<vz[i].second << ", ";
    }
    cout << endl;
    cout << sys->nz << endl;*/
    /***********************************************************************************************/
    sys->xlim1 = sys->xn[0];
    sys->xlim2 = sys->xn[sys->nx - 1];
    sys->ylim1 = sys->yn[0];
    sys->ylim2 = sys->yn[sys->ny - 1];
    sys->zlim1 = sys->zn[0];
    sys->zlim2 = sys->zn[sys->nz - 1];
    sys->N_cell_x = sys->nx - 1;
    sys->N_cell_y = sys->ny - 1;
    sys->N_cell_z = sys->nz - 1;

    sys->N_edge_s = sys->N_cell_y*(sys->N_cell_x + 1) + sys->N_cell_x*(sys->N_cell_y + 1);
    sys->N_edge_v = (sys->N_cell_x + 1)*(sys->N_cell_y + 1);
    sys->N_edge = sys->N_edge_s*(sys->N_cell_z + 1) + sys->N_edge_v*(sys->N_cell_z);

    sys->N_node_s = sys->N_edge_v;
    sys->N_node = (sys->N_cell_z + 1)*sys->N_node_s;

    sys->N_patch_s = sys->N_cell_x*sys->N_cell_y;
    sys->N_patch_v = (sys->N_cell_x + 1)*sys->N_cell_y + (sys->N_cell_y + 1)*sys->N_cell_x;
    sys->N_patch = sys->N_patch_s*(sys->N_cell_z + 1) + sys->N_patch_v*sys->N_cell_z;

    sys->markEdge = (int*)calloc(sys->N_edge, sizeof(int));   // mark which edge is inside the conductor
    sys->markNode = (int*)calloc(sys->N_node, sizeof(int));   // mark which node is inside the conductor
    double xc, yc;

    for (i = 0; i < sys->numCdtRow; i++){
        //cout << polyIn((0 + 4.9e-7) / 2, -1.7e-7, sys, 2) << endl;
        numNode = (xi[sys->conductorIn[i].xmax] - xi[sys->conductorIn[i].xmin] + 1)
            *(yi[sys->conductorIn[i].ymax] - yi[sys->conductorIn[i].ymin] + 1)
            *(zi[sys->conductorIn[i].zmax] - zi[sys->conductorIn[i].zmin] + 1);
        //cout << sys->conductorIn[i].xmax << " " << sys->conductorIn[i].xmin << " " << sys->conductorIn[i].ymax << " " << sys->conductorIn[i].ymin << " " << sys->conductorIn[i].zmax << " " << sys->conductorIn[i].zmin << endl;
        sys->conductorIn[i].cdtInNode = (int*)malloc(numNode*sizeof(int));
        sys->conductorIn[i].numNode = 0;

        for (j = xi[sys->conductorIn[i].xmin]; j <= xi[sys->conductorIn[i].xmax]; j++){
            for (k = yi[sys->conductorIn[i].ymin]; k <= yi[sys->conductorIn[i].ymax]; k++){
                if (polyIn(sys->xn[j], sys->yn[k], sys, i)){
                    for (m = zi[sys->conductorIn[i].zmin]; m < zi[sys->conductorIn[i].zmax]; m++){
                        //cout << sys->xn[j] << " " << sys->yn[k] << endl;
                        sys->conductorIn[i].cdtInNode[sys->conductorIn[i].numNode] = m*sys->N_node_s + (sys->N_cell_y + 1)*j + k;
                        sys->conductorIn[i].numNode++;
                        sys->markNode[m*sys->N_node_s + (sys->N_cell_y + 1)*j + k] = 1;
                        if (sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] == 0){   // set the z direction markEdge
                            sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] = 1;
                            //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k << endl;
                        }
                    }
                    sys->conductorIn[i].cdtInNode[sys->conductorIn[i].numNode] = m*sys->N_node_s + (sys->N_cell_y + 1)*j + k;
                    sys->conductorIn[i].numNode++;
                    sys->markNode[m*sys->N_node_s + (sys->N_cell_y + 1)*j + k] = 1;
                }

            }
        }
        for (j = xi[sys->conductorIn[i].xmin]; j <= xi[sys->conductorIn[i].xmax]; j++){   // set the y direction markEdge
            for (k = yi[sys->conductorIn[i].ymin]; k < yi[sys->conductorIn[i].ymax]; k++){
                for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
                    xc = sys->xn[j];
                    yc = (sys->yn[k] + sys->yn[k + 1]) / 2;
                    if (polyIn(xc, yc, sys, i)){
                        sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k] = 1;
                        //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
                        //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
                    }
                }
            }
        }
        for (j = yi[sys->conductorIn[i].ymin]; j <= yi[sys->conductorIn[i].ymax]; j++){    // set the x direction markEdge
            for (k = xi[sys->conductorIn[i].xmin]; k < xi[sys->conductorIn[i].xmax]; k++){
                for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
                    xc = (sys->xn[k] + sys->xn[k + 1]) / 2;
                    yc = sys->yn[j];
                    //cout << xc << " " << yc << "" << i << " " << polyIn(xc, yc, sys, i) << endl;
                    if (polyIn(xc, yc, sys, i)){
                        sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j] = 1;
                        //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
                        //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
                    }
                }
            }
        }
    }

    /* construct edgelink */
    int eno;
    sys->edgelink = (int*)malloc(2 * sizeof(int)*sys->N_edge);
    for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++){
        for (i = 1; i <= sys->N_cell_x + 1; i++){    //edge along y axis
            for (j = 1; j <= sys->N_cell_y; j++){
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (i - 1)*sys->N_cell_y + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;

            }
        }
        for (i = 1; i <= sys->N_cell_x; i++){    //edge along x axis
            for (j = 1; j <= sys->N_cell_y + 1; j++){
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1)*sys->N_cell_y + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + i*(sys->N_cell_y + 1) + j - 1;

            }
        }
    }
    for (lyr = 1; lyr <= sys->N_cell_z; lyr++){    // edge along z axis
        for (i = 1; i <= sys->N_cell_x + 1; i++){
            for (j = 1; j <= sys->N_cell_y + 1; j++){
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = lyr*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
            }
        }
    }

    /* construct nodepos */
    int nno;
    sys->nodepos = (double*)malloc(sizeof(double)*sys->N_node * 3);   //N_node rows and 3 columns, input row by row
    for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++){
        for (i = 1; i <= sys->N_cell_x + 1; i++){
            for (j = 1; j <= sys->N_cell_y + 1; j++){
                nno = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->nodepos[(nno - 1) * 3 + 1 - 1] = sys->xn[i - 1];
                sys->nodepos[(nno - 1) * 3 + 2 - 1] = sys->yn[j - 1];
                sys->nodepos[(nno - 1) * 3 + 3 - 1] = sys->zn[lyr - 1];
            }
        }
    }

    /* construct nodeEdge */
    double leng;
    vector<pair<int, double> > a;
    for (i = 0; i < sys->N_edge; i++){
        sys->nodeEdge.push_back(a);
    }
    for (i = 0; i < sys->N_edge; i++){
        leng = pow((sys->nodepos[sys->edgelink[i * 2] * 3] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 1] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 1]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 2] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 2]), 2);
        leng = sqrt(leng);
        sys->nodeEdge[sys->edgelink[i * 2]].push_back(make_pair(i, 1 / leng));
        sys->nodeEdge[sys->edgelink[i * 2 + 1]].push_back(make_pair(i, -1 / leng));
    }
    
    /* implement dfs */
    //cout <<"Number of nodes: " << sys->N_node << endl;
    int *visited;
    vector<int> st;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    count = 0;
    for (i = 0; i < sys->N_node; i++){
        if (sys->markNode[i] == 0){
            continue;
        }
        else{
            if (visited[i] != 0){
                continue;
            }
            else{
                st.clear();
                st.push_back(i);
                count++;
                visited[i] = count;
                while (!st.empty()){ 
                    mark = 0;
                    for (j = 0; j < sys->nodeEdge[st.back()].size(); j++){
                        if (sys->markEdge[sys->nodeEdge[st.back()][j].first] != 0){
                            if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] == 0)){
                                visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] = count;
                                st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]);
                                mark = 1;
                                break;
                            }
                            else if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] == 0)){
                                visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] = count;
                                st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]);
                                mark = 1;
                                break;
                            }
                        }
                    }
                    if (mark == 0){
                        st.pop_back();
                    }
                }
            }
        }
    }
    for (i = 0; i < sys->N_edge; i++){
        if (sys->markEdge[i] != 0 && visited[sys->edgelink[i * 2]] == visited[sys->edgelink[i * 2 + 1]] && visited[sys->edgelink[i * 2]] != 0){
            sys->markEdge[i] = visited[sys->edgelink[i * 2 + 1]];    // Mark the edge with each color for different conductors
        }
    }
    for (i = 0; i < sys->N_node; i++){
        sys->markNode[i] = visited[i];
    }
    
    /* Construct each conductor */
    sys->numCdt = count;
    sys->conductor = (fdtdCdt*)malloc(sys->numCdt * sizeof(fdtdCdt));
    sys->cdtNumNode = (int*)calloc(sys->numCdt, sizeof(int));
    for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            sys->cdtNumNode[visited[i] - 1]++;
        }
    }
    for (i = 0; i < sys->numCdt; i++){
        sys->conductor[i].node = (int*)malloc(sizeof(int) * sys->cdtNumNode[i]);
        sys->conductor[i].cdtNodeind = 0;
    }
    for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            sys->conductor[visited[i] - 1].node[sys->conductor[visited[i] - 1].cdtNodeind] = i;
            sys->conductor[visited[i] - 1].cdtNodeind++;
        }
    }

    free(visited);
    visited = NULL;


    /* set markCell */
    vector<int> aa;
    vector<double> bb;
    sys->markCell = (int*)calloc(sys->N_cell_x * sys->N_cell_y * sys->N_cell_z, sizeof(int));
    for (i = 0; i < sys->N_edge; i++)
    {
        sys->edgeCell.push_back(aa);
        sys->edgeCellArea.push_back(bb);
    }
    int cell;
    for (i = 0; i < sys->N_cell_z; i++){
        for (j = 0; j < sys->N_cell_x; j++){
            for (k = 0; k < sys->N_cell_y; k++){
                mark = 1;
                cell = i * sys->N_patch_s + j * sys->N_cell_y + k;
                count = sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k];

                // y axis
                sys->edgeCell[(i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k)].push_back(cell);
                sys->edgeCellArea[(i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k)].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] == 0){
                    mark = 0;
                }

                sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back(cell);
                sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] != count){
                    mark = 0;
                }

                sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k].push_back(cell);
                sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
                if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] != count){
                    mark = 0;
                }

                sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back(cell);
                sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
                if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] != count){
                    mark = 0;
                }

                // x axis
                sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back(cell);
                sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back((sys->yn[k + 1] - sys->yn[k]) * (sys->zn[i + 1] - sys->zn[i]));
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] != count){
                    mark = 0;
                }

                sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back(cell);
                sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back((sys->yn[k + 1] - sys->yn[k])*(sys->zn[i + 1] - sys->zn[i])); 
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] != count){
                    mark = 0;
                }

                sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back(cell);
                sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back((sys->yn[k + 1] - sys->yn[k])*(sys->zn[i + 1] - sys->zn[i])); 
                if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] != count){
                    mark = 0;
                }

                sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back(cell);
                sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back((sys->yn[k + 1] - sys->yn[k]) * (sys->zn[i + 1] - sys->zn[i]));
                if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] != count){
                    mark = 0;
                }

                // z axis
                sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k].push_back(cell);
                sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k] != count){
                    mark = 0;
                }

                sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1].push_back(cell);
                sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1] != count){
                    mark = 0;
                }

                sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k].push_back(cell);
                sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k] != count){
                    mark = 0;
                }

                sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1].push_back(cell);
                sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
                if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1] != count){
                    mark = 0;
                }

                if (mark == 1){
                    sys->markCell[cell] = 1;
                }
            }
        }
    }


    return 0;
}

int matrixConstruction(fdtdMesh *sys){
    int i, j;
    ofstream outfile;
    
    /* construct D_eps */
    sys->eps = (double*)malloc(sizeof(double)*sys->N_edge);
    for (i = 0; i < sys->N_edge; i++){
        if (i < sys->N_edge_s){
            sys->eps[i] = sys->stackEpsn[0] * EPSILON0;
            continue;
        }
        else{
            sys->eps[i] = sys->stackEpsn[(i - sys->N_edge_s) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
        }
    }
    outfile.open("eps.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_edge; i++){
        outfile << sys->eps[i] << endl;
    }
    outfile.close();

    /* construct D_sig */
    sys->sig = (double*)calloc(sys->N_edge, sizeof(double));
    double a, b;
    for (i = 0; i < sys->N_edge; i++){
        if (sys->markEdge[i] != 0){
            a = 0;
            b = 0;
            for (j = 0; j < sys->edgeCell[i].size(); j++){
                a += sys->markCell[sys->edgeCell[i][j]] * sys->edgeCellArea[i][j];
                b += sys->edgeCellArea[i][j];
            }

            sys->sig[i] = (a / b) * SIGMA;
        }
    }

    outfile.open("sig.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_edge; i++){
        outfile << sys->sig[i] << endl;
    }
    outfile.close();

    return 0;
}

int portSet(fdtdMesh* sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi){

    char s[FDTD_MAXC];
    char *word[FDTD_MAXC];
    int j, i, l, m;
    vector<int> edge;

    for (i = 0; i < sys->numPorts; i++)
    {
        sys->portCoor[i].node = (int*)malloc((zi[sys->portCoor[i].z2] - zi[sys->portCoor[i].z1] + 1)*(xi[sys->portCoor[i].x2] - xi[sys->portCoor[i].x1] + 1)*(yi[sys->portCoor[i].y2] - yi[sys->portCoor[i].y1] + 1)*sizeof(int));
        sys->portCoor[i].nodenum = (zi[sys->portCoor[i].z2] - zi[sys->portCoor[i].z1] + 1)*(xi[sys->portCoor[i].x2] - xi[sys->portCoor[i].x1] + 1)*(yi[sys->portCoor[i].y2] - yi[sys->portCoor[i].y1] + 1);
        int n = 0;
        for (j = zi[sys->portCoor[i].z1]; j <= zi[sys->portCoor[i].z2]; j++){
            for (l = xi[sys->portCoor[i].x1]; l <= xi[sys->portCoor[i].x2]; l++){
                for (m = yi[sys->portCoor[i].y1]; m <= yi[sys->portCoor[i].y2]; m++){

                    sys->portCoor[i].node[n] = j * sys->N_node_s + l * (sys->N_cell_y + 1) + m;
                    sys->portNno.insert(j * sys->N_node_s + l * (sys->N_cell_y + 1) + m);

                    n++;

                }
            }
        }
        sys->portCoor[i].portCnd = sys->markNode[sys->portCoor[i].node[0]];
    }
    

    return 0;
}

bool polyIn(double x, double y, fdtdMesh *sys, int inPoly){
    int npol;
    int i, j, k;
    int isCond = 0;
    double disMin = 1e-10;

    npol = sys->conductorIn[inPoly].numVert;

    for (i = 0, j = npol - 1; i < npol; j = i++){
        if ((abs(y - sys->conductorIn[inPoly].y[j]) < disMin && abs(y - sys->conductorIn[inPoly].y[i]) < disMin &&
            ((x >= sys->conductorIn[inPoly].x[j] && x <= sys->conductorIn[inPoly].x[i]) ||
            (x >= sys->conductorIn[inPoly].x[i] && x <= sys->conductorIn[inPoly].x[j])))){
            return true;
        }
        else if (abs(x - sys->conductorIn[inPoly].x[j]) < disMin && abs(x - sys->conductorIn[inPoly].x[i]) < disMin &&
            ((y >= sys->conductorIn[inPoly].y[j] && y <= sys->conductorIn[inPoly].y[i]) ||
            (y >= sys->conductorIn[inPoly].y[i] && y <= sys->conductorIn[inPoly].y[j]))){
            return true;
        }
        else if ((abs(sys->conductorIn[inPoly].y[i] - sys->conductorIn[inPoly].y[j]) > disMin &&
            (((sys->conductorIn[inPoly].y[i] <= y) && (y < sys->conductorIn[inPoly].y[j])) ||
            ((sys->conductorIn[inPoly].y[j] <= y) && (y < sys->conductorIn[inPoly].y[i])))) &&
            (x < (sys->conductorIn[inPoly].x[j] - sys->conductorIn[inPoly].x[i]) * (y - sys->conductorIn[inPoly].y[i]) /
            (sys->conductorIn[inPoly].y[j] - sys->conductorIn[inPoly].y[i]) + sys->conductorIn[inPoly].x[i])){
            isCond = !isCond;
        }
    }
    return isCond;
}

int fdtdStringWord(char *s, char *word[]){
    int wno;                      /* word char counter */
    int ctr;                       /* char flag         */

    ctr = 1;
    wno = 0;
    for (; *s != '\0'; s++)
        if (*s == ' ' || *s == '\t' || *s == ',' || *s == '(' ||
            *s == ')' || *s == '\n') {
        *s = '\0';
        ctr = 1;
        }
        else {
            *s = (char)toupper(*s);
            if (ctr == 1) {
                word[wno++] = s;
                ctr = 0;
            }
        }
        return(wno);
}

double fdtdGetValue(char *str)
{
    double       value;                            /* converted true value      */

    sscanf(str, "%lf", &value);

    return value;

}

// compareString function necessary to treat C strings as equal if terminated OR newline encountered
int compareString(char *a, char *b)
{
    int i = 0;
    while ((a[i] != '\n' && a[i] != '\0') && (b[i] != '\n' && b[i] != '\0')){
        if (a[i] != b[i]){
            return 0;
        }
        i++;
    }
    if ((!(a[i] >= 'a' && a[i] <= 'z') && !(a[i] >= 'A' && a[i] <= 'Z') && !(a[i] >= '0' && a[i] <= '9')) && (!(b[i] >= 'a' && b[i] <= 'z') && !(b[i] >= 'A' && b[i] <= 'Z') && !(b[i] >= '0' && b[i] <= '9'))){
        return 1;
    }
    else{
        return 0;
    }
}

// Print fdtdMesh information
void fdtdMesh::print()
{
    // Print
    cout << "Contents of fdtdMesh" << endl;
    cout << " Length unit: " << this->lengthUnit << " m" << endl;
    cout << " Frequency sweep parameters: " << endl;
    cout << "  Frequency unit: " << this->freqUnit << " Hz" << endl;
    cout << "  Starting frequency: " << this->freqStart << " * " << this->freqUnit << " Hz" << endl;
    cout << "  Ending frequency: " << this->freqEnd << " * " << this->freqUnit << " Hz" << endl;
    cout << "  Number of frequencies: " << this->nfreq << endl;
    cout << "  Frequency scaling: " << this->freqScale << endl;
    cout << " Number of periods: " << this->ix << " in x-dir, and " << this->iy << " in y-dir" << endl;
    cout << " Node coordinate information:" << endl;
    cout << "  Number of nodes along direction: " << this->nx << " in x-dir, " << this->ny << " in y-dir, and " << this->nz << " in z-dir" << endl;
    cout << "  Node arrays: xn exists (" << (this->xn != nullptr) << "), yn exists (" << (this->yn != nullptr) << "), and zn exists (" << (this->zn != nullptr) << ")" << endl;
    cout << "  Node arrays: xnu exists (" << (this->xnu != nullptr) << "), ynu exists (" << (this->ynu != nullptr) << "), and znu exists (" << (this->znu != nullptr) << ")" << endl;
    cout << " Mesh cell information:" << endl;
    cout << "  Number in each direction: " << this->N_cell_x << " in x-dir, " << this->N_cell_y << " in y-dir, and " << this->N_cell_z << " in z-dir" << endl;
    cout << "  Factor in each direction: " << this->factor_x << " in x-dir, " << this->factor_y << " in y-dir, and " << this->factor_z << " in z-dir" << endl;
    cout << "  Limits in x-direction: " << this->xlim1 << " * " << this->lengthUnit << " m to " << this->xlim2 << " * " << this->lengthUnit << " m" << endl;
    cout << "  Limits in y-direction: " << this->ylim1 << " * " << this->lengthUnit << " m to " << this->ylim2 << " * " << this->lengthUnit << " m" << endl;
    cout << "  Limits in z-direction: " << this->zlim1 << " * " << this->lengthUnit << " m to " << this->zlim2 << " * " << this->lengthUnit << " m" << endl;
    cout << " Mesh edge information:" << endl;
    cout << "  Number of edges: " << this->N_edge << endl;
    cout << "  Number of edges_s: " << this->N_edge_s << endl;
    cout << "  Number of edges_v: " << this->N_edge_v << endl;
    cout << " Mesh node information:" << endl;
    cout << "  Number of nodes: " << this->N_node << endl;
    cout << "  Number of nodes_s: " << this->N_node_s << endl;
    cout << " Mesh patch information:" << endl;
    cout << "  Number of patches: " << this->N_patch << endl;
    cout << "  Number of patches_s: " << this->N_patch_s << endl;
    cout << "  Number of patches_v: " << this->N_patch_v << endl;
    cout << " Mesh field location information:" << endl;
    cout << "  nodepos array exists (" << (this->nodepos != nullptr) << ")" << endl;
    cout << "  Epoints array exists (" << (this->Epoints != nullptr) << ")" << endl;
    cout << "  edgelink array exists (" << (this->edgelink != nullptr) << ")" << endl;
    cout << "  Hpoints array exists (" << (this->Hpoints != nullptr) << ")" << endl;
    cout << "  nodeEdge vector has size " << this->nodeEdge.size() << endl;
    cout << " PEC information:" << endl;
    cout << "  Boundary node 1 array exists (" << (this->bd_node1 != nullptr) << ")" << endl;
    cout << "  Boundary node 2 array exists (" << (this->bd_node2 != nullptr) << ")" << endl;
    cout << "  Boundary edge array exists (" << (this->bd_node1 != nullptr) << ")" << endl;
    cout << " Layer stack material parameters:" << endl;
    cout << "  Number of layers in stack: " << this->numStack << endl;
    cout << "  Permittivity array exists (" << (this->stackEps != nullptr) << ")" << endl;
    cout << "  Beginning coordinate array exists (" << (this->stackBegCoor != nullptr) << ")" << endl;
    cout << "  Ending coordinate array exists (" << (this->stackEndCoor != nullptr) << ")" << endl;
    cout << "  Layer name vector has size " << this->stackName.size() << endl;
    cout << "  Other permittivity array eps exists (" << (this->eps != nullptr) << ")" << endl;
    cout << "  stackEpsn array exists (" << (this->stackEpsn != nullptr) << ")" << endl;
    cout << "  Stack conductor marker array exists (" << (this->stackCdtMark != nullptr) << ")" << endl;
    cout << " Conductor parameter information:" << endl;
    cout << "  conductorIn vector has size " << this->conductorIn.size() << endl;
    cout << "  Number of conductor rows: " << this->numCdtRow << endl;
    cout << "  Current number of conductor rows: " << this->numCdt << endl;
    cout << "  Edge marker array exists (" << (this->markEdge != nullptr) << ")" << endl;
    cout << "  Cell marker array exists (" << (this->markCell != nullptr) << ")" << endl;
    cout << "  cdtNumNode array exists (" << (this->cdtNumNode != nullptr) << ")" << endl;
    cout << "  sig array exist (" << (this->sig != nullptr) << ")" << endl;
    cout << "  conductor array exists (" << (this->conductor != nullptr) << ")" << endl;
    cout << "  Node marker array exists (" << (this->markNode != nullptr) << ")" << endl;
    cout << "  edgeCell vector has size " << this->edgeCell.size() << endl;
    cout << "  edgeCellArea vector has size " << this->edgeCellArea.size() << endl;
    cout << " Patch information exists (" << (this->patch != nullptr) << ")" << endl;
    cout << " Boundary information exists (" << (this->bound != nullptr) << ")" << endl;
    cout << " V0c information:" << endl;
    cout << "  v0cRowId vector has size " << this->v0cRowId.size() << endl;
    cout << "  v0cColId vector has size " << this->v0cColId.size() << endl;
    cout << "  v0cColIdo vector has size " << this->v0cColIdo.size() << endl;
    cout << "  v0cval vector has size " << this->v0cval.size() << endl;
    cout << "  v0cvalo vector has size " << this->v0cvalo.size() << endl;
    cout << "  v0c2RowId vector has size " << this->v0c2RowId.size() << endl;
    cout << "  v0c2ColId vector has size " << this->v0c2ColId.size() << endl;
    cout << "  v0c2ColIdo vector has size " << this->v0c2ColIdo.size() << endl;
    cout << "  v0c2val vector has size " << this->v0c2val.size() << endl;
    cout << "  v0caRowId vector has size " << this->v0caRowId.size() << endl;
    cout << "  v0caColId vector has size " << this->v0caColId.size() << endl;
    cout << "  v0caColIdo vector has size " << this->v0caColIdo.size() << endl;
    cout << "  v0caval vector has size " << this->v0caval.size() << endl;
    cout << "  v0cavalo vector has size " << this->v0cavalo.size() << endl;
    cout << " V0d_a information:" << endl;
    cout << "  v0d1aRowId vector has size " << this->v0d1aRowId.size() << endl;
    cout << "  v0d1aColId vector has size " << this->v0d1aColId.size() << endl;
    cout << "  v0d1aval vector has size " << this->v0d1aval.size() << endl;
    cout << "  v0d1avalo vector has size " << this->v0d1avalo.size() << endl;
    cout << "  v0d2aRowId vector has size " << this->v0d2aRowId.size() << endl;
    cout << "  v0d2aColId vector has size " << this->v0d2aColId.size() << endl;
    cout << "  v0d2aColIdo vector has size " << this->v0d2aColIdo.size() << endl;
    cout << "  v0d2aval vector has size " << this->v0d2aval.size() << endl;
    cout << " y0c2 information:" << endl;
    cout << "  y0c2RowId vector has size " << this->y0c2RowId.size() << endl;
    cout << "  y0c2ColId vector has size " << this->y0c2ColId.size() << endl;
    cout << "  y0c2val vector has size " << this->y0c2val.size() << endl;
    cout << " Arrays from technical paper information:" << endl;
    cout << "  v0c2y0c2 array exists (" << (this->v0c2y0c2 != nullptr) << ")" << endl;
    cout << "  v0c2y0c2o array exists (" << (this->v0c2y0c2o != nullptr) << ")" << endl;
    cout << "  yc array exists (" << (this->yc != nullptr) << ")" << endl;
    cout << "  v0cy0c array exists (" << (this->v0cy0c != nullptr) << ")" << endl;
    cout << " V0c'*D_sig*V0c infromation:" << endl;
    cout << "  AcRowId vector has size " << this->AcRowId.size() << endl;
    cout << "  AcColId vector has size " << this->AcColId.size() << endl;
    cout << "  Acval vector has size " << this->Acval.size() << endl;
    cout << "  AdRowId vector has size " << this->AdRowId.size() << endl;
    cout << "  AdColId vector has size " << this->AdColId.size() << endl;
    cout << "  Adval vector has size " << this->Adval.size() << endl;
    cout << "  crhsRowId vector has size " << this->crhsRowId.size() << endl;
    cout << "  crhsColId vector has size " << this->crhsColId.size() << endl;
    cout << "  crhsval vector has size " << this->crhsval.size() << endl;
    cout << "  crhs array exists (" << (this->crhs != nullptr) << ")" << endl;
    cout << " V0d_ information:" << endl;
    cout << "  v0d1RowId vector has size " << this->v0d1RowId.size() << endl;
    cout << "  v0d1ColId vector has size " << this->v0d1ColId.size() << endl;
    cout << "  v0d1val vector has size " << this->v0d1val.size() << endl;
    cout << "  v0d1valo vector has size " << this->v0d1valo.size() << endl;
    cout << "  v0d2RowId vector has size " << this->v0d2RowId.size() << endl;
    cout << "  v0d2ColId vector has size " << this->v0d2ColId.size() << endl;
    cout << "  v0d2ColIdo vector has size " << this->v0d2ColIdo.size() << endl;
    cout << "  v0d2val vector has size " << this->v0d2val.size() << endl;
    cout << "  v0d2valo vector has size " << this->v0d2valo.size() << endl;
    cout << "  yd array exists (" << (this->yd != nullptr) << ")" << endl;
    cout << " Solution storage information:" << endl;
    cout << "  y array exists (" << (this->y != nullptr) << ")" << endl;
    cout << "  x array exists (" << (this->x != nullptr) << ")" << endl;
    cout << " Port information:" << endl;
    cout << "  Number of ports: " << this->numPorts << endl;
    cout << "  portCoor array exists (" << (this->portCoor != nullptr) << ")" << endl;
    cout << "  portEdge vector has size " << this->portEdge.size() << endl;
    cout << "  portArea vector has size " << this->portArea.size() << endl;
    cout << "  portNno unordered set has values stored (" << (!(this->portNno.empty())) << ")" << endl;
    cout << " Current source array exists (" << (this->J != nullptr) << ")" << endl;
    cout << " Current V0c,s^T*I information:" << endl;
    cout << "  v0csJ array exists (" << (this->v0csJ != nullptr) << ")" << endl;
    cout << "  Y array exists (" << (this->Y != nullptr) << ")" << endl;
    cout << "------" << endl;
}
