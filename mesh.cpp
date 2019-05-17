<<<<<<< HEAD
//#include "stdafx.h"
#include "fdtd.h"



int meshAndMark(fdtdMesh *sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi, unordered_set<double> *portCoorx, unordered_set<double> *portCoory)
{
    int lyr;
    myint i, j, k, m;
    double upper, lower;
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
    double disMin = MINDIS;
    double disMaxx, disMaxy;   // the max discretization in x, y, z directions
    for (i = 0; i < sys->numCdtRow; i++) {
        numNode += sys->conductorIn[i].numVert;
    }
    //cout << numNode << endl;
    xOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    yOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    zOrigOld = (double*)calloc(2 * sys->numStack + 2 * sys->numPorts, sizeof(double));

    j = 0;
    for (i = 0; i < sys->numCdtRow; i++) {
        for (k = 0; k < sys->conductorIn[i].numVert; k++) {
            xOrigOld[j] = sys->conductorIn[i].x[k];
            yOrigOld[j] = sys->conductorIn[i].y[k];
            j++;
        }
    }

    for (i = 0; i < sys->numPorts; i++) {
        xOrigOld[j] = sys->portCoor[i].x1;
        yOrigOld[j] = sys->portCoor[i].y1;
        j++;
        xOrigOld[j] = sys->portCoor[i].x2;
        yOrigOld[j] = sys->portCoor[i].y2;
        j++;
    }
    j = 0;
    for (i = 0; i < sys->numStack; i++) {
        zOrigOld[j] = sys->stackBegCoor[i];
        j++;
        zOrigOld[j] = sys->stackEndCoor[i];
        j++;
    }
    for (i = 0; i < sys->numPorts; i++) {
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
    int xMaxInd = 10;
    disMaxx =  (xmax - xmin) / xMaxInd;
    //xMaxInd = (xmax - xmin) / disMaxx;

    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(xOrigOld[i] - xOrigOld[i - 1]) > disMin){
            sys->nx++;
        }
    }
    double *xn = (double*)calloc(numNode + 6 * sys->numPorts + xMaxInd, sizeof(double));
    xn[0] = xOrigOld[0];
    double temp = xn[0];
    j = 0;
    sys->nx = 1;
    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(xOrigOld[i] - temp) > disMin && abs(xOrigOld[i] - temp) <= disMaxx){
            j++;
            xn[j] = xOrigOld[i];
            temp = xn[j];
            sys->nx++;
        }
        else if (abs(xOrigOld[i] - temp) > disMin && abs(xOrigOld[i] - temp) > disMaxx){
            while (abs(xOrigOld[i] - temp) > disMaxx){
                j++;
                xn[j] = xn[j - 1] + disMaxx;
                temp = xn[j];
                sys->nx++;
            }
            if (abs(xOrigOld[i] - temp) > disMin){
                sys->nx++;
                temp = xOrigOld[i];
            }
            j++;
            xn[j] = xOrigOld[i];
        }
        else {
            j++;
            xn[j] = xOrigOld[i];
        }
    }
    int countx = j;
    //sort(xn, xn + countx + 1);


    //sys->xnu = (double*)calloc(sys->nx, sizeof(double));

    j = 0;
    //sys->xnu[0] = xn[0];
    //xi[sys->xnu[0]] = j;
    double first, second;
    temp = xn[0];
    for (i = 1; i <= countx; i++){    // set the discretization length around port to be equal
        if (abs(xn[i] - temp) > disMin){
            j++;
            temp = xn[i];
            //sys->xnu[j] = xn[i];
            //xi[sys->xnu[j]] = j;
        }
        else {
            //xi[xn[i]] = j;
        }
    }
    sys->nx = j + 1;




    /***************************************************************************/
    sort(yOrigOld, yOrigOld + numNode + 2 * sys->numPorts);
    
    sys->ny = 1;
    ymin = yOrigOld[0];
    ymax = yOrigOld[numNode + 2 * sys->numPorts - 1];
    int yMaxInd = 10;    // the max discretization of y is total / 120
    disMaxy =  (ymax - ymin) / yMaxInd;
    //yMaxInd = (ymax - ymin) / disMaxy;

    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(yOrigOld[i] - yOrigOld[i - 1]) > disMin){
            sys->ny++;
        }
    }

    double *yn = (double*)calloc(numNode + 6 * sys->numPorts + yMaxInd, sizeof(double));
    yn[0] = yOrigOld[0];
    j = 0;
    sys->ny = 1;
    
    temp = yn[0];
    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(yOrigOld[i] - temp) > disMin && abs(yOrigOld[i] - temp) <= disMaxy){
            j++;
            yn[j] = yOrigOld[i];
            temp = yn[j];
            //cout << yn[j] << " ";
            sys->ny++;
        }
        else if (abs(yOrigOld[i] - temp) > disMin && abs(yOrigOld[i] - temp) > disMaxy){
            while (abs(yOrigOld[i] - temp) > disMaxy){
                j++;
                yn[j] = yn[j - 1] + disMaxy;
                temp = yn[j];
                //cout << yn[j] << " " << yn[j - 1] << "    ";
                sys->ny++;
            }
            if (abs(yOrigOld[i] - temp) > disMin){
                sys->ny++;
                temp = yOrigOld[i];
            }
            j++;
            yn[j] = yOrigOld[i];

            //cout << yn[j] << " ";
        }
        else {
            j++;
            yn[j] = yOrigOld[i];
            //cout << yn[j] << " ";
        }
    }
    cout << endl;

    int county = j;
    //sort(yn, yn + county + 1);

    //sys->ynu = (double*)calloc(sys->ny, sizeof(double));


    j = 0;
    //sys->ynu[0] = yn[0];
    //yi[sys->ynu[0]] = j;
    temp = yn[0];
    for (i = 1; i <= county; i++){    // set the discretization length around port to be equal
        
        if (abs(yn[i] - temp) > disMin){
            j++;
            temp = yn[i];
            //sys->ynu[j] = yn[i];
            //yi[sys->ynu[j]] = j;
            
        }
        else {
            //yi[yn[i]] = j;
        }
    }

    sys->ny = j + 1;




    /********************************************************************************/
    sort(zOrigOld, zOrigOld + 2 * sys->numStack + 2 * sys->numPorts);
    sys->nz = 1;
    double disMinz = 1e-9;
    for (i = 1; i < 2 * sys->numStack + 2 * sys->numPorts; i++){
        if (abs(zOrigOld[i] - zOrigOld[i - 1]) > disMinz){
            sys->nz++;
        }
    }

    double *zn = (double*)calloc(2 * sys->numStack + 6 * sys->numPorts, sizeof(double));
    for (i = 0; i < 2 * sys->numStack + 2 * sys->numPorts; i++) {
        zn[i] = zOrigOld[i];
    }
    int countz = 2 * sys->numStack + 2 * sys->numPorts - 1;



    /*************************************************************************************/


    sort(xn, xn + countx + 1);
    xi.clear();
    sys->xn = (double*)calloc(sys->nx, sizeof(double));
    j = 0;
    sys->xn[0] = xn[0];
    temp = sys->xn[0];
    //xi[sys->xn[0]] = j;
    for (i = 1; i <= countx; i++){    // set the discretization length around port to be equal
        if (abs(xn[i] - temp) > disMin){
            j++;
            sys->xn[j] = xn[i];
            temp = sys->xn[j];
            xi[sys->xn[j]] = j;
        }
        else {
            xi[xn[i]] = j;
        }
    }
    sys->nx = j + 1;
    free(xn); xn = NULL;

    sort(yn, yn + county + 1);
    yi.clear();
    sys->yn = (double*)calloc(sys->ny, sizeof(double));
    j = 0;
    sys->yn[0] = yn[0];
    temp = sys->yn[0];
    yi[sys->yn[0]] = j;
    for (i = 1; i <= county; i++){    // set the discretization length around port to be equal

        if (abs(yn[i] - temp) > disMin){
            j++;
            sys->yn[j] = yn[i];
            temp = sys->yn[j];
            yi[sys->yn[j]] = j;
        }
        else {
            yi[yn[i]] = j;
        }
    }
    sys->ny = j + 1;
    free(yn); yn = NULL;

    sort(zn, zn + countz + 1);
    zi.clear();
    sys->zn = (double*)calloc(sys->nz, sizeof(double));
    j = 0;
    sys->zn[0] = zn[0];
    zi[sys->zn[0]] = j;
    for (i = 1; i <= countz; i++){    // set the discretization length around port to be equal
        if (abs(zn[i] - zn[i - 1]) > disMinz){
            j++;
            sys->zn[j] = zn[i];
            zi[sys->zn[j]] = j;
        }
        else {
            zi[zn[i]] = j;
        }
    }
    sys->nz = j + 1;
    free(zn); zn = NULL;

    sys->stackEpsn = (double*)calloc(sys->nz - 1, sizeof(double));
    i = 0;


    if (sys->stackBegCoor[0] == 0){
        j = 0;
        while (i < sys->nz - 1) {
            if ((sys->zn[i] + sys->zn[i + 1]) / 2 >= sys->stackBegCoor[j] && (sys->zn[i] + sys->zn[i + 1]) / 2 <= sys->stackEndCoor[j]) {
                sys->stackEpsn[i] = sys->stackEps[j];
                i++;
            }
            else {
                j++;
            }
        }
    }
    else {
        j = sys->numStack - 1;
        while (i < sys->nz - 1) {
            if ((sys->zn[i] + sys->zn[i + 1]) / 2 >= sys->stackBegCoor[j] && (sys->zn[i] + sys->zn[i + 1]) / 2 <= sys->stackEndCoor[j]) {
                sys->stackEpsn[i] = sys->stackEps[j];
                i++;
            }
            else {
                j--;
            }
        }
    }
    /*for (i = 0; i < sys->nz - 1; i++){
        cout << sys->stackEpsn[i] << endl;
        }*/
    for (i = 0; i < sys->nx; i++){
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
    cout << "\n" << endl;

    /***********************************************************************************************/

    sys->N_cell_x = sys->nx - (myint)1;
    sys->N_cell_y = sys->ny - (myint)1;
    sys->N_cell_z = sys->nz - (myint)1;


    sys->N_edge_s = sys->N_cell_y*(sys->N_cell_x + 1) + sys->N_cell_x*(sys->N_cell_y + 1);
    sys->N_edge_v = (sys->N_cell_x + 1)*(sys->N_cell_y + 1);
    sys->N_edge = sys->N_edge_s*(sys->N_cell_z + 1) + sys->N_edge_v*(sys->N_cell_z);

    sys->N_node_s = sys->N_edge_v;
    sys->N_node = (sys->N_cell_z + 1)*sys->N_node_s;

    sys->N_patch_s = sys->N_cell_x*sys->N_cell_y;
    sys->N_patch_v = (sys->N_cell_x + 1)*sys->N_cell_y + (sys->N_cell_y + 1)*sys->N_cell_x;
    sys->N_patch = sys->N_patch_s*(sys->N_cell_z + 1) + sys->N_patch_v*sys->N_cell_z;

    sys->markEdge = (myint*)calloc(sys->N_edge, sizeof(myint));   // mark which edge is inside the conductor
    sys->markNode = (myint*)calloc(sys->N_node, sizeof(myint));   // mark which node is inside the conductor

    cout << "N_edge = " << sys->N_edge << endl;
    cout << "N_node = " << sys->N_node << endl;
    cout << "N_cell_x = " << sys->N_cell_x << endl;
    cout << "N_cell_y = " << sys->N_cell_y << endl;
    cout << "N_cell_z = " << sys->N_cell_z << endl;
    double xc, yc;
    
    unordered_map<myint, myint> xrange;
    vector<myint> xcoorv;
    set<myint> xcoor;
    unordered_map<myint, set<myint>> xcoory;    // store the start coordinate of the range, the end can be checked from xrange
    myint ss, ee;
    myint y1, y2;
    myint l;
    clock_t tt = clock();
    int mark1;
    // Fast algorithm to find nodes inside conductors
    for (i = 0; i < sys->numCdtRow; i++){
        mark1 = 0;
        //cout << "Number of CdtRow is " << i << endl;
        for (j = 0; j < sys->conductorIn[i].numVert; j++){
            xcoor.insert(xi[sys->conductorIn[i].x[j]]);
        }


        for (auto xcoori = xcoor.begin(); xcoori != xcoor.end(); ++xcoori){
            xcoorv.push_back(*xcoori);
        }
        

        for (j = 0; j < xcoorv.size() - 1; j++){
            mark1 = 1;    // the x coordinates are more than 1
            xrange[xcoorv[j]] = xcoorv[j + 1];
            
        }
        if (xcoorv.size() == 1){    // If it has only one value
            xrange[xcoorv[0]] = xcoorv[0];
        }
        

        for (j = 0; j < sys->conductorIn[i].numVert - 1; j++){
            if (sys->conductorIn[i].y[j] == sys->conductorIn[i].y[j + 1]){
                if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[j + 1]){
                    ss = xi[sys->conductorIn[i].x[j]];
                    ee = xi[sys->conductorIn[i].x[j + 1]];
                    if (ss == ee && mark1 == 1){
                        continue;
                    }
                    while (xrange[ss] <= ee){
                        xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                        if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
                else if (sys->conductorIn[i].x[j] > sys->conductorIn[i].x[j + 1]){
                    ss = xi[sys->conductorIn[i].x[j + 1]];
                    ee = xi[sys->conductorIn[i].x[j]];
                    if (ss == ee && mark1 == 1){
                        continue;
                    }
                    while (xrange[ss] <= ee){
                        xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                        if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
            }
        }
        if (sys->conductorIn[i].y[j] == sys->conductorIn[i].y[0]){
            if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[0]){
                ss = xi[sys->conductorIn[i].x[j]];
                ee = xi[sys->conductorIn[i].x[0]];
                if (ss == ee && mark1 == 1){
                    continue;
                }
                while (xrange[ss] <= ee){
                    xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                    if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                        break;
                    }
                    ss = xrange[ss];
                }
            }
            else if (sys->conductorIn[i].x[j] > sys->conductorIn[i].x[0]){
                ss = xi[sys->conductorIn[i].x[0]];
                ee = xi[sys->conductorIn[i].x[j]];
                if (ss == ee && mark1 == 1){
                    continue;
                }
                while (xrange[ss] <= ee){
                    xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                    if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                        break;
                    }
                    ss = xrange[ss];
                }
            }
        }
        /*for (auto xcooryi : xcoory){
            cout << xcooryi.first << " : ";
            for (auto xcooryii : xcooryi.second){
                cout << xcooryii << " ";
            }
            cout << endl;
        }*/
        for (auto xcooryi : xcoory){
            mark = 0;
            for (j = xcooryi.first; j < xrange[xcooryi.first]; j++){
                mark = 0;
                for (auto xrangey : xcooryi.second){
                    mark++;
                    if (mark % 2 == 1){
                        y1 = xrangey;
                        continue;
                    }
                    else if (mark % 2 == 0){
                        y2 = xrangey;
                        
                        for (l = zi[sys->conductorIn[i].zmin]; l <= zi[sys->conductorIn[i].zmax]; l++){
                            for (k = y1; k < y2; k++){
                                sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] = i + 1;
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] = i + 1;
                                if (l != zi[sys->conductorIn[i].zmax]){
                                    sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                                }
                            }
                            sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                            sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] = i + 1;
                            if (l != zi[sys->conductorIn[i].zmax]){
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                            }
                        }
                        
                    }
                }
            }
            mark = 0;
            for (auto xrangey : xcooryi.second){
                mark++;
                if (mark % 2 == 1){
                    y1 = xrangey;
                    continue;
                }
                else if (mark % 2 == 0){
                    y2 = xrangey;
                    for (l = zi[sys->conductorIn[i].zmin]; l <= zi[sys->conductorIn[i].zmax]; l++){
                        for (k = y1; k < y2; k++){
                            sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                            sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] = i + 1;
                            if (l != zi[sys->conductorIn[i].zmax]){
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                            }
                        }
                        sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                        if (l != zi[sys->conductorIn[i].zmax]){
                            sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                        }
                    }
                }
            }
        }
        xrange.clear();
        xcoor.clear();
        xcoorv.clear();
        xcoory.clear();
    }
    cout << "The time to mark Edge and Node is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << endl;
    //for (i = 0; i < sys->numCdtRow; i++){
    //    //cout << polyIn((0 + 4.9e-7) / 2, -1.7e-7, sys, 2) << endl;
    //    numNode = (xi[sys->conductorIn[i].xmax] - xi[sys->conductorIn[i].xmin] + 1)
    //        *(yi[sys->conductorIn[i].ymax] - yi[sys->conductorIn[i].ymin] + 1)
    //        *(zi[sys->conductorIn[i].zmax] - zi[sys->conductorIn[i].zmin] + 1);
    //    //cout << sys->conductorIn[i].xmax << " " << sys->conductorIn[i].xmin << " " << sys->conductorIn[i].ymax << " " << sys->conductorIn[i].ymin << " " << sys->conductorIn[i].zmax << " " << sys->conductorIn[i].zmin << endl;
    //    sys->conductorIn[i].cdtInNode = (myint*)malloc(numNode*sizeof(myint));
    //    sys->conductorIn[i].numNode = 0;
    //    

    //    for (j = xi[sys->conductorIn[i].xmin]; j <= xi[sys->conductorIn[i].xmax]; j++){
    //        for (k = yi[sys->conductorIn[i].ymin]; k <= yi[sys->conductorIn[i].ymax]; k++){
    //            if (polyIn(sys->xn[j], sys->yn[k], sys, i)){
    //                for (m = zi[sys->conductorIn[i].zmin]; m < zi[sys->conductorIn[i].zmax]; m++){
    //                    //cout << sys->xn[j] << " " << sys->yn[k] << endl;
    //                    sys->conductorIn[i].cdtInNode[sys->conductorIn[i].numNode] = m*sys->N_node_s + (sys->N_cell_y + 1)*j + k;
    //                    sys->conductorIn[i].numNode++;
    //                    sys->markNode[m*sys->N_node_s + (sys->N_cell_y + 1)*j + k] = 1;
    //                    if (sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] == 0){   // set the z direction markEdge
    //                        sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] = i + 1;
    //                        //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k << endl;
    //                    }
    //                }
    //                sys->conductorIn[i].cdtInNode[sys->conductorIn[i].numNode] = m*sys->N_node_s + (sys->N_cell_y + 1)*j + k;
    //                sys->conductorIn[i].numNode++;
    //                sys->markNode[m*sys->N_node_s + (sys->N_cell_y + 1)*j + k] = 1;
    //            }

    //        }
    //    }
    //    for (j = xi[sys->conductorIn[i].xmin]; j <= xi[sys->conductorIn[i].xmax]; j++){   // set the y direction markEdge
    //        for (k = yi[sys->conductorIn[i].ymin]; k < yi[sys->conductorIn[i].ymax]; k++){
    //            for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
    //                xc = sys->xn[j];
    //                yc = (sys->yn[k] + sys->yn[k + 1]) / 2;
    //                if (polyIn(xc, yc, sys, i)){
    //                    sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k] = i + 1;
    //                    //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
    //                    //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
    //                }
    //            }
    //        }
    //    }
    //    for (j = yi[sys->conductorIn[i].ymin]; j <= yi[sys->conductorIn[i].ymax]; j++){    // set the x direction markEdge
    //        for (k = xi[sys->conductorIn[i].xmin]; k < xi[sys->conductorIn[i].xmax]; k++){
    //            for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
    //                xc = (sys->xn[k] + sys->xn[k + 1]) / 2;
    //                yc = sys->yn[j];
    //                //cout << xc << " " << yc << "" << i << " " << polyIn(xc, yc, sys, i) << endl;
    //                if (polyIn(xc, yc, sys, i)){
    //                    sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j] = i + 1;
    //                    //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
    //                    //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
    //                }
    //            }
    //        }
    //    }
    //}

    for (i = 0; i < sys->N_edge_s; i++){    // the lower plane
        sys->markEdge[i] = sys->numCdtRow + 1;
    }
    for (i = 0; i < sys->N_node_s; i++){
        sys->markNode[i] = 1;
    }
    for (i = sys->N_edge - sys->N_edge_s; i < sys->N_edge; i++){    // the upper plane
        sys->markEdge[i] = sys->numCdtRow + 2;
    }
    for (i = sys->N_node - sys->N_node_s; i < sys->N_node; i++){
        sys->markNode[i] = 1;
    }

    /* construct edgelink */
    myint eno;
    sys->edgelink = (myint*)malloc(2 * sizeof(myint)*sys->N_edge);
    for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++) {
        for (i = 1; i <= sys->N_cell_x + 1; i++) {    //edge along y axis
            for (j = 1; j <= sys->N_cell_y; j++) {
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (i - 1)*sys->N_cell_y + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;

            }
        }
        for (i = 1; i <= sys->N_cell_x; i++) {    //edge along x axis
            for (j = 1; j <= sys->N_cell_y + 1; j++) {
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1)*sys->N_cell_y + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + i*(sys->N_cell_y + 1) + j - 1;

            }
        }
    }
    for (lyr = 1; lyr <= sys->N_cell_z; lyr++) {    // edge along z axis
        for (i = 1; i <= sys->N_cell_x + 1; i++) {
            for (j = 1; j <= sys->N_cell_y + 1; j++) {
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = lyr*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
            }
        }
    }

    /* construct nodepos */
    myint nno;

    sys->nodepos = (double*)malloc(sizeof(double)*sys->N_node * 3);   //N_node rows and 3 columns, input row by row
    for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++) {
        for (i = 1; i <= sys->N_cell_x + 1; i++) {
            for (j = 1; j <= sys->N_cell_y + 1; j++) {
                nno = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->nodepos[(nno - 1) * 3 + 1 - 1] = sys->xn[i - 1];
                sys->nodepos[(nno - 1) * 3 + 2 - 1] = sys->yn[j - 1];
                sys->nodepos[(nno - 1) * 3 + 3 - 1] = sys->zn[lyr - 1];
            }
        }
    }

    /* construct nodeEdge */
    double leng;
    vector<pair<myint, double> > a;
    for (i = 0; i < sys->N_node; i++) {
        sys->nodeEdge.push_back(a);
        sys->nodeEdgea.push_back(a);
    }
    for (i = 0; i < sys->N_edge; i++) {
        leng = pow((sys->nodepos[sys->edgelink[i * 2] * 3] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 1] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 1]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 2] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 2]), 2);
        leng = sqrt(leng);
        sys->nodeEdge[sys->edgelink[i * 2]].push_back(make_pair(i, 1 / leng));
        sys->nodeEdge[sys->edgelink[i * 2 + 1]].push_back(make_pair(i, -1 / leng));
    }
    cout << "nodeEdge is done\n";
    int ix, iy, iz;
    iz = sys->N_cell_z;
    for (ix = 0; ix < sys->nx; ix++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, -1 / (sys->zn[iz] - sys->zn[iz - 1])));
        }
    }
    iz = 0;
    for (ix = 0; ix < sys->nx; ix++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->zn[iz + 1] - sys->zn[iz])));
        }
    }
    for (iz = 1; iz < sys->N_cell_z; iz++) {
        for (ix = 0; ix < sys->nx; ix++) {
            for (iy = 0; iy < sys->ny; iy++) {
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, -2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
            }
        }
    }
    ix = sys->N_cell_x;
    for (iz = 0; iz < sys->nz; iz++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -1 / (sys->xn[ix] - sys->xn[ix - 1])));
        }
    }
    ix = 0;
    for (iz = 0; iz < sys->nz; iz++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->xn[ix + 1] - sys->xn[ix])));
        }
    }
    for (ix = 1; ix < sys->N_cell_x; ix++) {
        for (iz = 0; iz < sys->nz; iz++) {
            for (iy = 0; iy < sys->ny; iy++) {
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
            }
        }
    }
    iy = sys->N_cell_y;
    for (iz = 0; iz < sys->nz; iz++) {
        for (ix = 0; ix < sys->nx; ix++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -1 / (sys->yn[iy] - sys->yn[iy - 1])));
        }
    }
    iy = 0;
    for (iz = 0; iz < sys->nz; iz++) {
        for (ix = 0; ix < sys->nx; ix++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 1 / (sys->yn[iy + 1] - sys->yn[iy])));
        }
    }
    for (iy = 1; iy < sys->N_cell_y; iy++) {
        for (iz = 0; iz < sys->nz; iz++) {
            for (ix = 0; ix < sys->nx; ix++) {
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
            }
        }
    }

    /* implement dfs */
    //cout <<"Number of nodes: " << sys->N_node << endl;
    myint *visited;
    vector<int> st;
    unordered_set<int> base;
    visited = (myint*)calloc(sys->N_node, sizeof(myint));
    count = 0;

    for (i = 0; i < sys->N_node; i++) {
        if (sys->markNode[i] == 0) {
            continue;
        }
        else {
            if (visited[i] != 0) {
                continue;
            }
            else {
                st.clear();
                st.push_back(i);
                count++;
                visited[i] = count;
                sys->cond2condIn.push_back(base);
                while (!st.empty()) {
                    mark = 0;
                    for (j = 0; j < sys->nodeEdge[st.back()].size(); j++) {
                        if (sys->markEdge[sys->nodeEdge[st.back()][j].first] != 0) {
                            if (sys->cond2condIn[count - 1].find(sys->markEdge[sys->nodeEdge[st.back()][j].first]) == sys->cond2condIn[count - 1].end()) {
                                sys->cond2condIn[count - 1].insert(sys->markEdge[sys->nodeEdge[st.back()][j].first]);
                            }
                            if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] == 0)) {
                                visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] = count;
                                st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]);
                                mark = 1;

                                break;
                            }
                            else if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] == 0)) {
                                visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] = count;
                                st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]);
                                mark = 1;

                                break;
                            }

                        }
                    }
                    if (mark == 0) {
                        st.pop_back();
                    }
                }
            }
        }
    }
    
    cout << "Number of conductors is " << count << endl;
    /*for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            for (j = 0; j < sys->nodeEdge[i].size(); j++){
                if (sys->markEdge[sys->nodeEdge[i][j].first] != 0 && sys->cond2condIn[visited[i] - 1].find(sys->markEdge[sys->nodeEdge[i][j].first]) == sys->cond2condIn[visited[i] - 1].end()){
                    sys->cond2condIn[visited[i] - 1].insert(sys->markEdge[sys->nodeEdge[i][j].first]);
                }
            }
        }
    }*/


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
    /*for (i = 0; i < sys->numCdt; i++) {
        cout << "Cnt " << i << " has condIn as: ";
        for (auto ci : sys->cond2condIn[i]) {
            cout << ci << " ";
        }
        cout << endl;
    }*/
    sys->conductor = (fdtdCdt*)malloc(sys->numCdt * sizeof(fdtdCdt));
    sys->cdtNumNode = (myint*)calloc(sys->numCdt, sizeof(myint));
    for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            sys->cdtNumNode[visited[i] - 1]++;
        }
    }
    for (i = 0; i < sys->numCdt; i++){
        sys->conductor[i].node = (myint*)malloc(sizeof(myint) * sys->cdtNumNode[i]);
        sys->conductor[i].cdtNodeind = 0;
        sys->conductor[i].markPort = 0;
    }
    int portground_count = 0;
    for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            sys->conductor[visited[i] - 1].node[sys->conductor[visited[i] - 1].cdtNodeind] = i;
            if ((i < sys->N_node_s) && sys->conductor[visited[i] - 1].markPort != -1){
                /*portground_count++;
                if (portground_count <= 1)*/
                    sys->conductor[visited[i] - 1].markPort = -1;    // this conductor connect to the lower PEC
            }
            else if ((i >= sys->N_node - sys->N_node_s) && sys->conductor[visited[i] - 1].markPort != -2){
                sys->conductor[visited[i] - 1].markPort = -2;    // this conductor connects to the upper PEC
            }
            sys->conductor[visited[i] - 1].cdtNodeind++;
        }
    }
    free(visited);
    visited = NULL;
    
    /* set markCell */
    //vector<int> aa;
    //vector<double> bb;
    //sys->markCell = (int*)calloc(sys->N_cell_x * sys->N_cell_y * sys->N_cell_z, sizeof(int));
    //for (i = 0; i < sys->N_edge; i++)
    //{
    //    sys->edgeCell.push_back(aa);
    //    sys->edgeCellArea.push_back(bb);
    //}
    //int cell;
    //for (i = 0; i < sys->N_cell_z; i++){
    //    for (j = 0; j < sys->N_cell_x; j++){
    //        for (k = 0; k < sys->N_cell_y; k++){
    //            mark = 1;
    //            cell = i * sys->N_patch_s + j * sys->N_cell_y + k;
    //            count = sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k];

    //            // y axis
    //            sys->edgeCell[(i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k)].push_back(cell);
    //            sys->edgeCellArea[(i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k)].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] == 0){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] != count){
    //                mark = 0;
    //            }

    //            // x axis
    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back((sys->yn[k + 1] - sys->yn[k]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back((sys->yn[k + 1] - sys->yn[k])*(sys->zn[i + 1] - sys->zn[i])); 
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back((sys->yn[k + 1] - sys->yn[k])*(sys->zn[i + 1] - sys->zn[i])); 
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back((sys->yn[k + 1] - sys->yn[k]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            // z axis
    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            if (mark == 1){
    //                sys->markCell[cell] = 1;
    //            }
    //        }
    //    }
    //}
    


    return 0;
}

int matrixConstruction(fdtdMesh *sys){
    int i, j;

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
    ofstream out;
    /*out.open("eps.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = sys->N_edge_s; i < sys->N_edge - sys->N_edge_s; i++){
        out << sys->eps[i] << endl;
    }
    out.close();*/

    /* construct D_sig */
    sys->sig = (double*)calloc(sys->N_edge, sizeof(double));
    double a, b;
    for (i = 0; i < sys->N_edge; i++){
        if (sys->markEdge[i] != 0){
            a = 0;
            b = 0;
            /*for (j = 0; j < sys->edgeCell[i].size(); j++){
                a += sys->markCell[sys->edgeCell[i][j]] * sys->edgeCellArea[i][j];
                b += sys->edgeCellArea[i][j];
            }*/
            sys->sig[i] = SIGMA;/*(a / b) * SIGMA;*/
        }
    }
    /*out.open("sig.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = sys->N_edge_s; i < sys->N_edge - sys->N_edge_s; i++){
        out << sys->sig[i] << endl;
    }
    out.close();*/

    sys->edgeCell.clear();
    sys->edgeCellArea.clear();
    //free(sys->markCell); sys->markCell = NULL;

    return 0;
}

int portSet(fdtdMesh* sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi){
    char s[FDTD_MAXC];
    char *word[FDTD_MAXC];
    int j, i, mark, l, m, k;
    int node;
    vector<int> edge;
    double a;
    double sideLen = 0.;

    for (i = 0; i < sys->numPorts; i++)
    {
        
        
        if (sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] != 0 && sys->conductor[sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] - 1].markPort > -1){
            sys->portCoor[i].portCnd = sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]];
            sys->conductor[sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] - 1].markPort = i + 1;    // markPort start from 1
        }
        else if (sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] != 0 && sys->conductor[sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] - 1].markPort > -1){
            sys->portCoor[i].portCnd = sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]];
            sys->conductor[sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] - 1].markPort = i + 1;    // markPort start from 1
        }
        
        //cout << sys->portCoor[i].portCnd << endl;
        edge.clear();
        if (sys->portCoor[i].x1 != sys->portCoor[i].x2){
            if (sys->portCoor[i].x1 < sys->portCoor[i].x2){
                for (j = xi[sys->portCoor[i].x1]; j < xi[sys->portCoor[i].x2]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*(sys->N_cell_x + 1) + j*(sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]);
                }
            }
            else{
                for (j = xi[sys->portCoor[i].x2]; j < xi[sys->portCoor[i].x1]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*(sys->N_cell_x + 1) + j*(sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]);
                }
            }
            a = 1;
            if (yi[sys->portCoor[i].y1] == 0){
                a *= (sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1);
            }
            else if (yi[sys->portCoor[i].y1] == sys->N_cell_y){
                a *= (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]);
            }
            else{
                a *= ((sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1) / 2 + (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]) / 2);
            }
            if (zi[sys->portCoor[i].z1] == 0){
                a *=(sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
            }
            else if (zi[sys->portCoor[i].z1] == sys->N_cell_z){
                a *= (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]);
            }
            else{
                a *= ((sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1) / 2 + (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]) / 2);
            }
            sys->portArea.push_back(a);
        }
        else if (sys->portCoor[i].y1 != sys->portCoor[i].y2){
            if (sys->portCoor[i].y1 < sys->portCoor[i].y2){
                for (j = yi[sys->portCoor[i].y1]; j < yi[sys->portCoor[i].y2]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*xi[sys->portCoor[i].x1] + j);
                }
            }
            else{
                for (j = yi[sys->portCoor[i].y2]; j < yi[sys->portCoor[i].y1]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*xi[sys->portCoor[i].x1] + j);
                }
            }
            a = 1;
            if (xi[sys->portCoor[i].x1] == 0){
                a *= (sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1);
            }
            else if (xi[sys->portCoor[i].x1] == sys->N_cell_x){
                a *= (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]);
            }
            else{
                a *= ((sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1) / 2 + (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]) / 2);
            }
            if (zi[sys->portCoor[i].z1] == 0){
                a *=(sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
            }
            else if (zi[sys->portCoor[i].z1] == sys->N_cell_z){
                a *= (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]);
            }
            else{
                a *= ((sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1) / 2 + (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]) / 2);
            }
            sys->portArea.push_back(a);
        }
        else if (sys->portCoor[i].z1 != sys->portCoor[i].z2){
            if (sys->portCoor[i].z1 < sys->portCoor[i].z2){
                for (j = zi[sys->portCoor[i].z1]; j < zi[sys->portCoor[i].z2]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*j + sys->N_edge_s + (sys->N_cell_y + 1)*xi[sys->portCoor[i].x1] + yi[sys->portCoor[i].y1]);
                }
            }
            else{
                for (j = zi[sys->portCoor[i].z2]; j < zi[sys->portCoor[i].z1]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*j + sys->N_edge_s + (sys->N_cell_y + 1)*xi[sys->portCoor[i].x1] + yi[sys->portCoor[i].y1]);
                }
            }
            a = 1;
            if (xi[sys->portCoor[i].x1] == 0){
                a *= (sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1);
            }
            else if (xi[sys->portCoor[i].x1] == sys->N_cell_x){
                a *= (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]);
            }
            else{
                a *= ((sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1) / 2 + (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]) / 2);
            }
            if (yi[sys->portCoor[i].y1] == 0){
                a *=(sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1);
            }
            else if (yi[sys->portCoor[i].y1] == sys->N_cell_y){
                a *= (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]);
            }
            else{
                a *= ((sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1) / 2 + (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]) / 2);
            }
            sys->portArea.push_back(a);

        }
        sys->portEdge.push_back(edge);
        
    }

    clock_t t1 = clock();
    sys->markProSide = (int*)calloc(sys->N_node, sizeof(int));
    double x1, x2, y1, y2;
    myint x1_ind, x2_ind, y1_ind, y2_ind, z1_ind, z2_ind;

    for (i = 0; i < sys->numPorts; i++){
        /*cout << "Value of sys->portCoor[i].portCnd - 1: " << sys->portCoor[i].portCnd - 1 << endl;
        cout << "Size of the sys->cond2condIn[sys->portCoor[i].portCnd - 1] unordered_set: " << sys->cond2condIn[sys->portCoor[i].portCnd - 1].size() << endl; */
        for (auto ci : sys->cond2condIn[sys->portCoor[i].portCnd - 1]){
            for (l = 0; l < sys->conductorIn[ci - 1].numVert - 1; l++){
                if (sys->conductorIn[ci - 1].x[l] == sys->conductorIn[ci - 1].x[l + 1]){
                    x1 = sys->conductorIn[ci - 1].x[l];
                    x2 = x1;
                    if (sys->conductorIn[ci - 1].y[l] < sys->conductorIn[ci - 1].y[l + 1]){
                        y1 = sys->conductorIn[ci - 1].y[l];
                        y2 = sys->conductorIn[ci - 1].y[l + 1];
                    }
                    else{
                        y1 = sys->conductorIn[ci - 1].y[l + 1];
                        y2 = sys->conductorIn[ci - 1].y[l];
                    }
                    x1_ind = xi[x1];
                    x2_ind = xi[x2];
                    y1_ind = yi[y1];
                    y2_ind = yi[y2];
                    while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                        x1_ind--;
                    }
                    x1_ind++;
                    while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                        x2_ind++;
                    }
                    x2_ind--;
                    while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                        y1_ind--;
                    }
                    y1_ind++;
                    while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                        y2_ind++;
                    }
                    y2_ind--;
                    z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                    z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                    for (k = z1_ind; k <= z2_ind; k++){
                        for (j = x1_ind; j <= x2_ind; j++){
                            for (m = y1_ind; m <= y2_ind; m++){
                                if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                    sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                            }
                        }
                    }
                }
                else{
                    y1 = sys->conductorIn[ci - 1].y[l];
                    y2 = y1;
                    if (sys->conductorIn[ci - 1].x[l] < sys->conductorIn[ci - 1].x[l + 1]){
                        x1 = sys->conductorIn[ci - 1].x[l];
                        x2 = sys->conductorIn[ci - 1].x[l + 1];
                    }
                    else{
                        x1 = sys->conductorIn[ci - 1].x[l + 1];
                        x2 = sys->conductorIn[ci - 1].x[l];
                    }
                    x1_ind = xi[x1];
                    x2_ind = xi[x2];
                    y1_ind = yi[y1];
                    y2_ind = yi[y2];
                    while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                        x1_ind--;
                    }
                    x1_ind++;
                    while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                        x2_ind++;
                    }
                    x2_ind--;
                    while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                        y1_ind--;
                    }
                    y1_ind++;
                    while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                        y2_ind++;
                    }
                    y2_ind--;
                    z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                    z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                    for (k = z1_ind; k <= z2_ind; k++){
                        for (j = x1_ind; j <= x2_ind; j++){
                            for (m = y1_ind; m <= y2_ind; m++){
                                if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                    sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                            }
                        }
                    }
                }
                
            }
            if (sys->conductorIn[ci - 1].x[l] == sys->conductorIn[ci - 1].x[0]){
                x1 = sys->conductorIn[ci - 1].x[l];
                x2 = x1;
                if (sys->conductorIn[ci - 1].y[l] < sys->conductorIn[ci - 1].y[0]){
                    y1 = sys->conductorIn[ci - 1].y[l];
                    y2 = sys->conductorIn[ci - 1].y[0];
                }
                else{
                    y1 = sys->conductorIn[ci - 1].y[0];
                    y2 = sys->conductorIn[ci - 1].y[l];
                }
                x1_ind = xi[x1];
                x2_ind = xi[x2];
                y1_ind = yi[y1];
                y2_ind = yi[y2];
                while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                    x1_ind--;
                }
                x1_ind++;
                while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                    x2_ind++;
                }
                x2_ind--;
                while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                    y1_ind--;
                }
                y1_ind++;
                while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                    y2_ind++;
                }
                y2_ind--;
                z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                for (k = z1_ind; k <= z2_ind; k++){
                    for (j = x1_ind; j <= x2_ind; j++){
                        for (m = y1_ind; m <= y2_ind; m++){
                            if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                        }
                    }
                }
            }
            else{
                y1 = sys->conductorIn[ci - 1].y[l];
                y2 = y1;
                if (sys->conductorIn[ci - 1].x[l] < sys->conductorIn[ci - 1].x[0]){
                    x1 = sys->conductorIn[ci - 1].x[l];
                    x2 = sys->conductorIn[ci - 1].x[0];
                }
                else{
                    x1 = sys->conductorIn[ci - 1].x[0];
                    x2 = sys->conductorIn[ci - 1].x[l];
                }
                x1_ind = xi[x1];
                x2_ind = xi[x2];
                y1_ind = yi[y1];
                y2_ind = yi[y2];
                while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                    x1_ind--;
                }
                x1_ind++;
                while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                    x2_ind++;
                }
                x2_ind--;
                while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                    y1_ind--;
                }
                y1_ind++;
                while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                    y2_ind++;
                }
                y2_ind--;
                z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                for (k = z1_ind; k <= z2_ind; k++){
                    for (j = x1_ind; j <= x2_ind; j++){
                        for (m = y1_ind; m <= y2_ind; m++){
                            if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                        }
                    }
                }
            }
        }
    }
    cout << "Time of finding side nodes is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;

    sys->conductorIn.clear();
    return 0;
}

// Is point (x,y) within the polygon?
bool polyIn(double x, double y, fdtdMesh *sys, int inPoly){
    int npol;
    int i, j, k;
    bool isCond = false;
    double disMin = 1.e-10;

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

// Print fdtdPort information
void fdtdPort::print()
{
    // Print
    cout << " ------" << endl;
    cout << " Contents of fdtdPort" << endl;
    if (this->x == nullptr)
    {
        cout << "  x array exists (" << (this->x != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  x array has size " << NELEMENT(this->x) << endl;
    }
    if (this->y == nullptr)
    {
        cout << "  y array exists (" << (this->y != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  y array has size " << NELEMENT(this->y) << endl;
    }
    if (this->z == nullptr)
    {
        cout << "  z array exists (" << (this->z != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  z array has size " << NELEMENT(this->z) << endl;
    }
    cout << "  Coordinates of two points for the source:" << endl;
    cout << "   x-coordinates: " << this->x1 << " m and " << this->x2 << " m" << endl;
    cout << "   y-coordinates: " << this->y1 << " m and " << this->y2 << " m" << endl;
    cout << "   z-coordinates: " << this->z1 << " m and " << this->z2 << " m" << endl;
    cout << "  Port direction: " << this->portDirection << endl;
    if (this->node == nullptr)
    {
        cout << "  Node array exists (" << (this->node != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Node array has size " << NELEMENT(this->node) << endl;
    }
    cout << "  Number of nodes: " << this->nodenum << endl;
    cout << "  portCnd: " << this->portCnd << endl;
    cout << " ------" << endl;
}

// Print fdtdMesh information
void fdtdMesh::print()
{
    // Print
    cout << "------" << endl;
    cout << "Contents of fdtdMesh" << endl;
    cout << " Length unit: " << this->lengthUnit << " m" << endl;
    cout << " Frequency sweep parameters: " << endl;
    cout << "  Frequency unit: " << this->freqUnit << " Hz" << endl;
    cout << "  Starting frequency: " << this->freqStart << " * " << this->freqUnit << " Hz" << endl;
    cout << "  Ending frequency: " << this->freqEnd << " * " << this->freqUnit << " Hz" << endl;
    cout << "  Number of frequencies: " << this->nfreq << endl;
    cout << "  Frequency scaling: " << this->freqScale << endl;
    cout << " Node coordinate information:" << endl;
    cout << "  Number of nodes along direction: " << this->nx << " in x-dir, " << this->ny << " in y-dir, and " << this->nz << " in z-dir" << endl;
    cout << "  Node arrays: xn exists (" << (this->xn != nullptr) << "), yn exists (" << (this->yn != nullptr) << "), and zn exists (" << (this->zn != nullptr) << ")" << endl;
    cout << "  Node arrays: xnu exists (" << (this->xnu != nullptr) << "), ynu exists (" << (this->ynu != nullptr) << "), and znu exists (" << (this->znu != nullptr) << ")" << endl;
    cout << " Mesh cell information:" << endl;
    cout << "  Number in each direction: " << this->N_cell_x << " in x-dir, " << this->N_cell_y << " in y-dir, and " << this->N_cell_z << " in z-dir" << endl;
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
    cout << "  nodeEdgea vector has size " << this->nodeEdgea.size() << endl;
    cout << " PEC information:" << endl;
    cout << "  Boundary node 1 array exists (" << (this->bd_node1 != nullptr) << ")" << endl;
    cout << "  Boundary node 2 array exists (" << (this->bd_node2 != nullptr) << ")" << endl;
    cout << "  Boundary edge array exists (" << (this->bd_edge != nullptr) << ")" << endl;
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
    cout << "  Number of isolated conductors: " << this->numCdt << endl;
    if (this->patch == nullptr)
    {
        cout << "  Edge marker array exists (" << (this->markEdge != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Edge marker array has size " << NELEMENT(this->markEdge) << endl;
    }
    cout << "  Cell marker array exists (" << (this->markCell != nullptr) << ")" << endl;
    cout << "  cdtNumNode array exists (" << (this->cdtNumNode != nullptr) << ")" << endl;
    cout << "  sig array exist (" << (this->sig != nullptr) << ")" << endl;
    cout << "  conductor array exists (" << (this->conductor != nullptr) << ")" << endl;
    cout << "  Node marker array exists (" << (this->markNode != nullptr) << ")" << endl;
    cout << "  edgeCell vector has size " << this->edgeCell.size() << endl;
    cout << "  edgeCellArea vector has size " << this->edgeCellArea.size() << endl;
    cout << "  acu_cnno vector has size " << this->acu_cnno.size() << endl;
    cout << "  cindex vector has size " << this->cindex.size() << endl;
    cout << "  exciteCdtLayer array exists (" << (this->exciteCdtLayer != nullptr) << ")" << endl;
    cout << "  cond2condIn vector has size " << this->cond2condIn.size() << endl;
    cout << "  markProSide array exists (" << (this->markProSide != nullptr) << ")" << endl;
    if (this->patch == nullptr)
    {
        cout << " Patch information exists (" << (this->patch != nullptr) << ")" << endl;
    }
    else
    {
        cout << " Patch information has size " << NELEMENT(this->patch) << endl;
    }
    if (this->bound == nullptr)
    {
        cout << " Boundary information exists (" << (this->bound != nullptr) << ")" << endl;
    }
    else
    {
        cout << " Boundary information has size " << NELEMENT(this->bound) << endl;
    }
    cout << " V0c information:" << endl;
    if (this->v0cval == nullptr)
    {
        cout << "  v0cRowId array exists (" << (this->v0cRowId != nullptr) << ")" << endl;
        cout << "  v0cColId array exists (" << (this->v0cColId != nullptr) << ")" << endl;
        cout << "  v0cColIdo array exists (" << (this->v0cColIdo != nullptr) << ")" << endl;
        cout << "  v0cval array exists (" << (this->v0cval != nullptr) << ")" << endl;
        cout << "  v0cvalo array exists (" << (this->v0cvalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0cRowId array has size " << NELEMENT(this->v0cRowId) << endl;
        cout << "  v0cColId array has size " << NELEMENT(this->v0cColId) << endl;
        cout << "  v0cColIdo array has size " << NELEMENT(this->v0cColIdo) << endl;
        cout << "  v0cval array has size " << NELEMENT(this->v0cval) << endl;
        cout << "  v0cvalo array has size " << NELEMENT(this->v0cvalo) << endl;
    }
    if (this->v0caval == nullptr)
    {
        cout << "  v0caRowId array exists (" << (this->v0caRowId != nullptr) << ")" << endl;
        cout << "  v0caColId array exists (" << (this->v0caColId != nullptr) << ")" << endl;
        cout << "  v0caColIdo array exists (" << (this->v0caColIdo != nullptr) << ")" << endl;
        cout << "  v0caval array exists (" << (this->v0caval != nullptr) << ")" << endl;
        cout << "  v0cavalo array exists (" << (this->v0cavalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0caRowId array has size " << NELEMENT(this->v0caRowId) << endl;
        cout << "  v0caColId array has size " << NELEMENT(this->v0caColId) << endl;
        cout << "  v0caColIdo array has size " << NELEMENT(this->v0caColIdo) << endl;
        cout << "  v0caval array has size " << NELEMENT(this->v0caval) << endl;
        cout << "  v0cavalo array has size " << NELEMENT(this->v0cavalo) << endl;
    }
    cout << " Arrays from technical paper information:" << endl;
    cout << "  v0c2y0c2 array exists (" << (this->v0c2y0c2 != nullptr) << ")" << endl;
    cout << "  v0c2y0c2o array exists (" << (this->v0c2y0c2o != nullptr) << ")" << endl;
    cout << "  yc array exists (" << (this->yc != nullptr) << ")" << endl;
    cout << "  v0cy0c array exists (" << (this->v0cy0c != nullptr) << ")" << endl;
    cout << " V0c' * D_sig * V0c information:" << endl;
    if (this->Acval == nullptr)
    {
        cout << "  AcRowId array exists (" << (this->AcRowId != nullptr) << ")" << endl;
        cout << "  AcRowId1 array exists (" << (this->AcRowId1 != nullptr) << ")" << endl;
        cout << "  AcColId array exists (" << (this->AcColId != nullptr) << ")" << endl;
        cout << "  Acval array exists (" << (this->Acval != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  AcRowId array has size " << NELEMENT(this->AcRowId) << endl;
        cout << "  AcRowId1 array has size " << NELEMENT(this->AcRowId1) << endl;
        cout << "  AcColId array has size " << NELEMENT(this->AcColId) << endl;
        cout << "  Acval array has size " << NELEMENT(this->Acval) << endl;
    }
    if (this->Adval == nullptr)
    {
        cout << "  AdRowId array exists (" << (this->AdRowId != nullptr) << ")" << endl;
        cout << "  AdRowId1 array exists (" << (this->AdRowId1 != nullptr) << ")" << endl;
        cout << "  AdColId array exists (" << (this->AdColId != nullptr) << ")" << endl;
        cout << "  Adval array exists (" << (this->Adval != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  AdRowId array has size " << NELEMENT(this->AdRowId) << endl;
        cout << "  AdRowId1 array has size " << NELEMENT(this->AdRowId1) << endl;
        cout << "  AdColId array has size " << NELEMENT(this->AdColId) << endl;
        cout << "  Adval array has size " << NELEMENT(this->Adval) << endl;
    }
    cout << "  crhs array exists (" << (this->crhs != nullptr) << ")" << endl;
    cout << " V0d_ information:" << endl;
    if (this->v0d1val == nullptr)
    {
        cout << "  v0d1RowId array exists (" << (this->v0d1RowId != nullptr) << ")" << endl;
        cout << "  v0d1ColId array exists (" << (this->v0d1ColId != nullptr) << ")" << endl;
        cout << "  v0d1ColIdo array exists (" << (this->v0d1ColIdo != nullptr) << ")" << endl;
        cout << "  v0d1val array exists (" << (this->v0d1val != nullptr) << ")" << endl;
        cout << "  v0d1valo array exists (" << (this->v0d1valo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d1RowId array has size " << NELEMENT(this->v0d1RowId) << endl;
        cout << "  v0d1ColId array has size " << NELEMENT(this->v0d1ColId) << endl;
        cout << "  v0d1ColIdo array has size " << NELEMENT(this->v0d1ColIdo) << endl;
        cout << "  v0d1val array has size " << NELEMENT(this->v0d1val) << endl;
        cout << "  v0d1valo array has size " << NELEMENT(this->v0d1valo) << endl;
    }
    if (this->v0d1aval == nullptr)
    {
        cout << "  v0d1aRowId array exists (" << (this->v0d1aRowId != nullptr) << ")" << endl;
        cout << "  v0d1aColId array exists (" << (this->v0d1aColId != nullptr) << ")" << endl;
        cout << "  v0d1aColIdo array exists (" << (this->v0d1aColIdo != nullptr) << ")" << endl;
        cout << "  v0d1aval array exists (" << (this->v0d1aval != nullptr) << ")" << endl;
        cout << "  v0d1avalo array exists (" << (this->v0d1avalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d1aRowId array has size " << NELEMENT(this->v0d1aRowId) << endl;
        cout << "  v0d1aColId array has size " << NELEMENT(this->v0d1aColId) << endl;
        cout << "  v0d1aColIdo array has size " << NELEMENT(this->v0d1aColIdo) << endl;
        cout << "  v0d1aval array has size " << NELEMENT(this->v0d1aval) << endl;
        cout << "  v0d1avalo array has size " << NELEMENT(this->v0d1avalo) << endl;
    }
    if (this->v0d2val == nullptr)
    {
        cout << "  v0d2RowId array exists (" << (this->v0d2RowId != nullptr) << ")" << endl;
        cout << "  v0d2ColId array exists (" << (this->v0d2ColId != nullptr) << ")" << endl;
        cout << "  v0d2ColIdo array exists (" << (this->v0d2ColIdo != nullptr) << ")" << endl;
        cout << "  v0d2val array exists (" << (this->v0d2val != nullptr) << ")" << endl;
        cout << "  v0d2valo array exists (" << (this->v0d2valo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d2RowId array has size " << NELEMENT(this->v0d2RowId) << endl;
        cout << "  v0d2ColId array has size " << NELEMENT(this->v0d2ColId) << endl;
        cout << "  v0d2ColIdo array has size " << NELEMENT(this->v0d2ColIdo) << endl;
        cout << "  v0d2val array has size " << NELEMENT(this->v0d2val) << endl;
        cout << "  v0d2valo array has size " << NELEMENT(this->v0d2valo) << endl;
    }
    if (this->v0d2aval == nullptr)
    {
        cout << "  v0d2aRowId array exists (" << (this->v0d2aRowId != nullptr) << ")" << endl;
        cout << "  v0d2aColId array exists (" << (this->v0d2aColId != nullptr) << ")" << endl;
        cout << "  v0d2aColIdo array exists (" << (this->v0d2aColIdo != nullptr) << ")" << endl;
        cout << "  v0d2aval array exists (" << (this->v0d2aval != nullptr) << ")" << endl;
        cout << "  v0d2avalo array exists (" << (this->v0d2avalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d2aRowId array has size " << NELEMENT(this->v0d2aRowId) << endl;
        cout << "  v0d2aColId array has size " << NELEMENT(this->v0d2aColId) << endl;
        cout << "  v0d2aColIdo array has size " << NELEMENT(this->v0d2aColIdo) << endl;
        cout << "  v0d2aval array has size " << NELEMENT(this->v0d2aval) << endl;
        cout << "  v0d2avalo array has size " << NELEMENT(this->v0d2avalo) << endl;
    }
    cout << "  yd array exists (" << (this->yd != nullptr) << ")" << endl;
    cout << " S matrix information:" << endl;
    cout << "  SRowId array exists (" << (this->SRowId != nullptr) << ")" << endl;
    cout << "  SColId array exists (" << (this->SColId != nullptr) << ")" << endl;
    cout << "  Sval array exists (" << (this->Sval != nullptr) << ")" << endl;
    cout << "  Value array for S matrix has length " << this->leng_S << endl;
    cout << " Solution storage information:" << endl;
    if (this->y == nullptr)
    {
        cout << "  y array exists (" << (this->y != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  y array has size " << NELEMENT(this->y) << endl;
    }
    if (this->x == nullptr)
    {
        cout << "  x array exists (" << (this->x != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  x array has size " << NELEMENT(this->x) << endl;
    }
    cout << " Port information:" << endl;
    cout << "  Number of ports: " << this->numPorts << endl;
    if (this->portCoor == nullptr)
    {
        cout << "  portCoor array exists (" << (this->portCoor != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  portCoor array has size " << NELEMENT(this->portCoor) << endl;
        cout << "   Showing information about first port   vvvvvv" << endl;
        this->portCoor[0].print();
    }
    cout << "  portEdge vector has size " << this->portEdge.size() << endl;
    cout << "  portArea vector has size " << this->portArea.size() << endl;
    cout << "  portNno unordered set has values stored (" << (!(this->portNno.empty())) << ")" << endl;
    if (this->J == nullptr)
    {
        cout << " Current source array exists (" << (this->J != nullptr) << ")" << endl;
    }
    else
    {
        cout << " Current source array has size " << NELEMENT(this->J) << endl;
    }
    cout << " Current V0c,s^T*I information:" << endl;
    cout << "  v0csJ array exists (" << (this->v0csJ != nullptr) << ")" << endl;
    cout << "  Y array exists (" << (this->Y != nullptr) << ")" << endl;
    cout << "------" << endl;
}
=======
//#include "stdafx.h"
#include "fdtd.h"



int meshAndMark(fdtdMesh *sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi, unordered_set<double> *portCoorx, unordered_set<double> *portCoory)
{
    int lyr;
    myint i, j, k, m;
    double upper, lower;
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
    double disMin = MINDIS;
    double disMaxx, disMaxy;   // the max discretization in x, y, z directions
    for (i = 0; i < sys->numCdtRow; i++) {
        numNode += sys->conductorIn[i].numVert;
    }
    //cout << numNode << endl;
    xOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    yOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    zOrigOld = (double*)calloc(2 * sys->numStack + 2 * sys->numPorts, sizeof(double));

    j = 0;
    for (i = 0; i < sys->numCdtRow; i++) {
        for (k = 0; k < sys->conductorIn[i].numVert; k++) {
            xOrigOld[j] = sys->conductorIn[i].x[k];
            yOrigOld[j] = sys->conductorIn[i].y[k];
            j++;
        }
    }

    for (i = 0; i < sys->numPorts; i++) {
        xOrigOld[j] = sys->portCoor[i].x1;
        yOrigOld[j] = sys->portCoor[i].y1;
        j++;
        xOrigOld[j] = sys->portCoor[i].x2;
        yOrigOld[j] = sys->portCoor[i].y2;
        j++;
    }
    j = 0;
    for (i = 0; i < sys->numStack; i++) {
        zOrigOld[j] = sys->stackBegCoor[i];
        j++;
        zOrigOld[j] = sys->stackEndCoor[i];
        j++;
    }
    for (i = 0; i < sys->numPorts; i++) {
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
    int xMaxInd = 2;
    disMaxx = 1e-5;// (xmax - xmin) / xMaxInd;
    xMaxInd = (xmax - xmin) / disMaxx;

    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(xOrigOld[i] - xOrigOld[i - 1]) > disMin){
            sys->nx++;
        }
    }
    double *xn = (double*)calloc(numNode + 6 * sys->numPorts + xMaxInd, sizeof(double));
    xn[0] = xOrigOld[0];
    double temp = xn[0];
    j = 0;
    sys->nx = 1;
    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(xOrigOld[i] - temp) > disMin && abs(xOrigOld[i] - temp) <= disMaxx){
            j++;
            xn[j] = xOrigOld[i];
            temp = xn[j];
            sys->nx++;
        }
        else if (abs(xOrigOld[i] - temp) > disMin && abs(xOrigOld[i] - temp) > disMaxx){
            while (abs(xOrigOld[i] - temp) > disMaxx){
                j++;
                xn[j] = xn[j - 1] + disMaxx;
                temp = xn[j];
                sys->nx++;
            }
            if (abs(xOrigOld[i] - temp) > disMin){
                sys->nx++;
                temp = xOrigOld[i];
            }
            j++;
            xn[j] = xOrigOld[i];
        }
        else {
            j++;
            xn[j] = xOrigOld[i];
        }
    }
    int countx = j;
    //sort(xn, xn + countx + 1);


    //sys->xnu = (double*)calloc(sys->nx, sizeof(double));

    j = 0;
    //sys->xnu[0] = xn[0];
    //xi[sys->xnu[0]] = j;
    double first, second;
    temp = xn[0];
    for (i = 1; i <= countx; i++){    // set the discretization length around port to be equal
        if (abs(xn[i] - temp) > disMin){
            j++;
            temp = xn[i];
            //sys->xnu[j] = xn[i];
            //xi[sys->xnu[j]] = j;
        }
        else {
            //xi[xn[i]] = j;
        }
    }
    sys->nx = j + 1;




    /***************************************************************************/
    sort(yOrigOld, yOrigOld + numNode + 2 * sys->numPorts);
    
    sys->ny = 1;
    ymin = yOrigOld[0];
    ymax = yOrigOld[numNode + 2 * sys->numPorts - 1];
    int yMaxInd = 2;    // the max discretization of y is total / 120
    disMaxy = 1e-5;// (ymax - ymin) / yMaxInd;
    yMaxInd = (ymax - ymin) / disMaxy;

    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(yOrigOld[i] - yOrigOld[i - 1]) > disMin){
            sys->ny++;
        }
    }

    double *yn = (double*)calloc(numNode + 6 * sys->numPorts + yMaxInd, sizeof(double));
    yn[0] = yOrigOld[0];
    j = 0;
    sys->ny = 1;
    
    temp = yn[0];
    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(yOrigOld[i] - temp) > disMin && abs(yOrigOld[i] - temp) <= disMaxy){
            j++;
            yn[j] = yOrigOld[i];
            temp = yn[j];
            //cout << yn[j] << " ";
            sys->ny++;
        }
        else if (abs(yOrigOld[i] - temp) > disMin && abs(yOrigOld[i] - temp) > disMaxy){
            while (abs(yOrigOld[i] - temp) > disMaxy){
                j++;
                yn[j] = yn[j - 1] + disMaxy;
                temp = yn[j];
                //cout << yn[j] << " " << yn[j - 1] << "    ";
                sys->ny++;
            }
            if (abs(yOrigOld[i] - temp) > disMin){
                sys->ny++;
                temp = yOrigOld[i];
            }
            j++;
            yn[j] = yOrigOld[i];

            //cout << yn[j] << " ";
        }
        else {
            j++;
            yn[j] = yOrigOld[i];
            //cout << yn[j] << " ";
        }
    }
    //cout << endl;

    int county = j;
    //sort(yn, yn + county + 1);

    //sys->ynu = (double*)calloc(sys->ny, sizeof(double));


    j = 0;
    //sys->ynu[0] = yn[0];
    //yi[sys->ynu[0]] = j;
    temp = yn[0];
    for (i = 1; i <= county; i++){    // set the discretization length around port to be equal
        
        if (abs(yn[i] - temp) > disMin){
            j++;
            temp = yn[i];
            //sys->ynu[j] = yn[i];
            //yi[sys->ynu[j]] = j;
            
        }
        else {
            //yi[yn[i]] = j;
        }
    }

    sys->ny = j + 1;




    /********************************************************************************/
    sort(zOrigOld, zOrigOld + 2 * sys->numStack + 2 * sys->numPorts);
    sys->nz = 1;
    double disMinz = 1e-9;
    for (i = 1; i < 2 * sys->numStack + 2 * sys->numPorts; i++){
        if (abs(zOrigOld[i] - zOrigOld[i - 1]) > disMinz){
            sys->nz++;
        }
    }

    double *zn = (double*)calloc(2 * sys->numStack + 6 * sys->numPorts, sizeof(double));
    for (i = 0; i < 2 * sys->numStack + 2 * sys->numPorts; i++) {
        zn[i] = zOrigOld[i];
    }
    int countz = 2 * sys->numStack + 2 * sys->numPorts - 1;



    /*************************************************************************************/


    sort(xn, xn + countx + 1);
    xi.clear();
    sys->xn = (double*)calloc(sys->nx, sizeof(double));
    j = 0;
    sys->xn[0] = xn[0];
    temp = sys->xn[0];
    //xi[sys->xn[0]] = j;
    for (i = 1; i <= countx; i++){    // set the discretization length around port to be equal
        if (abs(xn[i] - temp) > disMin){
            j++;
            sys->xn[j] = xn[i];
            temp = sys->xn[j];
            xi[sys->xn[j]] = j;
        }
        else {
            xi[xn[i]] = j;
        }
    }
    sys->nx = j + 1;
    free(xn); xn = NULL;

    sort(yn, yn + county + 1);
    yi.clear();
    sys->yn = (double*)calloc(sys->ny, sizeof(double));
    j = 0;
    sys->yn[0] = yn[0];
    temp = sys->yn[0];
    yi[sys->yn[0]] = j;
    for (i = 1; i <= county; i++){    // set the discretization length around port to be equal

        if (abs(yn[i] - temp) > disMin){
            j++;
            sys->yn[j] = yn[i];
            temp = sys->yn[j];
            yi[sys->yn[j]] = j;
        }
        else {
            yi[yn[i]] = j;
        }
    }
    sys->ny = j + 1;
    free(yn); yn = NULL;

    sort(zn, zn + countz + 1);
    zi.clear();
    sys->zn = (double*)calloc(sys->nz, sizeof(double));
    j = 0;
    sys->zn[0] = zn[0];
    zi[sys->zn[0]] = j;
    for (i = 1; i <= countz; i++){    // set the discretization length around port to be equal
        if (abs(zn[i] - zn[i - 1]) > disMinz){
            j++;
            sys->zn[j] = zn[i];
            zi[sys->zn[j]] = j;
        }
        else {
            zi[zn[i]] = j;
        }
    }
    sys->nz = j + 1;
    free(zn); zn = NULL;

    sys->stackEpsn = (double*)calloc(sys->nz - 1, sizeof(double));
    i = 0;


    if (sys->stackBegCoor[0] == 0){
        j = 0;
        while (i < sys->nz - 1) {
            if ((sys->zn[i] + sys->zn[i + 1]) / 2 >= sys->stackBegCoor[j] && (sys->zn[i] + sys->zn[i + 1]) / 2 <= sys->stackEndCoor[j]) {
                sys->stackEpsn[i] = sys->stackEps[j];
                i++;
            }
            else {
                j++;
            }
        }
    }
    else {
        j = sys->numStack - 1;
        while (i < sys->nz - 1) {
            if ((sys->zn[i] + sys->zn[i + 1]) / 2 >= sys->stackBegCoor[j] && (sys->zn[i] + sys->zn[i + 1]) / 2 <= sys->stackEndCoor[j]) {
                sys->stackEpsn[i] = sys->stackEps[j];
                i++;
            }
            else {
                j--;
            }
        }
    }
    /*for (i = 0; i < sys->nz - 1; i++){
        cout << sys->stackEpsn[i] << endl;
        }*/
    /*for (i = 0; i < sys->nx; i++){
        cout << sys->xn[i] << " ";
    }
    cout << "\n" << endl;*/
    /*for (i = 0; i < sys->ny; i++){
        cout << sys->yn[i] << " ";
    }
    cout << "\n" << endl;*/
    /*for (i = 0; i < sys->nz; i++){
        cout << sys->zn[i] << " ";
    }
    cout << "\n" << endl;*/

    /***********************************************************************************************/

    sys->N_cell_x = sys->nx - (myint)1;
    sys->N_cell_y = sys->ny - (myint)1;
    sys->N_cell_z = sys->nz - (myint)1;


    sys->N_edge_s = sys->N_cell_y*(sys->N_cell_x + 1) + sys->N_cell_x*(sys->N_cell_y + 1);
    sys->N_edge_v = (sys->N_cell_x + 1)*(sys->N_cell_y + 1);
    sys->N_edge = sys->N_edge_s*(sys->N_cell_z + 1) + sys->N_edge_v*(sys->N_cell_z);

    sys->N_node_s = sys->N_edge_v;
    sys->N_node = (sys->N_cell_z + 1)*sys->N_node_s;

    sys->N_patch_s = sys->N_cell_x*sys->N_cell_y;
    sys->N_patch_v = (sys->N_cell_x + 1)*sys->N_cell_y + (sys->N_cell_y + 1)*sys->N_cell_x;
    sys->N_patch = sys->N_patch_s*(sys->N_cell_z + 1) + sys->N_patch_v*sys->N_cell_z;

    sys->markEdge = (myint*)calloc(sys->N_edge, sizeof(myint));   // mark which edge is inside the conductor
    sys->markNode = (myint*)calloc(sys->N_node, sizeof(myint));   // mark which node is inside the conductor

    /*cout << "N_edge = " << sys->N_edge << endl;
    cout << "N_node = " << sys->N_node << endl;
    cout << "N_cell_x = " << sys->N_cell_x << endl;
    cout << "N_cell_y = " << sys->N_cell_y << endl;
    cout << "N_cell_z = " << sys->N_cell_z << endl;*/
    double xc, yc;
    
    unordered_map<myint, myint> xrange;
    vector<myint> xcoorv;
    set<myint> xcoor;
    unordered_map<myint, set<myint>> xcoory;    // store the start coordinate of the range, the end can be checked from xrange
    myint ss, ee;
    myint y1, y2;
    myint l;
    clock_t tt = clock();
    int mark1;
    // Fast algorithm to find nodes inside conductors
    for (i = 0; i < sys->numCdtRow; i++){
        mark1 = 0;
        //cout << "Number of CdtRow is " << i << endl;
        for (j = 0; j < sys->conductorIn[i].numVert; j++){
            xcoor.insert(xi[sys->conductorIn[i].x[j]]);
        }


        for (auto xcoori = xcoor.begin(); xcoori != xcoor.end(); ++xcoori){
            xcoorv.push_back(*xcoori);
        }
        

        for (j = 0; j < xcoorv.size() - 1; j++){
            mark1 = 1;    // the x coordinates are more than 1
            xrange[xcoorv[j]] = xcoorv[j + 1];
            
        }
        if (xcoorv.size() == 1){    // If it has only one value
            xrange[xcoorv[0]] = xcoorv[0];
        }
        

        for (j = 0; j < sys->conductorIn[i].numVert - 1; j++){
            if (sys->conductorIn[i].y[j] == sys->conductorIn[i].y[j + 1]){
                if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[j + 1]){
                    ss = xi[sys->conductorIn[i].x[j]];
                    ee = xi[sys->conductorIn[i].x[j + 1]];
                    if (ss == ee && mark1 == 1){
                        continue;
                    }
                    while (xrange[ss] <= ee){
                        xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                        if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
                else if (sys->conductorIn[i].x[j] > sys->conductorIn[i].x[j + 1]){
                    ss = xi[sys->conductorIn[i].x[j + 1]];
                    ee = xi[sys->conductorIn[i].x[j]];
                    if (ss == ee && mark1 == 1){
                        continue;
                    }
                    while (xrange[ss] <= ee){
                        xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                        if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                            break;
                        }
                        ss = xrange[ss];
                    }
                }
            }
        }
        if (sys->conductorIn[i].y[j] == sys->conductorIn[i].y[0]){
            if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[0]){
                ss = xi[sys->conductorIn[i].x[j]];
                ee = xi[sys->conductorIn[i].x[0]];
                if (ss == ee && mark1 == 1){
                    continue;
                }
                while (xrange[ss] <= ee){
                    xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                    if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                        break;
                    }
                    ss = xrange[ss];
                }
            }
            else if (sys->conductorIn[i].x[j] > sys->conductorIn[i].x[0]){
                ss = xi[sys->conductorIn[i].x[0]];
                ee = xi[sys->conductorIn[i].x[j]];
                if (ss == ee && mark1 == 1){
                    continue;
                }
                while (xrange[ss] <= ee){
                    xcoory[ss].insert(yi[sys->conductorIn[i].y[j]]);
                    if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){
                        break;
                    }
                    ss = xrange[ss];
                }
            }
        }
        /*for (auto xcooryi : xcoory){
            cout << xcooryi.first << " : ";
            for (auto xcooryii : xcooryi.second){
                cout << xcooryii << " ";
            }
            cout << endl;
        }*/
        for (auto xcooryi : xcoory){
            mark = 0;
            for (j = xcooryi.first; j < xrange[xcooryi.first]; j++){
                mark = 0;
                for (auto xrangey : xcooryi.second){
                    mark++;
                    if (mark % 2 == 1){
                        y1 = xrangey;
                        continue;
                    }
                    else if (mark % 2 == 0){
                        y2 = xrangey;
                        
                        for (l = zi[sys->conductorIn[i].zmin]; l <= zi[sys->conductorIn[i].zmax]; l++){
                            for (k = y1; k < y2; k++){
                                sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] = i + 1;
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] = i + 1;
                                if (l != zi[sys->conductorIn[i].zmax]){
                                    sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                                }
                            }
                            sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                            sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] = i + 1;
                            if (l != zi[sys->conductorIn[i].zmax]){
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                            }
                        }
                        
                    }
                }
            }
            mark = 0;
            for (auto xrangey : xcooryi.second){
                mark++;
                if (mark % 2 == 1){
                    y1 = xrangey;
                    continue;
                }
                else if (mark % 2 == 0){
                    y2 = xrangey;
                    for (l = zi[sys->conductorIn[i].zmin]; l <= zi[sys->conductorIn[i].zmax]; l++){
                        for (k = y1; k < y2; k++){
                            sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                            sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] = i + 1;
                            if (l != zi[sys->conductorIn[i].zmax]){
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                            }
                        }
                        sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                        if (l != zi[sys->conductorIn[i].zmax]){
                            sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                        }
                    }
                }
            }
        }
        xrange.clear();
        xcoor.clear();
        xcoorv.clear();
        xcoory.clear();
    }
    cout << "The time to mark Edge and Node is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << endl;
    //for (i = 0; i < sys->numCdtRow; i++){
    //    //cout << polyIn((0 + 4.9e-7) / 2, -1.7e-7, sys, 2) << endl;
    //    numNode = (xi[sys->conductorIn[i].xmax] - xi[sys->conductorIn[i].xmin] + 1)
    //        *(yi[sys->conductorIn[i].ymax] - yi[sys->conductorIn[i].ymin] + 1)
    //        *(zi[sys->conductorIn[i].zmax] - zi[sys->conductorIn[i].zmin] + 1);
    //    //cout << sys->conductorIn[i].xmax << " " << sys->conductorIn[i].xmin << " " << sys->conductorIn[i].ymax << " " << sys->conductorIn[i].ymin << " " << sys->conductorIn[i].zmax << " " << sys->conductorIn[i].zmin << endl;
    //    sys->conductorIn[i].cdtInNode = (myint*)malloc(numNode*sizeof(myint));
    //    sys->conductorIn[i].numNode = 0;
    //    

    //    for (j = xi[sys->conductorIn[i].xmin]; j <= xi[sys->conductorIn[i].xmax]; j++){
    //        for (k = yi[sys->conductorIn[i].ymin]; k <= yi[sys->conductorIn[i].ymax]; k++){
    //            if (polyIn(sys->xn[j], sys->yn[k], sys, i)){
    //                for (m = zi[sys->conductorIn[i].zmin]; m < zi[sys->conductorIn[i].zmax]; m++){
    //                    //cout << sys->xn[j] << " " << sys->yn[k] << endl;
    //                    sys->conductorIn[i].cdtInNode[sys->conductorIn[i].numNode] = m*sys->N_node_s + (sys->N_cell_y + 1)*j + k;
    //                    sys->conductorIn[i].numNode++;
    //                    sys->markNode[m*sys->N_node_s + (sys->N_cell_y + 1)*j + k] = 1;
    //                    if (sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] == 0){   // set the z direction markEdge
    //                        sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] = i + 1;
    //                        //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k << endl;
    //                    }
    //                }
    //                sys->conductorIn[i].cdtInNode[sys->conductorIn[i].numNode] = m*sys->N_node_s + (sys->N_cell_y + 1)*j + k;
    //                sys->conductorIn[i].numNode++;
    //                sys->markNode[m*sys->N_node_s + (sys->N_cell_y + 1)*j + k] = 1;
    //            }

    //        }
    //    }
    //    for (j = xi[sys->conductorIn[i].xmin]; j <= xi[sys->conductorIn[i].xmax]; j++){   // set the y direction markEdge
    //        for (k = yi[sys->conductorIn[i].ymin]; k < yi[sys->conductorIn[i].ymax]; k++){
    //            for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
    //                xc = sys->xn[j];
    //                yc = (sys->yn[k] + sys->yn[k + 1]) / 2;
    //                if (polyIn(xc, yc, sys, i)){
    //                    sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k] = i + 1;
    //                    //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
    //                    //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
    //                }
    //            }
    //        }
    //    }
    //    for (j = yi[sys->conductorIn[i].ymin]; j <= yi[sys->conductorIn[i].ymax]; j++){    // set the x direction markEdge
    //        for (k = xi[sys->conductorIn[i].xmin]; k < xi[sys->conductorIn[i].xmax]; k++){
    //            for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
    //                xc = (sys->xn[k] + sys->xn[k + 1]) / 2;
    //                yc = sys->yn[j];
    //                //cout << xc << " " << yc << "" << i << " " << polyIn(xc, yc, sys, i) << endl;
    //                if (polyIn(xc, yc, sys, i)){
    //                    sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j] = i + 1;
    //                    //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
    //                    //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
    //                }
    //            }
    //        }
    //    }
    //}

    for (i = 0; i < sys->N_edge_s; i++){    // the lower plane
        sys->markEdge[i] = sys->numCdtRow + 1;
    }
    for (i = 0; i < sys->N_node_s; i++){
        sys->markNode[i] = 1;
    }
    for (i = sys->N_edge - sys->N_edge_s; i < sys->N_edge; i++){    // the upper plane
        sys->markEdge[i] = sys->numCdtRow + 2;
    }
    for (i = sys->N_node - sys->N_node_s; i < sys->N_node; i++){
        sys->markNode[i] = 1;
    }

    /* construct edgelink */
    myint eno;
    sys->edgelink = (myint*)malloc(2 * sizeof(myint)*sys->N_edge);
    for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++) {
        for (i = 1; i <= sys->N_cell_x + 1; i++) {    //edge along y axis
            for (j = 1; j <= sys->N_cell_y; j++) {
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (i - 1)*sys->N_cell_y + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;

            }
        }
        for (i = 1; i <= sys->N_cell_x; i++) {    //edge along x axis
            for (j = 1; j <= sys->N_cell_y + 1; j++) {
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1)*sys->N_cell_y + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + i*(sys->N_cell_y + 1) + j - 1;

            }
        }
    }
    for (lyr = 1; lyr <= sys->N_cell_z; lyr++) {    // edge along z axis
        for (i = 1; i <= sys->N_cell_x + 1; i++) {
            for (j = 1; j <= sys->N_cell_y + 1; j++) {
                eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
                sys->edgelink[(eno - 1) * 2 + 2 - 1] = lyr*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
            }
        }
    }

    /* construct nodepos */
    myint nno;

    sys->nodepos = (double*)malloc(sizeof(double)*sys->N_node * 3);   //N_node rows and 3 columns, input row by row
    for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++) {
        for (i = 1; i <= sys->N_cell_x + 1; i++) {
            for (j = 1; j <= sys->N_cell_y + 1; j++) {
                nno = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;
                sys->nodepos[(nno - 1) * 3 + 1 - 1] = sys->xn[i - 1];
                sys->nodepos[(nno - 1) * 3 + 2 - 1] = sys->yn[j - 1];
                sys->nodepos[(nno - 1) * 3 + 3 - 1] = sys->zn[lyr - 1];
            }
        }
    }

    /* construct nodeEdge */
    double leng;
    vector<pair<myint, double> > a;
    for (i = 0; i < sys->N_node; i++) {
        sys->nodeEdge.push_back(a);
        sys->nodeEdgea.push_back(a);
    }
    for (i = 0; i < sys->N_edge; i++) {
        leng = pow((sys->nodepos[sys->edgelink[i * 2] * 3] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 1] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 1]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 2] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 2]), 2);
        leng = sqrt(leng);
        sys->nodeEdge[sys->edgelink[i * 2]].push_back(make_pair(i, 1 / leng));
        sys->nodeEdge[sys->edgelink[i * 2 + 1]].push_back(make_pair(i, -1 / leng));
    }
    cout << "nodeEdge is done\n";
    int ix, iy, iz;
    iz = sys->N_cell_z;
    for (ix = 0; ix < sys->nx; ix++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, -1 / (sys->zn[iz] - sys->zn[iz - 1])));
        }
    }
    iz = 0;
    for (ix = 0; ix < sys->nx; ix++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->zn[iz + 1] - sys->zn[iz])));
        }
    }
    for (iz = 1; iz < sys->N_cell_z; iz++) {
        for (ix = 0; ix < sys->nx; ix++) {
            for (iy = 0; iy < sys->ny; iy++) {
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, -2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
            }
        }
    }
    ix = sys->N_cell_x;
    for (iz = 0; iz < sys->nz; iz++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -1 / (sys->xn[ix] - sys->xn[ix - 1])));
        }
    }
    ix = 0;
    for (iz = 0; iz < sys->nz; iz++) {
        for (iy = 0; iy < sys->ny; iy++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->xn[ix + 1] - sys->xn[ix])));
        }
    }
    for (ix = 1; ix < sys->N_cell_x; ix++) {
        for (iz = 0; iz < sys->nz; iz++) {
            for (iy = 0; iy < sys->ny; iy++) {
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
            }
        }
    }
    iy = sys->N_cell_y;
    for (iz = 0; iz < sys->nz; iz++) {
        for (ix = 0; ix < sys->nx; ix++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -1 / (sys->yn[iy] - sys->yn[iy - 1])));
        }
    }
    iy = 0;
    for (iz = 0; iz < sys->nz; iz++) {
        for (ix = 0; ix < sys->nx; ix++) {
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 1 / (sys->yn[iy + 1] - sys->yn[iy])));
        }
    }
    for (iy = 1; iy < sys->N_cell_y; iy++) {
        for (iz = 0; iz < sys->nz; iz++) {
            for (ix = 0; ix < sys->nx; ix++) {
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
            }
        }
    }

    /* implement dfs */
    //cout <<"Number of nodes: " << sys->N_node << endl;
    myint *visited;
    vector<int> st;
    unordered_set<int> base;
    visited = (myint*)calloc(sys->N_node, sizeof(myint));
    count = 0;

    for (i = 0; i < sys->N_node; i++) {
        if (sys->markNode[i] == 0) {
            continue;
        }
        else {
            if (visited[i] != 0) {
                continue;
            }
            else {
                st.clear();
                st.push_back(i);
                count++;
                visited[i] = count;
                sys->cond2condIn.push_back(base);
                while (!st.empty()) {
                    mark = 0;
                    for (j = 0; j < sys->nodeEdge[st.back()].size(); j++) {
                        if (sys->markEdge[sys->nodeEdge[st.back()][j].first] != 0) {
                            if (sys->cond2condIn[count - 1].find(sys->markEdge[sys->nodeEdge[st.back()][j].first]) == sys->cond2condIn[count - 1].end()) {
                                sys->cond2condIn[count - 1].insert(sys->markEdge[sys->nodeEdge[st.back()][j].first]);
                            }
                            if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] == 0)) {
                                visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] = count;
                                st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]);
                                mark = 1;

                                break;
                            }
                            else if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] == 0)) {
                                visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] = count;
                                st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]);
                                mark = 1;

                                break;
                            }

                        }
                    }
                    if (mark == 0) {
                        st.pop_back();
                    }
                }
            }
        }
    }
    
    cout << "Number of conductors is " << count << endl;
    /*for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            for (j = 0; j < sys->nodeEdge[i].size(); j++){
                if (sys->markEdge[sys->nodeEdge[i][j].first] != 0 && sys->cond2condIn[visited[i] - 1].find(sys->markEdge[sys->nodeEdge[i][j].first]) == sys->cond2condIn[visited[i] - 1].end()){
                    sys->cond2condIn[visited[i] - 1].insert(sys->markEdge[sys->nodeEdge[i][j].first]);
                }
            }
        }
    }*/


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
    /*for (i = 0; i < sys->numCdt; i++) {
        cout << "Cnt " << i << " has condIn as: ";
        for (auto ci : sys->cond2condIn[i]) {
            cout << ci << " ";
        }
        cout << endl;
    }*/
    sys->conductor = (fdtdCdt*)malloc(sys->numCdt * sizeof(fdtdCdt));
    sys->cdtNumNode = (myint*)calloc(sys->numCdt, sizeof(myint));
    for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            sys->cdtNumNode[visited[i] - 1]++;
        }
    }
    for (i = 0; i < sys->numCdt; i++){
        sys->conductor[i].node = (myint*)malloc(sizeof(myint) * sys->cdtNumNode[i]);
        sys->conductor[i].cdtNodeind = 0;
        sys->conductor[i].markPort = 0;
    }
    int portground_count = 0;
    for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            sys->conductor[visited[i] - 1].node[sys->conductor[visited[i] - 1].cdtNodeind] = i;
            if ((i < sys->N_node_s) && sys->conductor[visited[i] - 1].markPort != -1){
                /*portground_count++;
                if (portground_count <= 1)*/
                    sys->conductor[visited[i] - 1].markPort = -1;    // this conductor connect to the lower PEC
            }
            else if ((i >= sys->N_node - sys->N_node_s) && sys->conductor[visited[i] - 1].markPort != -2){
                sys->conductor[visited[i] - 1].markPort = -2;    // this conductor connects to the upper PEC
            }
            sys->conductor[visited[i] - 1].cdtNodeind++;
        }
    }
    free(visited);
    visited = NULL;
    
    /* set markCell */
    //vector<int> aa;
    //vector<double> bb;
    //sys->markCell = (int*)calloc(sys->N_cell_x * sys->N_cell_y * sys->N_cell_z, sizeof(int));
    //for (i = 0; i < sys->N_edge; i++)
    //{
    //    sys->edgeCell.push_back(aa);
    //    sys->edgeCellArea.push_back(bb);
    //}
    //int cell;
    //for (i = 0; i < sys->N_cell_z; i++){
    //    for (j = 0; j < sys->N_cell_x; j++){
    //        for (k = 0; k < sys->N_cell_y; k++){
    //            mark = 1;
    //            cell = i * sys->N_patch_s + j * sys->N_cell_y + k;
    //            count = sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k];

    //            // y axis
    //            sys->edgeCell[(i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k)].push_back(cell);
    //            sys->edgeCellArea[(i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k)].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] == 0){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + j * (sys->N_cell_y) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (j + 1) * (sys->N_cell_y) + k] != count){
    //                mark = 0;
    //            }

    //            // x axis
    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back((sys->yn[k + 1] - sys->yn[k]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back((sys->yn[k + 1] - sys->yn[k])*(sys->zn[i + 1] - sys->zn[i])); 
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k].push_back((sys->yn[k + 1] - sys->yn[k])*(sys->zn[i + 1] - sys->zn[i])); 
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1].push_back((sys->yn[k + 1] - sys->yn[k]) * (sys->zn[i + 1] - sys->zn[i]));
    //            if (sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[(i + 1) * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y) * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            // z axis
    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j *(sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k] != count){
    //                mark = 0;
    //            }

    //            sys->edgeCell[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1].push_back(cell);
    //            sys->edgeCellArea[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1].push_back((sys->xn[j + 1] - sys->xn[j]) * (sys->yn[k + 1] - sys->yn[k]));
    //            if (sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1] == 0 || sys->markEdge[i * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (j + 1) *(sys->N_cell_y + 1) + k + 1] != count){
    //                mark = 0;
    //            }

    //            if (mark == 1){
    //                sys->markCell[cell] = 1;
    //            }
    //        }
    //    }
    //}
    


    return 0;
}

int matrixConstruction(fdtdMesh *sys){
    int i, j;

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

    /* construct D_sig */
    sys->sig = (double*)calloc(sys->N_edge, sizeof(double));
    double a, b;
    for (i = 0; i < sys->N_edge; i++){
        if (sys->markEdge[i] != 0){
            a = 0;
            b = 0;
            /*for (j = 0; j < sys->edgeCell[i].size(); j++){
                a += sys->markCell[sys->edgeCell[i][j]] * sys->edgeCellArea[i][j];
                b += sys->edgeCellArea[i][j];
            }*/
            sys->sig[i] = SIGMA;/*(a / b) * SIGMA;*/
        }
    }

    sys->edgeCell.clear();
    sys->edgeCellArea.clear();
    //free(sys->markCell); sys->markCell = NULL;

    return 0;
}

int portSet(fdtdMesh* sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi){
    char s[FDTD_MAXC];
    char *word[FDTD_MAXC];
    int j, i, mark, l, m, k;
    int node;
    vector<int> edge;
    double a;
    double sideLen = 0.;

    for (i = 0; i < sys->numPorts; i++)
    {
        
        
        if (sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] != 0 && sys->conductor[sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] - 1].markPort > -1){
            sys->portCoor[i].portCnd = sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]];
            sys->conductor[sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] - 1].markPort = i + 1;    // markPort start from 1
        }
        else if (sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] != 0 && sys->conductor[sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] - 1].markPort > -1){
            sys->portCoor[i].portCnd = sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]];
            sys->conductor[sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] - 1].markPort = i + 1;    // markPort start from 1
        }
        /*for (j = 0; j < sys->cdtNumNode[sys->portCoor[i].portCnd - 1]; j++){
            sys->exciteCdtLayer[sys->conductor[sys->portCoor[i].portCnd - 1].node[j] / sys->N_node_s] = 1;
        }*/
        //cout << sys->portCoor[i].portCnd << endl;
        edge.clear();
        if (sys->portCoor[i].x1 != sys->portCoor[i].x2){
            if (sys->portCoor[i].x1 < sys->portCoor[i].x2){
                for (j = xi[sys->portCoor[i].x1]; j < xi[sys->portCoor[i].x2]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*(sys->N_cell_x + 1) + j*(sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]);
                }
            }
            else{
                for (j = xi[sys->portCoor[i].x2]; j < xi[sys->portCoor[i].x1]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*(sys->N_cell_x + 1) + j*(sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]);
                }
            }
            a = 1;
            if (yi[sys->portCoor[i].y1] == 0){
                a *= (sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1);
            }
            else if (yi[sys->portCoor[i].y1] == sys->N_cell_y){
                a *= (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]);
            }
            else{
                a *= ((sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1) / 2 + (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]) / 2);
            }
            if (zi[sys->portCoor[i].z1] == 0){
                a *=(sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
            }
            else if (zi[sys->portCoor[i].z1] == sys->N_cell_z){
                a *= (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]);
            }
            else{
                a *= ((sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1) / 2 + (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]) / 2);
            }
            sys->portArea.push_back(a);
        }
        else if (sys->portCoor[i].y1 != sys->portCoor[i].y2){
            if (sys->portCoor[i].y1 < sys->portCoor[i].y2){
                for (j = yi[sys->portCoor[i].y1]; j < yi[sys->portCoor[i].y2]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*xi[sys->portCoor[i].x1] + j);
                }
            }
            else{
                for (j = yi[sys->portCoor[i].y2]; j < yi[sys->portCoor[i].y1]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*zi[sys->portCoor[i].z1] + (sys->N_cell_y)*xi[sys->portCoor[i].x1] + j);
                }
            }
            a = 1;
            if (xi[sys->portCoor[i].x1] == 0){
                a *= (sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1);
            }
            else if (xi[sys->portCoor[i].x1] == sys->N_cell_x){
                a *= (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]);
            }
            else{
                a *= ((sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1) / 2 + (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]) / 2);
            }
            if (zi[sys->portCoor[i].z1] == 0){
                a *=(sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
            }
            else if (zi[sys->portCoor[i].z1] == sys->N_cell_z){
                a *= (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]);
            }
            else{
                a *= ((sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1) / 2 + (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]) / 2);
            }
            sys->portArea.push_back(a);
        }
        else if (sys->portCoor[i].z1 != sys->portCoor[i].z2){
            if (sys->portCoor[i].z1 < sys->portCoor[i].z2){
                for (j = zi[sys->portCoor[i].z1]; j < zi[sys->portCoor[i].z2]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*j + sys->N_edge_s + (sys->N_cell_y + 1)*xi[sys->portCoor[i].x1] + yi[sys->portCoor[i].y1]);
                }
            }
            else{
                for (j = zi[sys->portCoor[i].z2]; j < zi[sys->portCoor[i].z1]; j++){
                    edge.push_back((sys->N_edge_s + sys->N_edge_v)*j + sys->N_edge_s + (sys->N_cell_y + 1)*xi[sys->portCoor[i].x1] + yi[sys->portCoor[i].y1]);
                }
            }
            a = 1;
            if (xi[sys->portCoor[i].x1] == 0){
                a *= (sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1);
            }
            else if (xi[sys->portCoor[i].x1] == sys->N_cell_x){
                a *= (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]);
            }
            else{
                a *= ((sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1) / 2 + (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]) / 2);
            }
            if (yi[sys->portCoor[i].y1] == 0){
                a *=(sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1);
            }
            else if (yi[sys->portCoor[i].y1] == sys->N_cell_y){
                a *= (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]);
            }
            else{
                a *= ((sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1) / 2 + (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]) / 2);
            }
            sys->portArea.push_back(a);

        }
        sys->portEdge.push_back(edge);

    }

    clock_t t1 = clock();
    sys->markProSide = (int*)calloc(sys->N_node, sizeof(int));
    double x1, x2, y1, y2;
    myint x1_ind, x2_ind, y1_ind, y2_ind, z1_ind, z2_ind;

    for (i = 0; i < sys->numPorts; i++){
        /*cout << "Value of sys->portCoor[i].portCnd - 1: " << sys->portCoor[i].portCnd - 1 << endl;
        cout << "Size of the sys->cond2condIn[sys->portCoor[i].portCnd - 1] unordered_set: " << sys->cond2condIn[sys->portCoor[i].portCnd - 1].size() << endl; */
        for (auto ci : sys->cond2condIn[sys->portCoor[i].portCnd - 1]){
            for (l = 0; l < sys->conductorIn[ci - 1].numVert - 1; l++){
                if (sys->conductorIn[ci - 1].x[l] == sys->conductorIn[ci - 1].x[l + 1]){
                    x1 = sys->conductorIn[ci - 1].x[l];
                    x2 = x1;
                    if (sys->conductorIn[ci - 1].y[l] < sys->conductorIn[ci - 1].y[l + 1]){
                        y1 = sys->conductorIn[ci - 1].y[l];
                        y2 = sys->conductorIn[ci - 1].y[l + 1];
                    }
                    else{
                        y1 = sys->conductorIn[ci - 1].y[l + 1];
                        y2 = sys->conductorIn[ci - 1].y[l];
                    }
                    x1_ind = xi[x1];
                    x2_ind = xi[x2];
                    y1_ind = yi[y1];
                    y2_ind = yi[y2];
                    while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                        x1_ind--;
                    }
                    x1_ind++;
                    while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                        x2_ind++;
                    }
                    x2_ind--;
                    while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                        y1_ind--;
                    }
                    y1_ind++;
                    while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                        y2_ind++;
                    }
                    y2_ind--;
                    z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                    z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                    for (k = z1_ind; k <= z2_ind; k++){
                        for (j = x1_ind; j <= x2_ind; j++){
                            for (m = y1_ind; m <= y2_ind; m++){
                                if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                    sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                            }
                        }
                    }
                }
                else{
                    y1 = sys->conductorIn[ci - 1].y[l];
                    y2 = y1;
                    if (sys->conductorIn[ci - 1].x[l] < sys->conductorIn[ci - 1].x[l + 1]){
                        x1 = sys->conductorIn[ci - 1].x[l];
                        x2 = sys->conductorIn[ci - 1].x[l + 1];
                    }
                    else{
                        x1 = sys->conductorIn[ci - 1].x[l + 1];
                        x2 = sys->conductorIn[ci - 1].x[l];
                    }
                    x1_ind = xi[x1];
                    x2_ind = xi[x2];
                    y1_ind = yi[y1];
                    y2_ind = yi[y2];
                    while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                        x1_ind--;
                    }
                    x1_ind++;
                    while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                        x2_ind++;
                    }
                    x2_ind--;
                    while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                        y1_ind--;
                    }
                    y1_ind++;
                    while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                        y2_ind++;
                    }
                    y2_ind--;
                    z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                    z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                    for (k = z1_ind; k <= z2_ind; k++){
                        for (j = x1_ind; j <= x2_ind; j++){
                            for (m = y1_ind; m <= y2_ind; m++){
                                if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                    sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                            }
                        }
                    }
                }
                
            }
            if (sys->conductorIn[ci - 1].x[l] == sys->conductorIn[ci - 1].x[0]){
                x1 = sys->conductorIn[ci - 1].x[l];
                x2 = x1;
                if (sys->conductorIn[ci - 1].y[l] < sys->conductorIn[ci - 1].y[0]){
                    y1 = sys->conductorIn[ci - 1].y[l];
                    y2 = sys->conductorIn[ci - 1].y[0];
                }
                else{
                    y1 = sys->conductorIn[ci - 1].y[0];
                    y2 = sys->conductorIn[ci - 1].y[l];
                }
                x1_ind = xi[x1];
                x2_ind = xi[x2];
                y1_ind = yi[y1];
                y2_ind = yi[y2];
                while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                    x1_ind--;
                }
                x1_ind++;
                while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                    x2_ind++;
                }
                x2_ind--;
                while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                    y1_ind--;
                }
                y1_ind++;
                while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                    y2_ind++;
                }
                y2_ind--;
                z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                for (k = z1_ind; k <= z2_ind; k++){
                    for (j = x1_ind; j <= x2_ind; j++){
                        for (m = y1_ind; m <= y2_ind; m++){
                            if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                        }
                    }
                }
            }
            else{
                y1 = sys->conductorIn[ci - 1].y[l];
                y2 = y1;
                if (sys->conductorIn[ci - 1].x[l] < sys->conductorIn[ci - 1].x[0]){
                    x1 = sys->conductorIn[ci - 1].x[l];
                    x2 = sys->conductorIn[ci - 1].x[0];
                }
                else{
                    x1 = sys->conductorIn[ci - 1].x[0];
                    x2 = sys->conductorIn[ci - 1].x[l];
                }
                x1_ind = xi[x1];
                x2_ind = xi[x2];
                y1_ind = yi[y1];
                y2_ind = yi[y2];
                while (x1 - sys->xn[x1_ind] <= sideLen && x1_ind >= 0){
                    x1_ind--;
                }
                x1_ind++;
                while (sys->xn[x2_ind] - x2 <= sideLen && x2_ind < sys->nx){
                    x2_ind++;
                }
                x2_ind--;
                while (y1 - sys->yn[y1_ind] <= sideLen && y1_ind >= 0){
                    y1_ind--;
                }
                y1_ind++;
                while (sys->yn[y2_ind] - y2 <= sideLen && y2_ind < sys->ny){
                    y2_ind++;
                }
                y2_ind--;
                z1_ind = zi[sys->conductorIn[ci - 1].zmin];
                z2_ind = zi[sys->conductorIn[ci - 1].zmax];
                for (k = z1_ind; k <= z2_ind; k++){
                    for (j = x1_ind; j <= x2_ind; j++){
                        for (m = y1_ind; m <= y2_ind; m++){
                            if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0)
                                sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                        }
                    }
                }
            }
        }
    }
    cout << "Time of finding side nodes is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;

    sys->conductorIn.clear();
    return 0;
}

// Is point (x,y) within the polygon?
bool polyIn(double x, double y, fdtdMesh *sys, int inPoly){
    int npol;
    int i, j, k;
    bool isCond = false;
    double disMin = 1.e-10;

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

// Print fdtdPort information
void fdtdPort::print()
{
    // Print
    cout << " ------" << endl;
    cout << " Contents of fdtdPort" << endl;
    if (this->x == nullptr)
    {
        cout << "  x array exists (" << (this->x != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  x array has size " << NELEMENT(this->x) << endl;
    }
    if (this->y == nullptr)
    {
        cout << "  y array exists (" << (this->y != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  y array has size " << NELEMENT(this->y) << endl;
    }
    if (this->z == nullptr)
    {
        cout << "  z array exists (" << (this->z != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  z array has size " << NELEMENT(this->z) << endl;
    }
    cout << "  Coordinates of two points for the source:" << endl;
    cout << "   x-coordinates: " << this->x1 << " m and " << this->x2 << " m" << endl;
    cout << "   y-coordinates: " << this->y1 << " m and " << this->y2 << " m" << endl;
    cout << "   z-coordinates: " << this->z1 << " m and " << this->z2 << " m" << endl;
    cout << "  Port direction: " << this->portDirection << endl;
    if (this->node == nullptr)
    {
        cout << "  Node array exists (" << (this->node != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Node array has size " << NELEMENT(this->node) << endl;
    }
    cout << "  Number of nodes: " << this->nodenum << endl;
    cout << "  portCnd: " << this->portCnd << endl;
    cout << " ------" << endl;
}

// Print fdtdMesh information
void fdtdMesh::print()
{
    // Print
    cout << "------" << endl;
    cout << "Contents of fdtdMesh" << endl;
    cout << " Length unit: " << this->lengthUnit << " m" << endl;
    cout << " Frequency sweep parameters: " << endl;
    cout << "  Frequency unit: " << this->freqUnit << " Hz" << endl;
    cout << "  Starting frequency: " << this->freqStart << " * " << this->freqUnit << " Hz" << endl;
    cout << "  Ending frequency: " << this->freqEnd << " * " << this->freqUnit << " Hz" << endl;
    cout << "  Number of frequencies: " << this->nfreq << endl;
    cout << "  Frequency scaling: " << this->freqScale << endl;
    cout << " Node coordinate information:" << endl;
    cout << "  Number of nodes along direction: " << this->nx << " in x-dir, " << this->ny << " in y-dir, and " << this->nz << " in z-dir" << endl;
    cout << "  Node arrays: xn exists (" << (this->xn != nullptr) << "), yn exists (" << (this->yn != nullptr) << "), and zn exists (" << (this->zn != nullptr) << ")" << endl;
    cout << "  Node arrays: xnu exists (" << (this->xnu != nullptr) << "), ynu exists (" << (this->ynu != nullptr) << "), and znu exists (" << (this->znu != nullptr) << ")" << endl;
    cout << " Mesh cell information:" << endl;
    cout << "  Number in each direction: " << this->N_cell_x << " in x-dir, " << this->N_cell_y << " in y-dir, and " << this->N_cell_z << " in z-dir" << endl;
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
    cout << "  nodeEdgea vector has size " << this->nodeEdgea.size() << endl;
    cout << " PEC information:" << endl;
    cout << "  Boundary node 1 array exists (" << (this->bd_node1 != nullptr) << ")" << endl;
    cout << "  Boundary node 2 array exists (" << (this->bd_node2 != nullptr) << ")" << endl;
    cout << "  Boundary edge array exists (" << (this->bd_edge != nullptr) << ")" << endl;
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
    cout << "  Number of isolated conductors: " << this->numCdt << endl;
    if (this->patch == nullptr)
    {
        cout << "  Edge marker array exists (" << (this->markEdge != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Edge marker array has size " << NELEMENT(this->markEdge) << endl;
    }
    cout << "  Cell marker array exists (" << (this->markCell != nullptr) << ")" << endl;
    cout << "  cdtNumNode array exists (" << (this->cdtNumNode != nullptr) << ")" << endl;
    cout << "  sig array exist (" << (this->sig != nullptr) << ")" << endl;
    cout << "  conductor array exists (" << (this->conductor != nullptr) << ")" << endl;
    cout << "  Node marker array exists (" << (this->markNode != nullptr) << ")" << endl;
    cout << "  edgeCell vector has size " << this->edgeCell.size() << endl;
    cout << "  edgeCellArea vector has size " << this->edgeCellArea.size() << endl;
    cout << "  acu_cnno vector has size " << this->acu_cnno.size() << endl;
    cout << "  cindex vector has size " << this->cindex.size() << endl;
    cout << "  exciteCdtLayer array exists (" << (this->exciteCdtLayer != nullptr) << ")" << endl;
    cout << "  cond2condIn vector has size " << this->cond2condIn.size() << endl;
    cout << "  markProSide array exists (" << (this->markProSide != nullptr) << ")" << endl;
    if (this->patch == nullptr)
    {
        cout << " Patch information exists (" << (this->patch != nullptr) << ")" << endl;
    }
    else
    {
        cout << " Patch information has size " << NELEMENT(this->patch) << endl;
    }
    if (this->bound == nullptr)
    {
        cout << " Boundary information exists (" << (this->bound != nullptr) << ")" << endl;
    }
    else
    {
        cout << " Boundary information has size " << NELEMENT(this->bound) << endl;
    }
    cout << " V0c information:" << endl;
    if (this->v0cval == nullptr)
    {
        cout << "  v0cRowId array exists (" << (this->v0cRowId != nullptr) << ")" << endl;
        cout << "  v0cColId array exists (" << (this->v0cColId != nullptr) << ")" << endl;
        cout << "  v0cColIdo array exists (" << (this->v0cColIdo != nullptr) << ")" << endl;
        cout << "  v0cval array exists (" << (this->v0cval != nullptr) << ")" << endl;
        cout << "  v0cvalo array exists (" << (this->v0cvalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0cRowId array has size " << NELEMENT(this->v0cRowId) << endl;
        cout << "  v0cColId array has size " << NELEMENT(this->v0cColId) << endl;
        cout << "  v0cColIdo array has size " << NELEMENT(this->v0cColIdo) << endl;
        cout << "  v0cval array has size " << NELEMENT(this->v0cval) << endl;
        cout << "  v0cvalo array has size " << NELEMENT(this->v0cvalo) << endl;
    }
    if (this->v0caval == nullptr)
    {
        cout << "  v0caRowId array exists (" << (this->v0caRowId != nullptr) << ")" << endl;
        cout << "  v0caColId array exists (" << (this->v0caColId != nullptr) << ")" << endl;
        cout << "  v0caColIdo array exists (" << (this->v0caColIdo != nullptr) << ")" << endl;
        cout << "  v0caval array exists (" << (this->v0caval != nullptr) << ")" << endl;
        cout << "  v0cavalo array exists (" << (this->v0cavalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0caRowId array has size " << NELEMENT(this->v0caRowId) << endl;
        cout << "  v0caColId array has size " << NELEMENT(this->v0caColId) << endl;
        cout << "  v0caColIdo array has size " << NELEMENT(this->v0caColIdo) << endl;
        cout << "  v0caval array has size " << NELEMENT(this->v0caval) << endl;
        cout << "  v0cavalo array has size " << NELEMENT(this->v0cavalo) << endl;
    }
    cout << " Arrays from technical paper information:" << endl;
    cout << "  v0c2y0c2 array exists (" << (this->v0c2y0c2 != nullptr) << ")" << endl;
    cout << "  v0c2y0c2o array exists (" << (this->v0c2y0c2o != nullptr) << ")" << endl;
    cout << "  yc array exists (" << (this->yc != nullptr) << ")" << endl;
    cout << "  v0cy0c array exists (" << (this->v0cy0c != nullptr) << ")" << endl;
    cout << " V0c' * D_sig * V0c information:" << endl;
    if (this->Acval == nullptr)
    {
        cout << "  AcRowId array exists (" << (this->AcRowId != nullptr) << ")" << endl;
        cout << "  AcRowId1 array exists (" << (this->AcRowId1 != nullptr) << ")" << endl;
        cout << "  AcColId array exists (" << (this->AcColId != nullptr) << ")" << endl;
        cout << "  Acval array exists (" << (this->Acval != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  AcRowId array has size " << NELEMENT(this->AcRowId) << endl;
        cout << "  AcRowId1 array has size " << NELEMENT(this->AcRowId1) << endl;
        cout << "  AcColId array has size " << NELEMENT(this->AcColId) << endl;
        cout << "  Acval array has size " << NELEMENT(this->Acval) << endl;
    }
    if (this->Adval == nullptr)
    {
        cout << "  AdRowId array exists (" << (this->AdRowId != nullptr) << ")" << endl;
        cout << "  AdRowId1 array exists (" << (this->AdRowId1 != nullptr) << ")" << endl;
        cout << "  AdColId array exists (" << (this->AdColId != nullptr) << ")" << endl;
        cout << "  Adval array exists (" << (this->Adval != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  AdRowId array has size " << NELEMENT(this->AdRowId) << endl;
        cout << "  AdRowId1 array has size " << NELEMENT(this->AdRowId1) << endl;
        cout << "  AdColId array has size " << NELEMENT(this->AdColId) << endl;
        cout << "  Adval array has size " << NELEMENT(this->Adval) << endl;
    }
    cout << "  crhs array exists (" << (this->crhs != nullptr) << ")" << endl;
    cout << " V0d_ information:" << endl;
    if (this->v0d1val == nullptr)
    {
        cout << "  v0d1RowId array exists (" << (this->v0d1RowId != nullptr) << ")" << endl;
        cout << "  v0d1ColId array exists (" << (this->v0d1ColId != nullptr) << ")" << endl;
        cout << "  v0d1ColIdo array exists (" << (this->v0d1ColIdo != nullptr) << ")" << endl;
        cout << "  v0d1val array exists (" << (this->v0d1val != nullptr) << ")" << endl;
        cout << "  v0d1valo array exists (" << (this->v0d1valo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d1RowId array has size " << NELEMENT(this->v0d1RowId) << endl;
        cout << "  v0d1ColId array has size " << NELEMENT(this->v0d1ColId) << endl;
        cout << "  v0d1ColIdo array has size " << NELEMENT(this->v0d1ColIdo) << endl;
        cout << "  v0d1val array has size " << NELEMENT(this->v0d1val) << endl;
        cout << "  v0d1valo array has size " << NELEMENT(this->v0d1valo) << endl;
    }
    if (this->v0d1aval == nullptr)
    {
        cout << "  v0d1aRowId array exists (" << (this->v0d1aRowId != nullptr) << ")" << endl;
        cout << "  v0d1aColId array exists (" << (this->v0d1aColId != nullptr) << ")" << endl;
        cout << "  v0d1aColIdo array exists (" << (this->v0d1aColIdo != nullptr) << ")" << endl;
        cout << "  v0d1aval array exists (" << (this->v0d1aval != nullptr) << ")" << endl;
        cout << "  v0d1avalo array exists (" << (this->v0d1avalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d1aRowId array has size " << NELEMENT(this->v0d1aRowId) << endl;
        cout << "  v0d1aColId array has size " << NELEMENT(this->v0d1aColId) << endl;
        cout << "  v0d1aColIdo array has size " << NELEMENT(this->v0d1aColIdo) << endl;
        cout << "  v0d1aval array has size " << NELEMENT(this->v0d1aval) << endl;
        cout << "  v0d1avalo array has size " << NELEMENT(this->v0d1avalo) << endl;
    }
    if (this->v0d2val == nullptr)
    {
        cout << "  v0d2RowId array exists (" << (this->v0d2RowId != nullptr) << ")" << endl;
        cout << "  v0d2ColId array exists (" << (this->v0d2ColId != nullptr) << ")" << endl;
        cout << "  v0d2ColIdo array exists (" << (this->v0d2ColIdo != nullptr) << ")" << endl;
        cout << "  v0d2val array exists (" << (this->v0d2val != nullptr) << ")" << endl;
        cout << "  v0d2valo array exists (" << (this->v0d2valo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d2RowId array has size " << NELEMENT(this->v0d2RowId) << endl;
        cout << "  v0d2ColId array has size " << NELEMENT(this->v0d2ColId) << endl;
        cout << "  v0d2ColIdo array has size " << NELEMENT(this->v0d2ColIdo) << endl;
        cout << "  v0d2val array has size " << NELEMENT(this->v0d2val) << endl;
        cout << "  v0d2valo array has size " << NELEMENT(this->v0d2valo) << endl;
    }
    if (this->v0d2aval == nullptr)
    {
        cout << "  v0d2aRowId array exists (" << (this->v0d2aRowId != nullptr) << ")" << endl;
        cout << "  v0d2aColId array exists (" << (this->v0d2aColId != nullptr) << ")" << endl;
        cout << "  v0d2aColIdo array exists (" << (this->v0d2aColIdo != nullptr) << ")" << endl;
        cout << "  v0d2aval array exists (" << (this->v0d2aval != nullptr) << ")" << endl;
        cout << "  v0d2avalo array exists (" << (this->v0d2avalo != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  v0d2aRowId array has size " << NELEMENT(this->v0d2aRowId) << endl;
        cout << "  v0d2aColId array has size " << NELEMENT(this->v0d2aColId) << endl;
        cout << "  v0d2aColIdo array has size " << NELEMENT(this->v0d2aColIdo) << endl;
        cout << "  v0d2aval array has size " << NELEMENT(this->v0d2aval) << endl;
        cout << "  v0d2avalo array has size " << NELEMENT(this->v0d2avalo) << endl;
    }
    cout << "  yd array exists (" << (this->yd != nullptr) << ")" << endl;
    cout << " S matrix information:" << endl;
    cout << "  SRowId array exists (" << (this->SRowId != nullptr) << ")" << endl;
    cout << "  SColId array exists (" << (this->SColId != nullptr) << ")" << endl;
    cout << "  Sval array exists (" << (this->Sval != nullptr) << ")" << endl;
    cout << "  Value array for S matrix has length " << this->leng_S << endl;
    cout << " Solution storage information:" << endl;
    if (this->y == nullptr)
    {
        cout << "  y array exists (" << (this->y != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  y array has size " << NELEMENT(this->y) << endl;
    }
    if (this->x == nullptr)
    {
        cout << "  x array exists (" << (this->x != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  x array has size " << NELEMENT(this->x) << endl;
    }
    cout << " Port information:" << endl;
    cout << "  Number of ports: " << this->numPorts << endl;
    if (this->portCoor == nullptr)
    {
        cout << "  portCoor array exists (" << (this->portCoor != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  portCoor array has size " << NELEMENT(this->portCoor) << endl;
        cout << "   Showing information about first port   vvvvvv" << endl;
        this->portCoor[0].print();
    }
    cout << "  portEdge vector has size " << this->portEdge.size() << endl;
    cout << "  portArea vector has size " << this->portArea.size() << endl;
    cout << "  portNno unordered set has values stored (" << (!(this->portNno.empty())) << ")" << endl;
    if (this->J == nullptr)
    {
        cout << " Current source array exists (" << (this->J != nullptr) << ")" << endl;
    }
    else
    {
        cout << " Current source array has size " << NELEMENT(this->J) << endl;
    }
    cout << " Current V0c,s^T*I information:" << endl;
    cout << "  v0csJ array exists (" << (this->v0csJ != nullptr) << ")" << endl;
    cout << "  Y array exists (" << (this->Y != nullptr) << ")" << endl;
    cout << "------" << endl;
}
>>>>>>> 6f7fa1339b5be679bb5a99c8374394944978b561
