//#include "stdafx.h"
#include "fdtd.hpp"


int meshAndMark(fdtdMesh *sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi, unordered_set<double> *portCoorx, unordered_set<double> *portCoory)
{
    int lyr;
    myint i, j, k, m;
    double upper, lower;
    int status;
    int count;
    int mark;
    double xmin, xmax, ymin, ymax, xwid, ywid;
    clock_t tt = clock();

    /* Generate the mesh nodes based on conductorIn information */
    int numNode = 0;
    double *xOrigOld, *yOrigOld, *zOrigOld;

    double disMin = MINDISFRACXY * fmin(sys->xlim2 - sys->xlim1, sys->ylim2 - sys->ylim1) * sys->lengthUnit; // Minimum discretization retained in x- or y-directions after node merging is fraction of smaller of x-extent or y-extent

    for (i = 0; i < sys->numCdtRow; i++) {
        numNode += sys->conductorIn[i].numVert;
    }
    //cout << numNode << endl;
    xOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    yOrigOld = (double*)calloc(numNode + 2 * sys->numPorts, sizeof(double));
    zOrigOld = (double*)calloc(2 * sys->numStack + 2 * sys->numPorts, sizeof(double));
    double minLayerDist = sys->zlim2 - sys->zlim1;// Initialize smallest distance between layers as entire domain height (units included)

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
        if (sys->stackEndCoor[i] - sys->stackBegCoor[i] > 0)
            minLayerDist = fmin(minLayerDist, sys->stackEndCoor[i] - sys->stackBegCoor[i]); // Update smallest distance between layers as each layer processed (units included)
    }
    for (i = 0; i < sys->numPorts; i++) {
        zOrigOld[j] = sys->portCoor[i].z1;
        j++;
        zOrigOld[j] = sys->portCoor[i].z2;
        j++;
    }
    
    /*******************************************************************************************/
    /* Discretize domain in the x-direction */
    sort(xOrigOld, xOrigOld + numNode + 2 * sys->numPorts);
    sys->nx = 1;
    xmin = xOrigOld[0];
    xmax = xOrigOld[numNode + 2 * sys->numPorts - 1];
    double disMaxx = MAXDISFRACX * (sys->xlim2 - sys->xlim1) * sys->lengthUnit; // Maximum discretization distance in x-direction is fraction of x-extent
    int xMaxInd = (xmax - xmin) / disMaxx;

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
    /* Discretize domain in the y-direction */
    sort(yOrigOld, yOrigOld + numNode + 2 * sys->numPorts);
    sys->ny = 1;
    ymin = yOrigOld[0];
    ymax = yOrigOld[numNode + 2 * sys->numPorts - 1];
    double disMaxy = MAXDISFRACY * (sys->ylim2 - sys->ylim1) * sys->lengthUnit; // Maximum discretization distance in y-direction is fraction of y-extent
    int yMaxInd = (ymax - ymin) / disMaxy;

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
            sys->ny++;
        }
        else if (abs(yOrigOld[i] - temp) > disMin && abs(yOrigOld[i] - temp) > disMaxy){
            while (abs(yOrigOld[i] - temp) > disMaxy){
                j++;
                yn[j] = yn[j - 1] + disMaxy;
                temp = yn[j];
                sys->ny++;
            }
            if (abs(yOrigOld[i] - temp) > disMin){
                sys->ny++;
                temp = yOrigOld[i];
            }
            j++;
            yn[j] = yOrigOld[i];
        }
        else {
            j++;
            yn[j] = yOrigOld[i];
        }
    }

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
    /* Discretize domain in the z-direction */
    sort(zOrigOld, zOrigOld + 2 * sys->numStack + 2 * sys->numPorts);
    sys->nz = 1;
    double disMinz = minLayerDist * MINDISFRACZ; // Minimum discretization retained in z-direction after node merging is fraction of smallest distance between layers
    //double disMaxz = minLayerDist / MAXDISLAYERZ; // Maximum discretization distance in z-direction is fraction of closest distance between layers
    double *zn = (double*)calloc(2 * sys->numStack + 6 * sys->numPorts, sizeof(double));
    zn[0] = zOrigOld[0];
    j = 0;
    for (i = 1; i < 2 * sys->numStack + 2 * sys->numPorts; i++){
        if (abs(zOrigOld[i] - zOrigOld[i - 1]) > disMinz){
            sys->nz++;
        }
        j++;
        zn[j] = zOrigOld[i];
    }
    int countz = 2 * sys->numStack + 2 * sys->numPorts - 1;


    /*************************************************************************************/

    /* More discretization math */
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
    
    // Putting the layer permittivities in order of increasing z-coordinate?
    i = 0;
    if (sys->stackBegCoor[0] == sys->zn[0]){
        j = 0;
        sys->stackEpsn.push_back(sys->stackEps[j]);    // the stack eps with i < sys->N_edge_s
        while (i < sys->nz - 1) {
            if ((sys->zn[i] + sys->zn[i + 1]) / 2 >= sys->stackBegCoor[j] && (sys->zn[i] + sys->zn[i + 1]) / 2 <= sys->stackEndCoor[j]) {
                sys->stackEpsn.push_back(sys->stackEps[j]);
                i++;
            }
            else {
                j++;
            }
        }
    }
    else {
        j = sys->numStack - 1;
        sys->stackEpsn.push_back(sys->stackEps[j]);
        while (i < sys->nz - 1) {
            if ((sys->zn[i] + sys->zn[i + 1]) / 2 >= sys->stackBegCoor[j] && (sys->zn[i] + sys->zn[i + 1]) / 2 <= sys->stackEndCoor[j]) {
                sys->stackEpsn.push_back(sys->stackEps[j]);
                i++;
            }
            else {
                j--;
            }
        }
    }
    
    free(xOrigOld); xOrigOld = NULL;
    free(yOrigOld); yOrigOld = NULL;
    free(zOrigOld); zOrigOld = NULL;
    sys->stackEndCoor.clear();
    sys->stackBegCoor.clear();

#ifdef PRINT_NODE_COORD
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
#endif
    
    /***********************************************************************************************/

    /* Save counts of the final discretization */
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

    sys->markEdge = (int*)calloc(sys->N_edge, sizeof(int));   // mark which conductor index a given edge is inside
    //cout << "N_edge = " << sys->N_edge << ", sizeof(myint) = " << sizeof(myint) << ", their product is " << (size_t)(sys->N_edge * sizeof(myint)) << ", and sys->markEdge has size " << sizeof(sys->markEdge) << "/" << sizeof(sys->markEdge[0]) << endl;
    sys->markNode = (int*)calloc(sys->N_node, sizeof(int));   // mark which conductor index a given node is inside
    //cout << "N_node = " << sys->N_node << ", sizeof(myint) = " << sizeof(myint) << ", their product is " << sys->N_node * sizeof(myint) << ", and sys->markNode has size " << sizeof(sys->markNode) << "/" << sizeof(sys->markNode[0]) << endl;

#ifdef PRINT_VERBOSE_TIMING
    cout << "The time to read and assign x, y, z coordinates is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif

#ifdef PRINT_DIS_COUNT
    cout << endl;
    cout << "Discretization information: " << endl;
    cout << " disMin   = " << disMin << " m" << endl;
    cout << " disMinz  = " << disMinz << " m" << endl;
    cout << " disMaxx  = " << disMaxx << " m" << endl;
    cout << " disMaxy  = " << disMaxy << " m" << endl;
    //cout << " disMaxz  = " << disMaxz << " m" << endl;
    cout << endl;
    cout << " N_edge   = " << sys->N_edge << endl;
    cout << " N_node   = " << sys->N_node << endl;
    cout << " N_cell_x = " << sys->N_cell_x << endl;
    cout << " N_cell_y = " << sys->N_cell_y << endl;
    cout << " N_cell_z = " << sys->N_cell_z << endl;
    cout << endl;
#endif


    /* Discretization warnings */
    if (sys->N_cell_x <= 1)
    {
        cerr << "Failed to generate mesh with more than one element in the x-direction. Check value of disMin. Aborting now." << endl;
        return 1;
    }
    if (sys->N_cell_y <= 1)
    {
        cerr << "Failed to generate mesh with more than one element in the y-direction. Check value of disMin. Aborting now." << endl;
        return 1;
    }


    double xc, yc;
    myint xrange_max;
    unordered_map<myint, myint> xrange;
    vector<myint> xcoorv;
    set<myint> xcoor;
    unordered_map<myint, set<myint>> xcoory;    // store the start coordinate of the range, the end can be checked from xrange
    myint ss, ee;
    myint y1, y2, y3;
    myint ys, yl;
    myint x1, x2;
    myint l;
    tt = clock();
    int mark1, mini_k;
    double mini;
    // Fast algorithm to find nodes inside conductors
    
    for (i = 0; i < sys->numCdtRow; i++){
        xrange.clear();
        xcoor.clear();
        xcoorv.clear();
        xcoory.clear();
        mark1 = 0;
        //cout << "Number of CdtRow is " << i << endl;
        for (j = 0; j < sys->conductorIn[i].numVert - 1; j++){
            if (sys->conductorIn[i].x[j] != sys->conductorIn[i].x[j + 1] && sys->conductorIn[i].y[j] != sys->conductorIn[i].y[j + 1]){   // line with a non-zero slope
                if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[j + 1]){
                    x1 = xi[sys->conductorIn[i].x[j]];   // smaller x
                    x2 = xi[sys->conductorIn[i].x[j + 1]];    // larger x
                }
                else{
                    x1 = xi[sys->conductorIn[i].x[j + 1]];
                    x2 = xi[sys->conductorIn[i].x[j]];
                }
                for (l = x1; l <= x2; l++){
                    xcoor.insert(l);
                }
            }
            else{
                xcoor.insert(xi[sys->conductorIn[i].x[j]]);
            }
        }
        if (sys->conductorIn[i].x[j] != sys->conductorIn[i].x[0] && sys->conductorIn[i].y[j] != sys->conductorIn[i].y[0]){   // line with a non-zero slope
            if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[0]){
                x1 = xi[sys->conductorIn[i].x[j]];   // smaller x
                x2 = xi[sys->conductorIn[i].x[0]];    // larger x
            }
            else{
                x1 = xi[sys->conductorIn[i].x[0]];
                x2 = xi[sys->conductorIn[i].x[j]];
            }
            for (l = x1; l <= x2; l++){
                xcoor.insert(l);
            }
        }
        else{
            xcoor.insert(xi[sys->conductorIn[i].x[j]]);
        }

        for (auto xcoori : xcoor){
            xcoorv.push_back(xcoori);
        }
        xrange_max = xcoorv.back();   // the maximal x coordinate
        
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
                        if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){   // if this range is the rightmost one, or xrange only has one value, break
                            
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
            else if (sys->conductorIn[i].y[j] != sys->conductorIn[i].y[j + 1] && sys->conductorIn[i].x[j] != sys->conductorIn[i].x[j + 1]){   // this edge is with a slope
                if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[j + 1]){
                    ss = xi[sys->conductorIn[i].x[j]];
                    ee = xi[sys->conductorIn[i].x[j + 1]];
                    if (ss == ee && mark1 == 1){
                        continue;
                    }
                    while (xrange[ss] <= ee){
                        y3 = (sys->conductorIn[i].x[j + 1] - sys->xn[ss]) / (sys->conductorIn[i].x[j + 1] - sys->conductorIn[i].x[j]) * sys->conductorIn[i].y[j] + (sys->xn[ss] - sys->conductorIn[i].x[j]) / (sys->conductorIn[i].x[j + 1] - sys->conductorIn[i].x[j]) * sys->conductorIn[i].y[j + 1];
                        if (sys->conductorIn[i].y[j] > sys->conductorIn[i].y[j + 1]){
                            y1 = yi[sys->conductorIn[i].y[j + 1]];
                            y2 = yi[sys->conductorIn[i].y[j]];
                        }
                        else{
                            y1 = yi[sys->conductorIn[i].y[j]];
                            y2 = yi[sys->conductorIn[i].y[j + 1]];
                        }
                        mini = DOUBLEMAX;
                        mini_k = y1;
                        for (l = y1; l <= y2; l++){    // find the closest y to y3
                            if (mini < abs(sys->yn[l] - y3)){
                                mini = abs(sys->yn[l] - y3);
                                mini_k = l;
                            }
                        }
                        xcoory[ss].insert(mini_k);
                        if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){   // if this range is the rightmost one, or xrange only has one value, break
                            
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
                        y3 = (sys->conductorIn[i].x[j] - sys->xn[ss]) / (sys->conductorIn[i].x[j] - sys->conductorIn[i].x[j + 1]) * sys->conductorIn[i].y[j + 1] + (sys->xn[ss] - sys->conductorIn[i].x[j + 1]) / (sys->conductorIn[i].x[j] - sys->conductorIn[i].x[j + 1]) * sys->conductorIn[i].y[j];
                        if (sys->conductorIn[i].y[j] > sys->conductorIn[i].y[j + 1]){
                            y1 = yi[sys->conductorIn[i].y[j + 1]];
                            y2 = yi[sys->conductorIn[i].y[j]];
                        }
                        else{
                            y1 = yi[sys->conductorIn[i].y[j]];
                            y2 = yi[sys->conductorIn[i].y[j + 1]];
                        }
                        mini = DOUBLEMAX;
                        mini_k = y1;
                        for (l = y1; l <= y2; l++){    // find the closest y to y3
                            if (mini < abs(sys->yn[l] - y3)){
                                mini = abs(sys->yn[l] - y3);
                                mini_k = l;
                            }
                        }
                        xcoory[ss].insert(mini_k);
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
        else if (sys->conductorIn[i].y[j] != sys->conductorIn[i].y[0] && sys->conductorIn[i].x[j] != sys->conductorIn[i].x[0]){   // this edge is with a slope
            if (sys->conductorIn[i].x[j] < sys->conductorIn[i].x[0]){
                ss = xi[sys->conductorIn[i].x[j]];
                ee = xi[sys->conductorIn[i].x[0]];
                if (ss == ee && mark1 == 1){
                    continue;
                }
                while (xrange[ss] <= ee){
                    y3 = (sys->conductorIn[i].x[0] - sys->xn[ss]) / (sys->conductorIn[i].x[0] - sys->conductorIn[i].x[j]) * sys->conductorIn[i].y[j] + (sys->xn[ss] - sys->conductorIn[i].x[j]) / (sys->conductorIn[i].x[0] - sys->conductorIn[i].x[j]) * sys->conductorIn[i].y[0];
                    if (sys->conductorIn[i].y[j] > sys->conductorIn[i].y[0]){
                        y1 = yi[sys->conductorIn[i].y[0]];
                        y2 = yi[sys->conductorIn[i].y[j]];
                    }
                    else{
                        y1 = yi[sys->conductorIn[i].y[j]];
                        y2 = yi[sys->conductorIn[i].y[0]];
                    }
                    mini = DOUBLEMAX;
                    mini_k = y1;
                    for (l = y1; l <= y2; l++){    // find the closest y to y3
                        if (mini < abs(sys->yn[l] - y3)){
                            mini = abs(sys->yn[l] - y3);
                            mini_k = l;
                        }
                    }
                    xcoory[ss].insert(mini_k);
                    if (xrange.find(xrange[ss]) == xrange.end() || mark1 == 0){   // if this range is the rightmost one, or xrange only has one value, break
                        
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
                    y3 = (sys->conductorIn[i].x[j] - sys->xn[ss]) / (sys->conductorIn[i].x[j] - sys->conductorIn[i].x[0]) * sys->conductorIn[i].y[0] + (sys->xn[ss] - sys->conductorIn[i].x[0]) / (sys->conductorIn[i].x[j] - sys->conductorIn[i].x[0]) * sys->conductorIn[i].y[j];
                    if (sys->conductorIn[i].y[j] > sys->conductorIn[i].y[0]){
                        y1 = yi[sys->conductorIn[i].y[0]];
                        y2 = yi[sys->conductorIn[i].y[j]];
                    }
                    else{
                        y1 = yi[sys->conductorIn[i].y[j]];
                        y2 = yi[sys->conductorIn[i].y[0]];
                    }
                    mini = DOUBLEMAX;
                    mini_k = y1;
                    for (l = y1; l <= y2; l++){    // find the closest y to y3
                        if (mini < abs(sys->yn[l] - y3)){
                            mini = abs(sys->yn[l] - y3);
                            mini_k = l;
                        }
                    }
                    xcoory[ss].insert(mini_k);
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
        
        if (mark1 == 0){    // only has one x
            y1 = sys->ny;
            y2 = 0;
            for (auto xcooryi : xcoory){    // only one xcooryi
                for (auto xrangey : xcooryi.second){
                    if (y1 > xrangey){
                        y1 = xrangey;
                    }
                    if (y2 < xrangey){
                        y2 = xrangey;
                    }
                }
                j = xcooryi.first;
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
            continue;
        }
        for (auto xcooryi : xcoory){
            mark = 0;
            x1 = xcooryi.first;
            x2 = xrange[xcooryi.first];
            for (j = x1; j < x2; j++){
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
                                if (x2 == xrange_max && j == x2 - 1){    // If this range is the rightmost range, include the rightmost point
                                    sys->markNode[l * sys->N_node_s + x2 * (sys->N_cell_y + 1) + k] = 1;
                                    sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + x2 * (sys->N_cell_y) + k] = i + 1;
                                    if (l != zi[sys->conductorIn[i].zmax]){
                                        sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + x2 * (sys->N_cell_y + 1) + k] = i + 1;
                                    }
                                }
                            }
                            sys->markNode[l * sys->N_node_s + j * (sys->N_cell_y + 1) + k] = 1;
                            sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + j * (sys->N_cell_y + 1) + k] = i + 1;
                            if (l != zi[sys->conductorIn[i].zmax]){
                                sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j * (sys->N_cell_y + 1) + k] = i + 1;
                            }
                            if (x2 == xrange_max && j == x2 - 1){    // If this range is the rightmost range, include the rightmost point
                                sys->markNode[l * sys->N_node_s + x2 * (sys->N_cell_y + 1) + k] = 1;
                                if (l != zi[sys->conductorIn[i].zmax]){
                                    sys->markEdge[l * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + x2 * (sys->N_cell_y + 1) + k] = i + 1;
                                }
                            }
                        }

                    }
                }
            }
            
        }
        
    }

#ifdef PRINT_VERBOSE_TIMING
    cout << "The time to mark Edge and Node is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif

    // Assign markers to nodes and edges beyond included conductors
    for (i = 0; i < sys->N_edge_s; i++){    // the lower PEC plane
        sys->markEdge[i] = sys->numCdtRow + 1;
    }
    for (i = 0; i < sys->N_node_s; i++){
        sys->markNode[i] = 1;
    }
    for (i = sys->N_edge - sys->N_edge_s; i < sys->N_edge; i++){    // the upper PEC plane
        sys->markEdge[i] = sys->numCdtRow + 2;
    }
    for (i = sys->N_node - sys->N_node_s; i < sys->N_node; i++){
        sys->markNode[i] = 1;
    }

    /* experimental: store PEC planes as CdtRow to avoid segfault*/
    /* fdtdOneCondct class lacks a parametrized constructor or GDSII layer numbers to make this easy */
    /* do NOT update sys->numCdtRow afterwards */
    double xOuter[4] = { sys->xlim2, sys->xlim2, sys->xlim1, sys->xlim1 }; // Lower-right around counter-clockwise
    double yOuter[4] = { sys->ylim1, sys->ylim2, sys->ylim2, sys->ylim1 }; // Lower-right around counter-clockwise
    int layerMin = 65536;
    int layerMax = -1;
    for (size_t indi = 0; indi < sys->conductorIn.size(); indi++)
    {
        if (sys->conductorIn[indi].layer < layerMin)
        {
            layerMin = sys->conductorIn[indi].layer;
        }
        else if (sys->conductorIn[indi].layer > layerMax)
        {
            layerMax = sys->conductorIn[indi].layer;
        }
    }
    sys->conductorIn.push_back(fdtdOneCondct()); // the lower PEC plane
    sys->conductorIn.back().numVert = 4;
    sys->conductorIn.back().xmax = sys->xlim2;
    sys->conductorIn.back().xmin = sys->xlim1;
    sys->conductorIn.back().ymax = sys->ylim2;
    sys->conductorIn.back().ymin = sys->ylim1;
    sys->conductorIn.back().x = xOuter;
    sys->conductorIn.back().y = yOuter;
    sys->conductorIn.back().layer = layerMin;
    sys->conductorIn.push_back(fdtdOneCondct()); // the upper PEC plane
    sys->conductorIn.back().numVert = 4;
    sys->conductorIn.back().xmax = sys->xlim2;
    sys->conductorIn.back().xmin = sys->xlim1;
    sys->conductorIn.back().ymax = sys->ylim2;
    sys->conductorIn.back().ymin = sys->ylim1;
    sys->conductorIn.back().x = xOuter;
    sys->conductorIn.back().y = yOuter;
    sys->conductorIn.back().layer = layerMax;

    /* construct edgelink, no need for edgelink */
    myint eno;
    //sys->edgelink = (myint*)malloc(2 * sizeof(myint)*sys->N_edge);
    //for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++) {
    //    for (i = 1; i <= sys->N_cell_x + 1; i++) {    //edge along y axis
    //        for (j = 1; j <= sys->N_cell_y; j++) {
    //            eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (i - 1)*sys->N_cell_y + j;
    //            sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
    //            sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;

    //        }
    //    }
    //    for (i = 1; i <= sys->N_cell_x; i++) {    //edge along x axis
    //        for (j = 1; j <= sys->N_cell_y + 1; j++) {
    //            eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1)*sys->N_cell_y + (i - 1)*(sys->N_cell_y + 1) + j;
    //            sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
    //            sys->edgelink[(eno - 1) * 2 + 2 - 1] = (lyr - 1)*sys->N_node_s + i*(sys->N_cell_y + 1) + j - 1;

    //        }
    //    }
    //}
    //for (lyr = 1; lyr <= sys->N_cell_z; lyr++) {    // edge along z axis
    //    for (i = 1; i <= sys->N_cell_x + 1; i++) {
    //        for (j = 1; j <= sys->N_cell_y + 1; j++) {
    //            eno = (lyr - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (i - 1)*(sys->N_cell_y + 1) + j;
    //            sys->edgelink[(eno - 1) * 2 + 1 - 1] = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
    //            sys->edgelink[(eno - 1) * 2 + 2 - 1] = lyr*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j - 1;
    //        }
    //    }
    //}

    /* construct nodepos, no need for nodepos */
    myint nno;

    //sys->nodepos = (double*)malloc(sizeof(double)*sys->N_node * 3);   //N_node rows and 3 columns, input row by row
    //for (lyr = 1; lyr <= sys->N_cell_z + 1; lyr++) {
    //    for (i = 1; i <= sys->N_cell_x + 1; i++) {
    //        for (j = 1; j <= sys->N_cell_y + 1; j++) {
    //            nno = (lyr - 1)*sys->N_node_s + (i - 1)*(sys->N_cell_y + 1) + j;
    //            sys->nodepos[(nno - 1) * 3 + 1 - 1] = sys->xn[i - 1];
    //            sys->nodepos[(nno - 1) * 3 + 2 - 1] = sys->yn[j - 1];
    //            sys->nodepos[(nno - 1) * 3 + 3 - 1] = sys->zn[lyr - 1];
    //        }
    //    }
    //}

    /* construct nodeEdge */
    double leng;
    myint inz, inx, iny;
    myint node1, node2;

    /* no need for nodeEdge and nodeEdgea */
    //vector<pair<myint, double> > a;
    //for (i = 0; i < sys->N_node; i++) {
    //    sys->nodeEdge.push_back(a);
    //    sys->nodeEdgea.push_back(a);
    //}
    //for (i = 0; i < sys->N_edge; i++) {
    //    status = compute_edgelink(sys, i, node1, node2);
    //    if (i % (sys->N_edge_s + sys->N_edge_v) >= sys->N_edge_s){    // this edge is along z axis
    //        inz = i / (sys->N_edge_s + sys->N_edge_v);
    //        leng = sys->zn[inz + 1] - sys->zn[inz];
    //    }
    //    else if (i % (sys->N_edge_s + sys->N_edge_v) >= (sys->N_cell_y) * (sys->N_cell_x + 1)){    // this edge is along x axis
    //        inx = ((i % (sys->N_edge_s + sys->N_edge_v)) - (sys->N_cell_y) * (sys->N_cell_x + 1)) / (sys->N_cell_y + 1);
    //        leng = sys->xn[inx + 1] - sys->xn[inx];
    //    }
    //    else{    // this edge is along y axis
    //        iny = (i % (sys->N_edge_s + sys->N_edge_v)) % sys->N_cell_y;
    //        leng = sys->yn[iny + 1] - sys->yn[iny];
    //    }

    //    /*leng = pow((sys->nodepos[sys->edgelink[i * 2] * 3] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3]), 2);
    //    leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 1] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 1]), 2);
    //    leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 2] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 2]), 2);
    //    leng = sqrt(leng);*/
    //    sys->nodeEdge[node1].push_back(make_pair(i, 1 / leng));
    //    sys->nodeEdge[node2].push_back(make_pair(i, -1 / leng));
    //}
    //cout << "nodeEdge is done" << endl;
    //int ix, iy, iz;
    //iz = sys->N_cell_z;
    //for (ix = 0; ix < sys->nx; ix++) {
    //    for (iy = 0; iy < sys->ny; iy++) {
    //        sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, -1 / (sys->zn[iz] - sys->zn[iz - 1])));
    //    }
    //}
    //iz = 0;
    //for (ix = 0; ix < sys->nx; ix++) {
    //    for (iy = 0; iy < sys->ny; iy++) {
    //        sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->zn[iz + 1] - sys->zn[iz])));
    //    }
    //}
    //for (iz = 1; iz < sys->N_cell_z; iz++) {
    //    for (ix = 0; ix < sys->nx; ix++) {
    //        for (iy = 0; iy < sys->ny; iy++) {
    //            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, -2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
    //            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
    //        }
    //    }
    //}
    //ix = sys->N_cell_x;
    //for (iz = 0; iz < sys->nz; iz++) {
    //    for (iy = 0; iy < sys->ny; iy++) {
    //        sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -1 / (sys->xn[ix] - sys->xn[ix - 1])));
    //    }
    //}
    //ix = 0;
    //for (iz = 0; iz < sys->nz; iz++) {
    //    for (iy = 0; iy < sys->ny; iy++) {
    //        sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->xn[ix + 1] - sys->xn[ix])));
    //    }
    //}
    //for (ix = 1; ix < sys->N_cell_x; ix++) {
    //    for (iz = 0; iz < sys->nz; iz++) {
    //        for (iy = 0; iy < sys->ny; iy++) {
    //            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
    //            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
    //        }
    //    }
    //}
    //iy = sys->N_cell_y;
    //for (iz = 0; iz < sys->nz; iz++) {
    //    for (ix = 0; ix < sys->nx; ix++) {
    //        sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -1 / (sys->yn[iy] - sys->yn[iy - 1])));
    //    }
    //}
    //iy = 0;
    //for (iz = 0; iz < sys->nz; iz++) {
    //    for (ix = 0; ix < sys->nx; ix++) {
    //        sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 1 / (sys->yn[iy + 1] - sys->yn[iy])));
    //    }
    //}
    //for (iy = 1; iy < sys->N_cell_y; iy++) {
    //    for (iz = 0; iz < sys->nz; iz++) {
    //        for (ix = 0; ix < sys->nx; ix++) {
    //            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
    //            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
    //        }
    //    }
    //}

    /* implement bfs */
    myint *visited;
    vector<int> st;
    unordered_set<int> base;
    visited = (myint*)calloc(sys->N_node, sizeof(myint));
    count = 0;
    queue<myint> qu;


/* bfs to figure out the disjoint conductors */

tt = clock();
for (i = 0; i < sys->N_node; i++) {
    if (sys->markNode[i] == 0) {
        continue;
    }
    else {
        if (visited[i] != 0) {
            continue;
        }
        else {
            
            qu.push(i);
            count++;
            visited[i] = count;
            //sys->cond2condIn.push_back(base);
            while (!qu.empty()) {
                mark = 0;
                inz = qu.front() / (sys->N_node_s);
                inx = (qu.front() % sys->N_node_s) / (sys->N_cell_y + 1);
                iny = (qu.front() % sys->N_node_s) % (sys->N_cell_y + 1);
                
                // how many edges this node connects to
                if (inz != 0){    // this node is not on the bottom plane
                    eno = (inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny;    // the node's lower edge
                    if (sys->markEdge[eno] != 0){
                        
                        status = compute_edgelink(sys, eno, node1, node2);    // compute_edgelink is used to get edge eno's two side's nodes node1, node2
                        if ((node1 != qu.front() && visited[node1] == 0)) {
                            visited[node1] = count;
                            qu.push(node1);
                        }
                        else if ((node2 != qu.front() && visited[node2] == 0)) {
                            visited[node2] = count;
                            qu.push(node2);
                        }
                    }
                }
                if (inz != sys->nz - 1){    // this node is not on the upper plane
                    eno = inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny;    // the node's lower edge
                    if (sys->markEdge[eno] != 0) {
                        
                        status = compute_edgelink(sys, eno, node1, node2);
                        if ((node1 != qu.front() && visited[node1] == 0)) {
                            visited[node1] = count;
                            qu.push(node1);
                        }
                        else if ((node2 != qu.front() && visited[node2] == 0)) {
                            visited[node2] = count;
                            qu.push(node2);
                        }

                    }
                }
                if (inx != 0){    // this node is not on the left plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (inx - 1) * (sys->N_cell_y + 1) + iny;
                    if (sys->markEdge[eno] != 0) {
                        status = compute_edgelink(sys, eno, node1, node2);
                        if ((node1 != qu.front() && visited[node1] == 0)) {
                            visited[node1] = count;
                            qu.push(node1);
                        }
                        else if ((node2 != qu.front() && visited[node2] == 0)) {
                            visited[node2] = count;
                            qu.push(node2);
                        }

                    }
                }
                if (inx != sys->nx - 1){    // this node is not on the right plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + inx * (sys->N_cell_y + 1) + iny;
                    if (sys->markEdge[eno] != 0) {

                        status = compute_edgelink(sys, eno, node1, node2);
                        if ((node1 != qu.front() && visited[node1] == 0)) {
                            visited[node1] = count;
                            qu.push(node1);
                        }
                        else if ((node2 != qu.front() && visited[node2] == 0)) {
                            visited[node2] = count;
                            qu.push(node2);
                        }

                    }
                }
                if (iny != 0){    // this node is not on the front plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny - 1;
                    if (sys->markEdge[eno] != 0) {
                        
                        status = compute_edgelink(sys, eno, node1, node2);
                        if ((node1 != qu.front() && visited[node1] == 0)) {
                            visited[node1] = count;
                            qu.push(node1);
                        }
                        else if ((node2 != qu.front() && visited[node2] == 0)) {
                            visited[node2] = count;
                            qu.push(node2);
                        }

                    }
                }
                if (iny != sys->ny - 1){    // this node is not on the back plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny;
                    if (sys->markEdge[eno] != 0) {
                        
                        status = compute_edgelink(sys, eno, node1, node2);
                        if ((node1 != qu.front() && visited[node1] == 0)) {
                            visited[node1] = count;
                            qu.push(node1);
                        }
                        else if ((node2 != qu.front() && visited[node2] == 0)) {
                            visited[node2] = count;
                            qu.push(node2);
                        }

                    }
                }
                qu.pop();
            }
        }
    }
}

#ifdef PRINT_VERBOSE_TIMING
    cout << "Time to find conductors is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
    cout << "Number of conductors is " << count << endl;
#endif

    /*for (i = 0; i < sys->N_node; i++){
    if (visited[i] != 0){
    for (j = 0; j < sys->nodeEdge[i].size(); j++){
    if (sys->markEdge[sys->nodeEdge[i][j].first] != 0 && sys->cond2condIn[visited[i] - 1].find(sys->markEdge[sys->nodeEdge[i][j].first]) == sys->cond2condIn[visited[i] - 1].end()){
    sys->cond2condIn[visited[i] - 1].insert(sys->markEdge[sys->nodeEdge[i][j].first]);
    }
    }
    }
    }*/

    tt = clock();
    for (i = 0; i < sys->N_node; i++){
        sys->markNode[i] = visited[i];
        
    }

#ifdef PRINT_VERBOSE_TIMING
    cout << "Time to assign markEdge and markNode is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif

    tt = clock();

    /* Construct each isolated conductor */
    sys->numCdt = count;
    /*for (i = 0; i < sys->numCdt; i++) {
    cout << "Cnt " << i << " has condIn as: ";
    for (auto ci : sys->cond2condIn[i]) {
    cout << ci << " ";
    }
    cout << endl;
    }*/
    tt = clock();
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
                sys->conductor[visited[i] - 1].markPort = -1;    // this conductor connects to the lower PEC
            }
            else if ((i >= sys->N_node - sys->N_node_s) && sys->conductor[visited[i] - 1].markPort != -2 && sys->conductor[visited[i] - 1].markPort != -1){
                sys->conductor[visited[i] - 1].markPort = -2;    // this conductor connects to the upper PEC
            }
            sys->conductor[visited[i] - 1].cdtNodeind++;
        }
    }
    
    myint indPortNode1, indPortNode2;
    set<int> cond;    // the port conductors
#ifdef PRINT_VERBOSE_TIMING
    cout << "Time to assign conductor parameter is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif
    cout << "Number of conductors is " << sys->numCdt << endl;

    tt = clock();
    for (i = 0; i < sys->numPorts; i++){
        indPortNode1 = sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]];
        indPortNode2 = sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]];
        
        if (cond.find(indPortNode1) == cond.end()){
            cond.insert(indPortNode1);
            for (j = 0; j < sys->cdtNumNode[indPortNode1 - 1]; j++){
                inz = sys->conductor[indPortNode1 - 1].node[j] / sys->N_node_s;
                inx = ((sys->conductor[indPortNode1 - 1].node[j]) % sys->N_node_s) / (sys->N_cell_y + 1);
                iny = ((sys->conductor[indPortNode1 - 1].node[j]) % sys->N_node_s) % (sys->N_cell_y + 1);
                if (inz != 0){    // this node is not on the bottom plane
                    eno = (inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny;    // the node's lower edge
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()){
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
                if (inz != sys->nz - 1){    // this node is not on the upper plane
                    eno = inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny;    // the node's lower edge
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
                if (inx != 0){    // this node is not on the left plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (inx - 1) * (sys->N_cell_y + 1) + iny;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);

                    }
                }
                if (inx != sys->nx - 1){    // this node is not on the right plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + inx * (sys->N_cell_y + 1) + iny;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
                if (iny != 0){    // this node is not on the front plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny - 1;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);

                    }
                }
                if (iny != sys->ny - 1){    // this node is not on the back plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
            }
        }
        
        if (cond.find(indPortNode2) == cond.end()){
            cond.insert(indPortNode2);
            for (j = 0; j < sys->cdtNumNode[indPortNode2 - 1]; j++){
                inz = sys->conductor[indPortNode2 - 1].node[j] / sys->N_node_s;
                inx = ((sys->conductor[indPortNode2 - 1].node[j]) % sys->N_node_s) / (sys->N_cell_y + 1);
                iny = ((sys->conductor[indPortNode2 - 1].node[j]) % sys->N_node_s) % (sys->N_cell_y + 1);
                if (inz != 0){    // this node is not on the bottom plane
                    eno = (inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny;    // the node's lower edge
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()){
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
                if (inz != sys->nz - 1){    // this node is not on the upper plane
                    eno = inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny;    // the node's lower edge
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
                if (inx != 0){    // this node is not on the left plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (inx - 1) * (sys->N_cell_y + 1) + iny;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);

                    }
                }
                if (inx != sys->nx - 1){    // this node is not on the right plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + inx * (sys->N_cell_y + 1) + iny;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
                if (iny != 0){    // this node is not on the front plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny - 1;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);

                    }
                }
                if (iny != sys->ny - 1){    // this node is not on the back plane
                    eno = inz *(sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny;
                    if (sys->markEdge[eno] != 0 && sys->cond2condIn.find(sys->markEdge[eno]) == sys->cond2condIn.end()) {
                        sys->cond2condIn.insert(sys->markEdge[eno]);
                    }
                }
            }
        }
    }
    cond.clear();
    cout << "cond2condIn is set sucessfully!" << endl;
#ifdef PRINT_VERBOSE_TIMING
    cout << "Time to construct sys->cond2condIn (size " <<  sys->cond2condIn.size() << ") is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif

    sys->outedge = 0;
    sys->inedge = 0;
    tt = clock();
    for (i = 0; i < sys->N_edge; i++){
        // find edge i's two nodes
        status = compute_edgelink(sys, i, node1, node2);

        if (sys->markEdge[i] != 0 && visited[node1] == visited[node2] && visited[node1] != 0){
            sys->markEdge[i] = visited[node2];    // Mark the edge with each color for different conductors
            
        }
        if (sys->markEdge[i] != 0){
            sys->inedge++;
        }
        else{
            sys->outedge++;
        }
    }

#ifdef PRINT_VERBOSE_TIMING
    cout << "Time to assign markEdge is " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
    cout << "The number of inside conductor edge is " << sys->inedge << endl;
    cout << "The number of outside conductor edge is " << sys->outedge << endl;
#endif

    free(visited);
    visited = NULL;

    /* set markCell */
#ifndef SKIP_MARK_CELL
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
#endif

    return 0;
}

int compute_edgelink(fdtdMesh *sys, myint eno, myint &node1, myint &node2){
    int inz, inx, iny;
    if (eno % (sys->N_edge_s + sys->N_edge_v) >= sys->N_edge_s){    // this edge is along z axis
        inz = eno / (sys->N_edge_s + sys->N_edge_v);
        inx = ((eno % (sys->N_edge_s + sys->N_edge_v)) - sys->N_edge_s) / (sys->N_cell_y + 1);
        iny = ((eno % (sys->N_edge_s + sys->N_edge_v)) - sys->N_edge_s) % (sys->N_cell_y + 1);
        node1 = inz * sys->N_node_s + (sys->N_cell_y + 1) * inx + iny;
        node2 = (inz + 1) * sys->N_node_s + (sys->N_cell_y + 1) * inx + iny;
    }
    else if (eno % (sys->N_edge_s + sys->N_edge_v) >= (sys->N_cell_y) * (sys->N_cell_x + 1)){    // this edge is along x axis
        inz = eno / (sys->N_edge_s + sys->N_edge_v);
        inx = ((eno % (sys->N_edge_s + sys->N_edge_v)) - (sys->N_cell_y) * (sys->N_cell_x + 1)) / (sys->N_cell_y + 1);
        iny = ((eno % (sys->N_edge_s + sys->N_edge_v)) - (sys->N_cell_y) * (sys->N_cell_x + 1)) % (sys->N_cell_y + 1);
        node1 = inz * sys->N_node_s + inx * (sys->N_cell_y + 1) + iny;
        node2 = inz * sys->N_node_s + (inx + 1) * (sys->N_cell_y + 1) + iny;
    }
    else{    // this edge is along y axis
        inz = eno / (sys->N_edge_s + sys->N_edge_v);
        inx = (eno % (sys->N_edge_s + sys->N_edge_v)) / sys->N_cell_y;
        iny = (eno % (sys->N_edge_s + sys->N_edge_v)) % sys->N_cell_y;
        node1 = inz * sys->N_node_s + inx * (sys->N_cell_y + 1) + iny;
        node2 = inz * sys->N_node_s + inx * (sys->N_cell_y + 1) + iny + 1;
    }
    return 0;
}

int matrixConstruction(fdtdMesh *sys){
    int i, j;
    //cout << sys->stackEpsn.size() << endl;
    /* construct D_eps */
    /*sys->eps = (double*)malloc(sizeof(double)*sys->N_edge);
    for (i = 0; i < sys->N_edge; i++){
        if (i < sys->N_edge_s){
            sys->eps[i] = sys->stackEpsn[0] * EPSILON0;
            continue;
        }
        else{
            sys->eps[i] = sys->stackEpsn[(i - sys->N_edge_s) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
        }
    }*/
    
    /* construct D_sig */
    /*sys->sig = (double*)calloc(sys->N_edge, sizeof(double));
    double a, b;
    for (i = 0; i < sys->N_edge; i++){
        if (sys->markEdge[i] != 0){
#ifndef SKIP_MARK_CELL
            a = 0;
            b = 0;
            for (j = 0; j < sys->edgeCell[i].size(); j++){
                a += sys->markCell[sys->edgeCell[i][j]] * sys->edgeCellArea[i][j];
                b += sys->edgeCellArea[i][j];
            }
            sys->sig[i] = a / b * SIGMA;
#else
            sys->sig[i] = SIGMA;
#endif
        }
    }*/
    

    sys->edgeCell.clear();
    sys->edgeCellArea.clear();
    //free(sys->markCell); sys->markCell = NULL;


    /*sys->stackEpsn.clear();*/
    

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
        //cout << "Reminder to send error if sys->portCoor[i] == 0 because program will fail later" << endl;
        myint indMarkNode1 = sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]];
        myint indMarkNode2 = sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]];
        
        if ((indMarkNode1 != 0) && (sys->conductor[indMarkNode1 - 1].markPort > -1)){ // Only supply end of port not on PEC
            sys->portCoor[i].portCnd = indMarkNode1;
            sys->conductor[indMarkNode1 - 1].markPort = i + 1;    // markPort start from 1
        }
        else if ((indMarkNode2 != 0) && (sys->conductor[indMarkNode2 - 1].markPort > -1)){ // Only return end of port not on PEC
            sys->portCoor[i].portCnd = indMarkNode2;
            sys->conductor[indMarkNode2 - 1].markPort = i + 1;    // markPort start from 1
        }
        else if (indMarkNode1 != 0 && indMarkNode2 != 0) {
            if (sys->conductor[indMarkNode1 - 1].markPort == -2){    // the upper PEC conductor is excited
                sys->portCoor[i].portCnd = indMarkNode1;
            }
            else if (sys->conductor[indMarkNode2 - 1].markPort == -2){
                sys->portCoor[i].portCnd = indMarkNode2;
            }
            //sys->portCoor[i].portCnd = indMarkNode2;
            /*sys->conductor[indMarkNode1 - 1].markPort = i + 1;
            sys->conductor[indMarkNode2 - 1].markPort = i + 1;*/
        }
        /*for (j = 0; j < sys->cdtNumNode[sys->portCoor[i].portCnd - 1]; j++){
        sys->exciteCdtLayer[sys->conductor[sys->portCoor[i].portCnd - 1].node[j] / sys->N_node_s] = 1;
        }*/
#ifdef PRINT_PORT_SET
        cout << "portSet() logic test: markNode on Point 1 (" << indMarkNode1 << "), port marker on corresponding isolated conductor for Point 1 (" << sys->conductor[indMarkNode1 - 1].markPort << "), markNode on Point 2 (" << indMarkNode2 << "), port marker on corresponding isolated conductor for Point 2 (" << sys->conductor[indMarkNode2 - 1].markPort << ")" << endl;
        cout << "Value of portCnd for port #" << i << ": " << sys->portCoor[i].portCnd << endl;
#endif
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
                a *= (sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
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
                a *= (sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
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
                a *= (sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1);
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
    if (sideLen != 0){

        for (auto ci : sys->cond2condIn){
            //cout << "ci = " << ci << " out of " << sys->cond2condIn[sys->portCoor[i].portCnd - 1].size() << " possible in unordered set" << endl;
            for (l = 0; l < sys->conductorIn[ci - 1].numVert - 1; l++){
                //cout << "l = " << l << " out of " << sys->conductorIn[ci - 1].numVert - 1 << endl;
                if (sys->conductorIn[ci - 1].x[l] == sys->conductorIn[ci - 1].x[l + 1]){
                    //cout << "Made it to first if statement in auto ci" << endl;
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
                    //cout << "Made it to first else statement in auto ci" << endl;
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
                //cout << "Made it to second if statement in auto ci" << endl;
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
                //cout << "Made it to second else statement in auto ci" << endl;
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
                //cout << "Did all index nonsense within second else statement in auto ci" << endl;
                for (k = z1_ind; k <= z2_ind; k++){
                    for (j = x1_ind; j <= x2_ind; j++){
                        for (m = y1_ind; m <= y2_ind; m++){
                            /*if (k * sys->N_node_s + j * (sys->N_cell_y + 1) + m >= sys->N_node) {
                                cout << "Checking index " << k * sys->N_node_s + j * (sys->N_cell_y + 1) + m << " of size-" << sys->N_node << " array markNode against 0" << endl;
                                }*/
                            if (sys->markNode[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] == 0) {
                                sys->markProSide[k * sys->N_node_s + j * (sys->N_cell_y + 1) + m] = 1;
                            }
                        }
                    }
                }
                //cout << "Did all markProSide nonsense within second else statement in auto ci" << endl;
            }
        }
    }

#ifdef PRINT_VERBOSE_TIMING
    cout << "Time of finding side nodes is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif

    sys->conductorIn.clear();
    return 0;
}

/* Is point (x,y) within the polygon? */
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
    //cout << "  nodepos array exists (" << (this->nodepos != nullptr) << ")" << endl;
    cout << "  Epoints array exists (" << (this->Epoints != nullptr) << ")" << endl;
   //cout << "  edgelink array exists (" << (this->edgelink != nullptr) << ")" << endl;
    cout << "  Hpoints array exists (" << (this->Hpoints != nullptr) << ")" << endl;
    //cout << "  nodeEdge vector has size " << this->nodeEdge.size() << endl;
    cout << "  nodeEdgea vector has size " << this->nodeEdgea.size() << endl;
    cout << " PEC information:" << endl;
    cout << "  Boundary node 1 array exists (" << (this->bd_node1 != nullptr) << ")" << endl;
    cout << "  Boundary node 2 array exists (" << (this->bd_node2 != nullptr) << ")" << endl;
    cout << "  Boundary edge array exists (" << (this->bd_edge != nullptr) << ")" << endl;
    cout << " Layer stack up parameters:" << endl;
    cout << "  Number of layers in stack: " << this->numStack << endl;
    cout << "  Permittivity vector has size " << this->stackEps.size() << endl;
    cout << "  Beginning z-coordinate vector has size " << this->stackBegCoor.size() << endl;
    cout << "  Ending z-coordinate vector has size " << this->stackEndCoor.size() << endl;
    cout << "  Layer name vector has size " << this->stackName.size() << endl;
    if (this->eps == nullptr)
    {
        cout << "  Edge permittivity array D_eps exists (" << (this->eps != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Edge permittivity array D_eps has size " << NELEMENT(this->eps) << endl;
    }
    cout << "  stackEpsn vector has size " << this->stackEpsn.size() << endl;
    if (this->stackCdtMark == nullptr)
    {
        cout << "  Stack conductor marker array exists (" << (this->stackCdtMark != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Stack conductor marker array has size " << NELEMENT(this->stackCdtMark) << endl;
    }
    cout << " Conductor parameter information:" << endl;
    cout << "  conductorIn vector has size " << this->conductorIn.size() << endl;
    cout << "  Number of conductor rows: " << this->numCdtRow << endl;
    cout << "  Number of isolated conductors: " << this->numCdt << endl;
    if (this->markEdge == nullptr)
    {
        cout << "  Edge marker array exists (" << (this->markEdge != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Edge marker array has size " << NELEMENT(this->markEdge) << endl;
    }
    if (this->markCell == nullptr)
    {
        cout << "  Cell marker array exists (" << (this->markCell != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Cell marker array has size " << NELEMENT(this->markCell) << endl;
    }
    if (this->cdtNumNode == nullptr)
    {
        cout << "  cdtNumNode array exists (" << (this->cdtNumNode != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  cdtNumNode array has size " << NELEMENT(this->cdtNumNode) << endl;
    }
    if (this->sig == nullptr)
    {
        cout << "  Edge conductivity array D_sig exists (" << (this->sig != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Edge conductivity array D_sig has size " << NELEMENT(this->sig) << endl;
    }
    if (this->conductor == nullptr)
    {
        cout << "  conductor array exists (" << (this->conductor != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  conductor array has size " << NELEMENT(this->conductor) << endl;
    }
    if (this->markNode == nullptr)
    {
        cout << "  Node marker array exists (" << (this->markNode != nullptr) << ")" << endl;
    }
    else
    {
        cout << "  Node marker array has size " << NELEMENT(this->markNode) << endl;
    }
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
    cout << "  x vector has size " << this->x.size() << endl;
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
