//#include "stdafx.h"
#include "fdtd.h"

static bool comp(pair<double, int> a, pair<double, int> b)
{
    return a.second <= b.second;
};

int readInput(const char *stackFile, fdtdMesh *sys, unordered_map<double, int> &xi, unordered_map<double, int> &yi, unordered_map<double, int> &zi){
    FILE *fp, *cfp;
    char fbase[FDTD_MAXC];
    char s[FDTD_MAXC];
    char *word[FDTD_MAXC];
    int lyr;
    int i;
    int j, k, m;
    double upper, lower;
    int status;
    int count;
    int mark;
    double xmin, xmax, ymin, ymax, xwid, ywid;
    unordered_set<double> portCoorx, portCoory;


    xmin = DOUBLEMAX;    // whole structure xmin
    xmax = DOUBLEMIN;    // whole structure xmax
    ymin = DOUBLEMAX;    // whole structure ymin
    ymax = DOUBLEMIN;    // whole structure ymax

    /*sys->ix = 2;
    sys->iy = 4;*/

    /*strcpy(fbase, argv[1]);
    strcat(fbase, ".sim_input");*/   //cause segmentation fault


    /*read polygon file*/
    lyr = 0;

    fp = fopen(stackFile, "r");
    /*cfp = fopen(polyFile, "r");*/

    if (fp == NULL){
        perror("Failed: ");
        return -2;
    }
    /*if (cfp == NULL){
        perror("Failed: ");
        return -2;
        }*/

    cout << "Begin reading the stack file!" << endl;
    while (fgets(s, FDTD_MAXC, fp) != NULL){
        fdtdStringWord(s, word);
        if (!strncmp(word[0], "LENGTHUNIT", 10)){
            sys->lengthUnit = fdtdGetValue(word[2]);
            break;
        }
    }

    while (fgets(s, FDTD_MAXC, fp) != NULL){
        fdtdStringWord(s, word);
        if (!strncmp(word[0], "FREQUENCY", 9)){
            if (fgets(s, FDTD_MAXC, fp) != NULL){
                fdtdStringWord(s, word);
                sys->freqUnit = fdtdGetValue(word[2]);
            }
            if (fgets(s, FDTD_MAXC, fp) != NULL){
                fdtdStringWord(s, word);
                sys->freqStart = fdtdGetValue(word[2]);
            }
            if (fgets(s, FDTD_MAXC, fp) != NULL){
                fdtdStringWord(s, word);
                sys->freqEnd = fdtdGetValue(word[2]);
            }
            if (fgets(s, FDTD_MAXC, fp) != NULL){
                fdtdStringWord(s, word);
                sys->nfreq = fdtdGetValue(word[2]);
            }
            break;
        }
    }
    
    while (fgets(s, FDTD_MAXC, fp) != NULL){
        fdtdStringWord(s, word);
        if (!strncmp(word[0], "NUMSTACK", 8)){
            sys->numStack = (int)fdtdGetValue(word[2]);
            break;
        }
    }

    lyr = 0;
    int lyr1 = 0;
    sys->stackEps = (double*)calloc(sys->numStack, sizeof(double));
    sys->stackBegCoor = (double*)calloc(sys->numStack, sizeof(double));
    sys->stackEndCoor = (double*)calloc(sys->numStack, sizeof(double));
    sys->stackCdtMark = (double*)calloc(sys->numStack, sizeof(double));
    //cout << sys->numStack << endl;
    while (fgets(s, FDTD_MAXC, fp) != NULL){
        fdtdStringWord(s, word);
        if (lyr >= sys->numStack){    //blank
            break;
        }
        if (fdtdGetValue(word[3]) == 0){
            lyr++;
            continue;
        }
        sys->stackEps[lyr1] = fdtdGetValue(word[6]);   // get the stack eps
        if (lyr1 == 0){
            sys->stackBegCoor[lyr1] = 0;
        }
        else{
            sys->stackBegCoor[lyr1] = sys->stackEndCoor[lyr1 - 1];
        }
        sys->stackEndCoor[lyr1] = sys->stackBegCoor[lyr1] + fdtdGetValue(word[3])*sys->lengthUnit;
        sys->stackName.push_back(word[0]);
        lyr1++;
        lyr++;
    }
    /*sys->numStack = lyr1;*/
    /*for (i = 0; i < sys->numStack; i++){
        cout << sys->stackBegCoor[i] << " " << sys->stackEndCoor[i] << " " << sys->stackEps[i] << endl;
        }*/

    while (fgets(s, FDTD_MAXC, fp) != NULL){
        fdtdStringWord(s, word);
        if (!strncmp(word[0], "NUMPORT", 7)){
            sys->numPorts = (int)fdtdGetValue(word[2]);
            break;
        }
    }

    int temp;
    lyr = 0;
    sys->portCoor = (fdtdPort*)malloc(sizeof(fdtdPort) * sys->numPorts);
    while (fgets(s, FDTD_MAXC, fp) != NULL){
        fdtdStringWord(s, word);
        if (lyr >= sys->numPorts)
            break;
        sys->portCoor[lyr].x1 = fdtdGetValue(word[0])*sys->lengthUnit;
        sys->portCoor[lyr].y1 = fdtdGetValue(word[1])*sys->lengthUnit;
        sys->portCoor[lyr].x2 = fdtdGetValue(word[3])*sys->lengthUnit;
        sys->portCoor[lyr].y2 = fdtdGetValue(word[4])*sys->lengthUnit;
        sys->portCoor[lyr].z1 = fdtdGetValue(word[2])*sys->lengthUnit;
        sys->portCoor[lyr].z2 = fdtdGetValue(word[5])*sys->lengthUnit;
        
        if (sys->portCoor[lyr].x1 == sys->portCoor[lyr].x2){
            portCoorx.insert(sys->portCoor[lyr].x1);    // the x coordinates on the yz plane is recorded
        }
        if (sys->portCoor[lyr].y1 == sys->portCoor[lyr].y2){
            portCoory.insert(sys->portCoor[lyr].y1);    // the y coordinates on the xz plane is recorded
        }
        if (sys->portCoor[lyr].x1 > sys->portCoor[lyr].x2){
            temp = sys->portCoor[lyr].x1;
            sys->portCoor[lyr].x1 = sys->portCoor[lyr].x2;
            sys->portCoor[lyr].x2 = temp;
        }
        if (sys->portCoor[lyr].y1 > sys->portCoor[lyr].y2){
            temp = sys->portCoor[lyr].y1;
            sys->portCoor[lyr].y1 = sys->portCoor[lyr].y2;
            sys->portCoor[lyr].y2 = temp;
        }
        if (sys->portCoor[lyr].z1 > sys->portCoor[lyr].z2){
            temp = sys->portCoor[lyr].z1;
            sys->portCoor[lyr].z1 = sys->portCoor[lyr].z2;
            sys->portCoor[lyr].z2 = temp;
        }
        sys->portCoor[lyr].portDirection = (int)fdtdGetValue(word[6]);
        lyr++;
    }

    
    
    cout << "Begin reading the conductor information!" << endl;
    //sys->conductorIn = (fdtdOneCondct*)malloc(sizeof(fdtdOneCondct)*sys->numCdtRow);
    i = 0;
    //k = 0;
    while (true){    //get one line and put each word in word arr
        /*count = 0;
        while (strPoints[k] != '\n'){
            s[count] = strPoints[k];
            count++;
            k++;
        }
        s[count] = strPoints[k];
        /*for (m = 0; m <= count; m++)
            cout << s[m];*/
        /*k++;*/
        fdtdStringWord(s, word);
        if (i >= sys->numCdtRow)    //if on this line there is nothing
            break;
        else{
            /*sys->conductorIn[i].numVert = (int)fdtdGetValue(word[0]);
            sys->conductorIn[i].xmax = DOUBLEMIN;
            sys->conductorIn[i].xmin = DOUBLEMAX;
            sys->conductorIn[i].ymax = DOUBLEMIN;
            sys->conductorIn[i].ymin = DOUBLEMAX;
            sys->conductorIn[i].x = (double*)calloc(sys->conductorIn[i].numVert, sizeof(double));
            sys->conductorIn[i].y = (double*)calloc(sys->conductorIn[i].numVert, sizeof(double));
            for (j = 0; j < sys->conductorIn[i].numVert; j++){
                if (fdtdGetValue(word[j * 2 + 2]) > sys->conductorIn[i].xmax){
                    sys->conductorIn[i].xmax = fdtdGetValue(word[j * 2 + 2]);
                }
                if (fdtdGetValue(word[j * 2 + 2]) < sys->conductorIn[i].xmin){
                    sys->conductorIn[i].xmin = fdtdGetValue(word[j * 2 + 2]);
                }
                if (fdtdGetValue(word[j * 2 + 3]) > sys->conductorIn[i].ymax){
                    sys->conductorIn[i].ymax = fdtdGetValue(word[j * 2 + 3]);
                }
                if (fdtdGetValue(word[j * 2 + 3]) < sys->conductorIn[i].ymin){
                    sys->conductorIn[i].ymin = fdtdGetValue(word[j * 2 + 3]);
                }
                sys->conductorIn[i].x[j] = fdtdGetValue(word[j * 2 + 2]);
                sys->conductorIn[i].y[j] = fdtdGetValue(word[j * 2 + 3]);
            }*/
            //lyr = (int)fdtdGetValue(word[1]);

            for (j = 0; j < sys->numStack; j++){
                if (to_string(sys->conductorIn[i].layer) == sys->stackName[j]){
                    sys->conductorIn[i].zmin = sys->stackBegCoor[j];
                    sys->conductorIn[i].zmax = sys->stackEndCoor[j];
                    sys->stackCdtMark[j] = 1;
                }
            }
        }
        i++;
    }
    

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
    int xMaxInd = 10;
    disMaxx = (xmax - xmin) / xMaxInd;

    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(xOrigOld[i] - xOrigOld[i - 1]) > disMin){
            sys->nx++;
        }
    }
    double *xn = (double*)calloc(numNode + 6 * sys->numPorts + xMaxInd, sizeof(double));
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


    //sys->xnu = (double*)calloc(sys->nx, sizeof(double));

    j = 0;
    //sys->xnu[0] = xn[0];
    //xi[sys->xnu[0]] = j;
    double first, second;
    for (i = 1; i <= countx; i++){    // set the discretization length around port to be equal
        if (abs(xn[i] - xn[i - 1]) > disMin){
            j++;
            //sys->xnu[j] = xn[i];
            //xi[sys->xnu[j]] = j;
        }
        else{
            //xi[xn[i]] = j;
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
    int yMaxInd = 10;    // the max discretization of y is total / 120
    disMaxy = (ymax - ymin) / yMaxInd;

    for (i = 1; i < numNode + 2 * sys->numPorts; i++){
        if (abs(yOrigOld[i] - yOrigOld[i - 1]) > disMin){
            sys->ny++;
        }
    }
    double *yn = (double*)calloc(numNode + 6 * sys->numPorts + yMaxInd, sizeof(double));
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

    //sys->ynu = (double*)calloc(sys->ny, sizeof(double));

    
    j = 0;
    //sys->ynu[0] = yn[0];
    //yi[sys->ynu[0]] = j;
    
    for (i = 1; i <= county; i++){    // set the discretization length around port to be equal
        if (abs(yn[i] - yn[i - 1]) > disMin){
            j++;
            //sys->ynu[j] = yn[i];
            //yi[sys->ynu[j]] = j;
        }
        else{
            //yi[yn[i]] = j;
        }
    }
    
    
    sys->ny = j + 1;
    /*vector<pair<double, int> > v(yi.begin(), yi.end());
    sort(v.begin(), v.end(), comp);
    for (i = 0; i < v.size(); i++){
    cout << v[i].first << " " << v[i].second << endl;
    }*/
    //cout << sys->ny << endl;

    

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

    //sys->znu = (double*)calloc(sys->nz, sizeof(double));
    //sys->znu[0] = zOrigOld[0];
    //zi[sys->znu[0]] = 0;
    j = 0;
    for (i = 1; i < 2 * sys->numStack + 2 * sys->numPorts; i++){
        if (abs(zOrigOld[i] - zOrigOld[i - 1]) > disMin){
            j++;
            //sys->znu[j] = zOrigOld[i];
            //zi[sys->znu[j]] = j;
        }
        else{
            //zi[zOrigOld[i]] = j;    // considering the double precision, put all the show up y in the map
        }
    }
    

    /*************************************************************************************/

    
    sort(xn, xn + countx + 1);
    xi.clear();
    sys->xn = (double*)calloc(sys->nx + 4 * sys->numPorts, sizeof(double));
    j = 0;
    sys->xn[0] = xn[0];
    //xi[sys->xn[0]] = j;
    for (i = 1; i <= countx; i++){    // set the discretization length around port to be equal
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
    free(xn); xn = NULL;
    
    sort(yn, yn + county + 1);
    yi.clear();
    sys->yn = (double*)calloc(sys->ny + 4 * sys->numPorts, sizeof(double));
    j = 0;
    sys->yn[0] = yn[0];
    yi[sys->yn[0]] = j;
    for (i = 1; i <= county; i++){    // set the discretization length around port to be equal
        
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
    free(yn); yn = NULL;
    
    sort(zn, zn + countz + 1);
    zi.clear();
    sys->zn = (double*)calloc(sys->nz + 4 * sys->numPorts, sizeof(double));
    j = 0;
    sys->zn[0] = zn[0];
    zi[sys->zn[0]] = j;
    for (i = 1; i <= countz; i++){    // set the discretization length around port to be equal
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
    free(zn); zn = NULL;
    
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
        cout << "Index: " << i << " out of " << sys->numCdtRow << endl;
        numNode = (xi[sys->conductorIn[i].xmax] - xi[sys->conductorIn[i].xmin] + 1)
            *(yi[sys->conductorIn[i].ymax] - yi[sys->conductorIn[i].ymin] + 1)
            *(zi[sys->conductorIn[i].zmax] - zi[sys->conductorIn[i].zmin] + 1);
        //cout << xi[sys->conductorIn[i].xmax] << " " << xi[sys->conductorIn[i].xmin] << " " << yi[sys->conductorIn[i].ymax] << " " << yi[sys->conductorIn[i].ymin] << " " << zi[sys->conductorIn[i].zmax] << " " << zi[sys->conductorIn[i].zmin] << " " << numNode << endl;
        sys->conductorIn[i].cdtInNode = (int*)malloc(numNode*sizeof(int));
        sys->conductorIn[i].numNode = 0;

        for (j = xi[sys->conductorIn[i].xmin]; j <= xi[sys->conductorIn[i].xmax]; j++){
            cout << "Second index part 1: " << j << " out of " << xi[sys->conductorIn[i].xmax] << endl;
            cout << "Starting point of third index: " << yi[sys->conductorIn[i].ymin] << endl;
            cout << "Ending point of third index: " << yi[sys->conductorIn[i].ymax] << endl;
            for (k = yi[sys->conductorIn[i].ymin]; k <= yi[sys->conductorIn[i].ymax]; k++){
                cout << "Third index part 1: " << k << endl;
                if (polyIn(sys->xn[j], sys->yn[k], sys, i)){
                    cout << "polyIn == True" << endl;
                    for (m = zi[sys->conductorIn[i].zmin]; m < zi[sys->conductorIn[i].zmax]; m++){
                        //cout << sys->xn[j] << " " << sys->yn[k] << endl;
                        sys->conductorIn[i].cdtInNode[sys->conductorIn[i].numNode] = m*sys->N_node_s + (sys->N_cell_y + 1)*j + k;
                        sys->conductorIn[i].numNode++;
                        sys->markNode[m*sys->N_node_s + (sys->N_cell_y + 1)*j + k] = 1;
                        if (sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] == 0){   // set the z direction markEdge
                            sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + j*(sys->N_cell_y + 1) + k] = i + 1;    // mark this edge's corresponding index in conductorIn
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
            cout << "Second index part 2: " << j << " out of " << xi[sys->conductorIn[i].xmax] << endl;
            for (k = yi[sys->conductorIn[i].ymin]; k < yi[sys->conductorIn[i].ymax]; k++){
                for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
                    xc = sys->xn[j];
                    yc = (sys->yn[k] + sys->yn[k + 1]) / 2;
                    if (polyIn(xc, yc, sys, i)){
                        sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k] = i + 1;    // mark this edge's corresponding index in conductorIn
                        //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
                        //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + j*(sys->N_cell_y) + k << endl;
                    }
                }
            }
        }
        for (j = yi[sys->conductorIn[i].ymin]; j <= yi[sys->conductorIn[i].ymax]; j++){    // set the x direction markEdge
            cout << "Second index part 3: " << j << " out of " << yi[sys->conductorIn[i].ymax] << endl;
            for (k = xi[sys->conductorIn[i].xmin]; k < xi[sys->conductorIn[i].xmax]; k++){
                for (m = zi[sys->conductorIn[i].zmin]; m <= zi[sys->conductorIn[i].zmax]; m++){
                    xc = (sys->xn[k] + sys->xn[k + 1]) / 2;
                    yc = sys->yn[j];
                    if (polyIn(xc, yc, sys, i)){
                        sys->markEdge[m * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j] = i + 1;    // mark this edge's corresponding index in conductorIn
                        //cout << zi[sys->conductorIn[i].zmin] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
                        //cout << zi[sys->conductorIn[i].zmax] * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + k * (sys->N_cell_y + 1) + j << endl;
                    }
                }
            }
        }
    }
    
    /*fclose(cfp);*/
    cout << "Trying to close file" << endl;
    fclose(fp);
    cout << "File closed" << endl;
    
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
    for (i = 0; i < sys->N_node; i++){
        sys->nodeEdge.push_back(a);
        sys->nodeEdgea.push_back(a);
    }
    for (i = 0; i < sys->N_edge; i++){
        leng = pow((sys->nodepos[sys->edgelink[i * 2] * 3] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 1] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 1]), 2);
        leng = leng + pow((sys->nodepos[sys->edgelink[i * 2] * 3 + 2] - sys->nodepos[sys->edgelink[i * 2 + 1] * 3 + 2]), 2);
        leng = sqrt(leng);
        sys->nodeEdge[sys->edgelink[i * 2]].push_back(make_pair(i, 1 / leng));
        sys->nodeEdge[sys->edgelink[i * 2 + 1]].push_back(make_pair(i, -1 / leng));
    }

    int ix, iy, iz;
    iz = sys->N_cell_z;
    for (ix = 0; ix < sys->nx; ix++){
        for (iy = 0; iy < sys->ny; iy++){
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, -1 / (sys->zn[iz] - sys->zn[iz - 1])));
        }
    }
    iz = 0;
    for (ix = 0; ix < sys->nx; ix++){
        for (iy = 0; iy < sys->ny; iy++){
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->zn[iz + 1] - sys->zn[iz])));
        }
    }
    for (iz = 1; iz < sys->N_cell_z; iz++){
        for (ix = 0; ix < sys->nx; ix++){
            for (iy = 0; iy < sys->ny; iy++){
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair((iz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy , -2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->zn[iz + 1] - sys->zn[iz - 1])));
            }
        }
    }
    ix = sys->N_cell_x;
    for (iz = 0; iz < sys->nz; iz++){
        for (iy = 0; iy < sys->ny; iy++){
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -1 / (sys->xn[ix] - sys->xn[ix - 1])));
        }
    }
    ix = 0;
    for (iz = 0; iz < sys->nz; iz++){
        for (iy = 0; iy < sys->ny; iy++){
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 1 / (sys->xn[ix + 1] - sys->xn[ix])));
        }
    }
    for (ix = 1; ix < sys->N_cell_x; ix++){
        for (iz = 0; iz < sys->nz; iz++){
            for (iy = 0; iy < sys->ny; iy++){
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (ix - 1) * (sys->N_cell_y + 1) + iy, -2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + ix * (sys->N_cell_y + 1) + iy, 2 / (sys->xn[ix + 1] - sys->xn[ix - 1])));
            }
        }
    }
    iy = sys->N_cell_y;
    for (iz = 0; iz < sys->nz; iz++){
        for (ix = 0; ix < sys->nx; ix++){
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -1 / (sys->yn[iy] - sys->yn[iy - 1])));
        }
    }
    iy = 0;
    for (iz = 0; iz < sys->nz; iz++){
        for (ix = 0; ix < sys->nx; ix++){
            sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 1 / (sys->yn[iy + 1] - sys->yn[iy])));
        }
    }
    for (iy = 1; iy < sys->N_cell_y; iy++){
        for (iz = 0; iz < sys->nz; iz++){
            for (ix = 0; ix < sys->nx; ix++){
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy - 1, -2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
                sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].push_back(make_pair(iz * (sys->N_edge_s + sys->N_edge_v) + ix * sys->N_cell_y + iy, 2 / (sys->yn[iy + 1] - sys->yn[iy - 1])));
            }
        }
    }
    
    /* implement dfs */
    //cout <<"Number of nodes: " << sys->N_node << endl;
    int *visited;
    vector<int> st;
    unordered_set<int> base;
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
                sys->cond2condIn.push_back(base);
                while (!st.empty()){ 
                    mark = 0;
                    for (j = 0; j < sys->nodeEdge[st.back()].size(); j++){
                        if (sys->markEdge[sys->nodeEdge[st.back()][j].first] != 0){
                            if (sys->cond2condIn[count - 1].find(sys->markEdge[sys->nodeEdge[st.back()][j].first]) == sys->cond2condIn[count - 1].end()){
                                sys->cond2condIn[count - 1].insert(sys->markEdge[sys->nodeEdge[st.back()][j].first]);
                            }
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
        sys->conductor[i].markPort = 0;
    }
    //sys->markCondt = 0;
    for (i = 0; i < sys->N_node; i++){
        if (visited[i] != 0){
            sys->conductor[visited[i] - 1].node[sys->conductor[visited[i] - 1].cdtNodeind] = i;
            /*if (sys->markCondt == 0 && sys->nodepos[i * 3] >= 0.0022562 && sys->nodepos[i * 3] <= 0.00240265 && sys->nodepos[i * 3 + 1] >= 0.0037381 && sys->nodepos[i * 3 + 1] <= 0.00386135){
                sys->markCondt = visited[i];
            }*/
            if ((i < sys->N_node_s || i >= sys->N_node - sys->N_node_s) && sys->conductor[visited[i] - 1].markPort != -1){
                sys->conductor[visited[i] - 1].markPort = -1;
            }
            sys->conductor[visited[i] - 1].cdtNodeind++;
        }
    }
    //cout << "The added conductor is " << sys->markCondt << endl;
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
    /*outfile.open("eps.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_edge; i++){
        outfile << sys->eps[i] << endl;
    }
    outfile.close();*/

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

    /*outfile.open("sig.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_edge; i++){
        outfile << sys->sig[i] << endl;
    }
    outfile.close();*/

    sys->edgeCell.clear();
    sys->edgeCellArea.clear();
    //free(sys->markCell); sys->markCell = NULL;

    return 0;
}

int portSet(fdtdMesh* sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi){
    ofstream outfile;
    char s[FDTD_MAXC];
    char *word[FDTD_MAXC];
    int j, i, mark, l, m;
    int node;
    vector<int> edge;
    double a;
    double sideLen = 0;
    
    cout << "numCdtRow is " << sys->numCdtRow << endl;
    sys->exciteCdtLayer = (int*)calloc(sys->nz, sizeof(int));
    for (i = 0; i < sys->numPorts; i++)
    {
        if (sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] != 0 && sys->conductor[sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] - 1].markPort != -1){
            sys->portCoor[i].portCnd = sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]];
            sys->conductor[sys->markNode[zi[sys->portCoor[i].z1] * sys->N_node_s + xi[sys->portCoor[i].x1] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y1]] - 1].markPort = i + 1;    // markPort start from 1
        }
        else if (sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] != 0 && sys->conductor[sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] - 1].markPort != -1){
            sys->portCoor[i].portCnd = sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]];
            sys->conductor[sys->markNode[zi[sys->portCoor[i].z2] * sys->N_node_s + xi[sys->portCoor[i].x2] * (sys->N_cell_y + 1) + yi[sys->portCoor[i].y2]] - 1].markPort = i + 1;    // markPort start from 1
        }

        /*for (j = 0; j < sys->cdtNumNode[sys->portCoor[i].portCnd - 1]; j++){
            sys->exciteCdtLayer[sys->conductor[sys->portCoor[i].portCnd - 1].node[j] / sys->N_node_s] = 1;
        }*/
        
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
                a = a * (sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1);
            }
            else if (yi[sys->portCoor[i].y1] == sys->N_cell_y){
                a = a * (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]);
            }
            else{
                a = a * ((sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1) / 2 + (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]) / 2);
            }
            if (zi[sys->portCoor[i].z1] == 0){
                a = a *(sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
            }
            else if (zi[sys->portCoor[i].z1] == sys->N_cell_z){
                a = a * (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]);
            }
            else{
                a = a * ((sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1) / 2 + (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]) / 2);
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
                a = a * (sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1);
            }
            else if (xi[sys->portCoor[i].x1] == sys->N_cell_x){
                a = a * (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]);
            }
            else{
                a = a * ((sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1) / 2 + (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]) / 2);
            }
            if (zi[sys->portCoor[i].z1] == 0){
                a = a *(sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1);
            }
            else if (zi[sys->portCoor[i].z1] == sys->N_cell_z){
                a = a * (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]);
            }
            else{
                a = a * ((sys->zn[zi[sys->portCoor[i].z1] + 1] - sys->portCoor[i].z1) / 2 + (sys->portCoor[i].z1 - sys->zn[zi[sys->portCoor[i].z1] - 1]) / 2);
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
                a = a * (sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1);
            }
            else if (xi[sys->portCoor[i].x1] == sys->N_cell_x){
                a = a * (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]);
            }
            else{
                a = a * ((sys->xn[xi[sys->portCoor[i].x1] + 1] - sys->portCoor[i].x1) / 2 + (sys->portCoor[i].x1 - sys->xn[xi[sys->portCoor[i].x1] - 1]) / 2);
            }
            if (yi[sys->portCoor[i].y1] == 0){
                a = a *(sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1);
            }
            else if (yi[sys->portCoor[i].y1] == sys->N_cell_y){
                a = a * (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]);
            }
            else{
                a = a * ((sys->yn[yi[sys->portCoor[i].y1] + 1] - sys->portCoor[i].y1) / 2 + (sys->portCoor[i].y1 - sys->yn[yi[sys->portCoor[i].y1] - 1]) / 2);
            }
            sys->portArea.push_back(a);

        }
        sys->portEdge.push_back(edge);

    }
    
    clock_t t1 = clock();
    sys->markProSide = (int*)calloc(sys->N_node_s, sizeof(int));
    double x1, x2, y1, y2;
    int x1_ind, x2_ind, y1_ind, y2_ind;
    for (i = 0; i < sys->numPorts; i++){
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
                    for (j = x1_ind; j <= x2_ind; j++){
                        for (m = y1_ind; m <= y2_ind; m++){
                            sys->markProSide[j * (sys->N_cell_y + 1) + m] = 1;
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
                    for (j = x1_ind; j <= x2_ind; j++){
                        for (m = y1_ind; m <= y2_ind; m++){
                            sys->markProSide[j * (sys->N_cell_y + 1) + m] = 1;
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
                for (j = x1_ind; j <= x2_ind; j++){
                    for (m = y1_ind; m <= y2_ind; m++){
                        sys->markProSide[j * (sys->N_cell_y + 1) + m] = 1;
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
                for (j = x1_ind; j <= x2_ind; j++){
                    for (m = y1_ind; m <= y2_ind; m++){
                        sys->markProSide[j * (sys->N_cell_y + 1) + m] = 1;
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
