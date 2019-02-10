//#include "stdafx.h"
#include <ctime>
#include "fdtd.h"

static bool comp(pair<double, int> a, pair<double, int> b)
{
    return a.first <= b.first;
};

int paraGenerator(fdtdMesh *sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi){
    int i, j, mark, k, l, n;
    int status;
    int count;
    int xcol;
    ofstream outfile1;
    vector<int> rowId;
    vector<int> colId;
    vector<double> val;
    vector<int> temp2;
    int inx, iny, inz;


    cout << "Number of edges: " << sys->N_edge << endl;
    cout << "Number of nodes in one layer: " << sys->N_node_s << endl;
    /* Merged conductor node */



    

    /* Construct V0d with row id, col id and its val */
    int leng_v0d1 = 0, v0d1num = 0;    // store the num of v0d1 vectors, which are nodes outside the conductors
    int leng_v0d1a = 0, v0d1anum = 0;
    int leng_Ad = 0;
    int *map = (int*)calloc(sys->N_node, (sizeof(int)));
    double block1_x, block1_y, block2_x, block2_y, sideLen;

    block1_x = (sys->xn[sys->nx - 1] - sys->xn[0]) / 10;
    block1_y = (sys->yn[sys->ny - 1] - sys->yn[0]) / 10;
    block2_x = (sys->xn[sys->nx - 1] - sys->xn[0]) / 50;
    block2_y = (sys->yn[sys->ny - 1] - sys->yn[0]) / 50;


    status = merge_v0d1(sys, block1_x, block1_y, block2_x, block2_y, v0d1num, leng_v0d1, v0d1anum, leng_v0d1a);
    sys->v0d1valo = (double*)malloc(v0d1num * sizeof(double));
    for (i = 0; i < v0d1num; i++)
        sys->v0d1valo[i] = sys->v0d1val[i];
    for (i = 0; i < v0d1num; i++){    // compute sqrt(D_eps)*V0d1
        sys->v0d1val[i] = sys->v0d1val[i] * sqrt(sys->eps[sys->v0d1RowId[i]]);
    }
    sys->v0d1avalo = (double*)malloc(v0d1anum * sizeof(double));
    for (i = 0; i < v0d1anum; i++)
        sys->v0d1avalo[i] = sys->v0d1aval[i];
    for (i = 0; i < v0d1anum; i++){
        sys->v0d1aval[i] = sys->v0d1aval[i] * sqrt(sys->eps[sys->v0d1aRowId[i]]);
    }

    status = matrixMulti(sys->v0d1aColId, sys->v0d1aRowId, sys->v0d1aval, v0d1anum, sys->v0d1RowId, sys->v0d1ColId, sys->v0d1val, v0d1num, sys, leng_Ad, 1);
    if (status != 0)
        return status;
    

    //for (i = sys->N_node_s; i < sys->N_node - sys->N_node_s; i++){    // the upper and lower planes are PEC
    //    
    //    if (sys->markNode[i] == 0){
    //        for (j = 0; j < sys->nodeEdge[i].size(); j++){
    //            v0d1num++;
    //        }
    //        leng_v0d1++;
    //        
    //        inz = i / sys->N_node_s;
    //        inx = (i - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    //        iny = i % (sys->N_cell_y + 1);
    //        if (inz != 0 && inz != sys->N_cell_z){
    //            v0d1anum++;
    //            if (sys->markNode[i - sys->N_node_s] == 0 && ((i - sys->N_node_s) >= sys->N_node_s && (i - sys->N_node_s) < sys->N_node - sys->N_node_s)){    // the bottom plane is set to be ground
    //                leng_Ad++;
    //            }
    //        }
    //        else if (inz == sys->N_cell_z){
    //            v0d1anum++;
    //            if (sys->markNode[i - sys->N_node_s] == 0 && ((i - sys->N_node_s) >= sys->N_node_s && (i - sys->N_node_s) < sys->N_node-sys->N_node_s)){
    //                leng_Ad++;
    //            }
    //        }

    //        if (iny != 0 && iny != sys->N_cell_y){
    //            v0d1anum++;
    //            v0d1anum++;
    //        }
    //        else if (iny == sys->N_cell_y){
    //            v0d1anum++;
    //        }
    //        else{
    //            v0d1anum++;
    //        }

    //        if (inx != 0 && inx != sys->N_cell_x){
    //            v0d1anum++;
    //            v0d1anum++;
    //        }
    //        else if (inx == sys->N_cell_x){
    //            v0d1anum++;
    //        }
    //        else{
    //            v0d1anum++;
    //        }

    //        if (inz != sys->N_cell_z && inz != 0){
    //            v0d1anum++;
    //        }
    //        else if (inz == 0){
    //            v0d1anum++;
    //        }
    //        
    //        leng_v0d1a++;

    //        if (inx != 0 && inx != sys->N_cell_x){
    //            if (sys->markNode[i - sys->N_cell_y - 1] == 0 && ((i - sys->N_cell_y - 1) >= sys->N_node_s && (i - sys->N_cell_y - 1) < sys->N_node - sys->N_node_s)){
    //                leng_Ad++;
    //            }
    //            if (iny != 0 && iny != sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //            }
    //            else{
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //            }
    //            if (sys->markNode[i + sys->N_cell_y + 1] == 0 && ((i + sys->N_cell_y + 1) >= sys->N_node_s && (i + sys->N_cell_y + 1) < sys->N_node - sys->N_node_s)){
    //                leng_Ad++;
    //            }
    //        }
    //        else if (inx == sys->N_cell_x){
    //            if (sys->markNode[i - sys->N_cell_y - 1] == 0 && ((i - sys->N_cell_y - 1) >= sys->N_node_s && (i - sys->N_cell_y - 1) < sys->N_node - sys->N_node_s)){
    //                leng_Ad++;
    //            }
    //            if (iny != 0 && iny != sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //            }
    //            else{
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //            }
    //        }
    //        else{
    //            if (iny != 0 && iny != sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //            }
    //            else{
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    leng_Ad++;
    //                }
    //                else{
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0){
    //                    leng_Ad++;
    //                }
    //            }
    //            if (sys->markNode[i + sys->N_cell_y + 1] == 0 && ((i + sys->N_cell_y + 1) >= sys->N_node_s && (i + sys->N_cell_y + 1) < sys->N_node - sys->N_node_s)){
    //                leng_Ad++;
    //            }
    //        }

    //        if (inz != 0 && inz != sys->N_cell_z){
    //            if (sys->markNode[i + sys->N_node_s] == 0 && ((i + sys->N_node_s) >= sys->N_node_s && (i + sys->N_node_s) < sys->N_node - sys->N_node_s)){
    //                leng_Ad++;
    //            }
    //        }
    //        else if (inz == 0){
    //            if (sys->markNode[i + sys->N_node_s] == 0 && ((i + sys->N_node_s) >= sys->N_node_s && (i + sys->N_node_s) < sys->N_node - sys->N_node_s)){
    //                leng_Ad++;
    //            }
    //        }

    //    }
    //}
    //
    //sys->AdRowId = (int*)malloc(leng_Ad * sizeof(int));
    //sys->AdColId = (int*)malloc(leng_Ad * sizeof(int));
    //sys->Adval = (double*)malloc(leng_Ad * sizeof(double));
    //sys->v0d1RowId = (int*)malloc(v0d1num * sizeof(int));
    //sys->v0d1ColId = (int*)malloc(v0d1num * sizeof(int));
    //sys->v0d1val = (double*)malloc(v0d1num * sizeof(double));
    //sys->v0d1aRowId = (int*)malloc(v0d1anum * sizeof(int));
    //sys->v0d1aColId = (int*)malloc(v0d1anum* sizeof(int));
    //sys->v0d1aval = (double*)malloc(v0d1anum * sizeof(double));
    //
    //v0d1num = 0;
    //v0d1anum = 0;
    //leng_v0d1 = 0;
    //leng_v0d1a = 0;
    //leng_Ad = 0;

    //for (i = sys->N_node_s; i < sys->N_node - sys->N_node_s; i++){    // the upper and lower planes are PEC
    //    
    //    if (sys->markNode[i] == 0){
    //        for (j = 0; j < sys->nodeEdge[i].size(); j++){
    //            sys->v0d1RowId[v0d1num] = (sys->nodeEdge[i][j].first);
    //            sys->v0d1ColId[v0d1num] = (leng_v0d1);
    //            sys->v0d1val[v0d1num] = (sys->nodeEdge[i][j].second);
    //            v0d1num++;
    //        }
    //        map[i] = leng_v0d1 + 1;
    //        leng_v0d1++;
    //        
    //        inz = i / sys->N_node_s;
    //        inx = (i - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    //        iny = i % (sys->N_cell_y + 1);
    //        if (inz != 0 && inz != sys->N_cell_z){
    //            sys->v0d1aRowId[v0d1anum] = ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //            v0d1anum++;
    //            if (sys->markNode[i - sys->N_node_s] == 0 && ((i - sys->N_node_s) >= sys->N_node_s && (i - sys->N_node_s) < sys->N_node - sys->N_node_s)){    // the bottom plane is set to be ground
    //                sys->AdRowId[leng_Ad] = i;
    //                sys->AdColId[leng_Ad] = i - sys->N_node_s;
    //                sys->Adval[leng_Ad] = -2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny];
    //                leng_Ad++;
    //            }
    //        }
    //        else if (inz == sys->N_cell_z){
    //            sys->v0d1aRowId[v0d1anum] = ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    //            v0d1anum++;
    //            if (sys->markNode[i - sys->N_node_s] == 0 && ((i - sys->N_node_s) >= sys->N_node_s && (i - sys->N_node_s) < sys->N_node-sys->N_node_s)){
    //                sys->AdRowId[leng_Ad] = i;
    //                sys->AdColId[leng_Ad] = (i - sys->N_node_s);
    //                sys->Adval[leng_Ad] = (-1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                leng_Ad++;
    //            }
    //        }

    //        if (iny != 0 && iny != sys->N_cell_y){
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //            v0d1anum++;
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //            v0d1anum++;
    //        }
    //        else if (iny == sys->N_cell_y){
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    //            v0d1anum++;
    //        }
    //        else{
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (1 / (sys->yn[iny + 1] - sys->yn[iny]));
    //            v0d1anum++;
    //        }

    //        if (inx != 0 && inx != sys->N_cell_x){
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //            v0d1anum++;
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //            v0d1anum++;
    //            
    //        }
    //        else if (inx == sys->N_cell_x){
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    //            v0d1anum++;
    //            
    //        }
    //        else{
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (1 / (sys->xn[inx + 1] - sys->xn[inx]));
    //            v0d1anum++;
    //        }

    //        if (inz != sys->N_cell_z && inz != 0){
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //            v0d1anum++;
    //        }
    //        else if (inz == 0){
    //            sys->v0d1aRowId[v0d1anum] = (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //            sys->v0d1aColId[v0d1anum] = (leng_v0d1a);
    //            sys->v0d1aval[v0d1anum] = (1 / (sys->zn[inz + 1] - sys->zn[inz]));
    //            v0d1anum++;
    //        }
    //        
    //        leng_v0d1a++;

    //        if (inx != 0 && inx != sys->N_cell_x){
    //            if (sys->markNode[i - sys->N_cell_y - 1] == 0 && ((i - sys->N_cell_y - 1) >= sys->N_node_s && (i - sys->N_cell_y - 1) < sys->N_node - sys->N_node_s)){
    //                sys->AdRowId[leng_Ad] = (i);
    //                sys->AdColId[leng_Ad] = (i - sys->N_cell_y - 1);
    //                sys->Adval[leng_Ad] = (-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]);
    //                leng_Ad++;
    //            }
    //            if (iny != 0 && iny != sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i - 1);
    //                    sys->Adval[leng_Ad] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz])* sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i + 1);
    //                    sys->Adval[leng_Ad] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i - 1);
    //                    sys->Adval[leng_Ad] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            else{
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i + 1);
    //                    sys->Adval[leng_Ad] = (-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            if (sys->markNode[i + sys->N_cell_y + 1] == 0 && ((i + sys->N_cell_y + 1) >= sys->N_node_s && (i + sys->N_cell_y + 1) < sys->N_node - sys->N_node_s)){
    //                sys->AdRowId[leng_Ad] = (i);
    //                sys->AdColId[leng_Ad] = (i + sys->N_cell_y + 1);
    //                sys->Adval[leng_Ad] = (-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]);
    //                leng_Ad++;
    //            }
    //        }
    //        else if (inx == sys->N_cell_x){
    //            if (sys->markNode[i - sys->N_cell_y - 1] == 0 && ((i - sys->N_cell_y - 1) >= sys->N_node_s && (i - sys->N_cell_y - 1) < sys->N_node - sys->N_node_s)){
    //                sys->AdRowId[leng_Ad] = (i);
    //                sys->AdColId[leng_Ad] = (i - sys->N_cell_y - 1);
    //                sys->Adval[leng_Ad] = (-1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]);
    //                leng_Ad++;
    //            }
    //            if (iny != 0 && iny != sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i - 1);
    //                    sys->Adval[leng_Ad] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i + 1);
    //                    sys->Adval[leng_Ad] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i - 1);
    //                    sys->Adval[leng_Ad] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            else{
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i + 1);
    //                    sys->Adval[leng_Ad] = (-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //        }
    //        else{
    //            if (iny != 0 && iny != sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i - 1);
    //                    sys->Adval[leng_Ad] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0 && ((i + 1) >= sys->N_node_s && (i + 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i + 1);
    //                    sys->Adval[leng_Ad] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                if (sys->markNode[i - 1] == 0 && ((i - 1) >= sys->N_node_s && (i - 1) < sys->N_node - sys->N_node_s)){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i - 1);
    //                    sys->Adval[leng_Ad] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
    //                    leng_Ad++;
    //                }
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            else{
    //                if (inz != 0 && inz != sys->N_cell_z){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
    //                        + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else if (inz == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                else{
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i);
    //                    sys->Adval[leng_Ad] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
    //                        + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
    //                        + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->eps[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
    //                    leng_Ad++;
    //                }
    //                if (sys->markNode[i + 1] == 0){
    //                    sys->AdRowId[leng_Ad] = (i);
    //                    sys->AdColId[leng_Ad] = (i + 1);
    //                    sys->Adval[leng_Ad] = (-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
    //                    leng_Ad++;
    //                }
    //            }
    //            if (sys->markNode[i + sys->N_cell_y + 1] == 0 && ((i + sys->N_cell_y + 1) >= sys->N_node_s && (i + sys->N_cell_y + 1) < sys->N_node - sys->N_node_s)){
    //                sys->AdRowId[leng_Ad] = (i);
    //                sys->AdColId[leng_Ad] = (i + sys->N_cell_y + 1);
    //                sys->Adval[leng_Ad] = (-1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]);
    //                leng_Ad++;
    //            }
    //        }

    //        if (inz != 0 && inz != sys->N_cell_z){
    //            if (sys->markNode[i + sys->N_node_s] == 0 && ((i + sys->N_node_s) >= sys->N_node_s && (i + sys->N_node_s) < sys->N_node - sys->N_node_s)){
    //                sys->AdRowId[leng_Ad] = (i);
    //                sys->AdColId[leng_Ad] = (i + sys->N_node_s);
    //                sys->Adval[leng_Ad] = (-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny]);
    //                leng_Ad++;
    //            }
    //        }
    //        else if (inz == 0){
    //            if (sys->markNode[i + sys->N_node_s] == 0 && ((i + sys->N_node_s) >= sys->N_node_s && (i + sys->N_node_s) < sys->N_node - sys->N_node_s)){
    //                sys->AdRowId[leng_Ad] = (i);
    //                sys->AdColId[leng_Ad] = (i + sys->N_node_s);
    //                sys->Adval[leng_Ad] = (-1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->eps[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny]);
    //                leng_Ad++;
    //            }
    //        }

    //    }
    //}

   
    //for (i = 0; i < leng_Ad; i++){
    //    sys->AdRowId[i] = map[sys->AdRowId[i]] - 1;
    //    sys->AdColId[i] = map[sys->AdColId[i]] - 1;

    //}
    //outfile1.open("map.txt", std::ofstream::out | std::ofstream::trunc);
    //for (i = 0; i < sys->N_node; i++){
    //    outfile1 << map[i] << endl;
    //}
    //outfile1.close();
    //free(map); map = NULL;

    //
    //sys->v0d1valo = (double*)malloc(v0d1num * sizeof(double));
    //for (i = 0; i < v0d1num; i++)
    //    sys->v0d1valo[i] = sys->v0d1val[i];
    //for (i = 0; i < v0d1num; i++){    // compute sqrt(D_eps)*V0d1
    //    sys->v0d1val[i] = sys->v0d1val[i] * sqrt(sys->eps[sys->v0d1RowId[i]]);
    //}
    //sys->v0d1avalo = (double*)malloc(v0d1anum * sizeof(double));
    //for (i = 0; i < v0d1anum; i++)
    //    sys->v0d1avalo[i] = sys->v0d1aval[i];
    //for (i = 0; i < v0d1anum; i++){
    //    sys->v0d1aval[i] = sys->v0d1aval[i] * sqrt(sys->eps[sys->v0d1aRowId[i]]);
    //}

    

    cout << "Length of V0d1 is " << leng_v0d1 << endl;
    cout << "The number of nonzeros in Ad is " << leng_Ad << endl;

    int leng_v0d2 = 0, v0d2num = 0;
    int leng_v0d2a = 0, v0d2anum = 0;
    int condNodeNum = 0;
    for (i = 0; i < sys->numCdt; i++){
        
        for (j = 0; j < sys->cdtNumNode[i]; j++){
            if (sys->conductor[i].node[j] < sys->N_node_s || sys->conductor[i].node[j] >= sys->N_node - sys->N_node_s){
                sys->conductor[i].markPort = -1;
                break;
            }
        }
    }
    if (sys->numCdt > 0){
        
        int *ind;
        int indnum;
        for (i = 0; i < sys->numCdt; i++){
            if (sys->conductor[i].markPort != -1){
                ind = (int*)malloc(sys->cdtNumNode[i] * sizeof(int));
                indnum = 0;
                for (j = 0; j < sys->cdtNumNode[i]; j++){
                    ind[indnum] = (sys->conductor[i].node[j]);
                    indnum++;
                }

                status = nodeAdd_count(ind, indnum, sys->N_edge, sys, v0d2num, leng_v0d2);    // used to generate v0d2
                if (status != 0)
                    return status;

                status = nodeAddAvg_count(ind, indnum, sys->N_edge, sys, v0d2anum, leng_v0d2a);    // used to generate v0d2a
                if (status != 0)
                    return status;

                free(ind); ind = NULL;
            }
        }
        
        sys->v0d2RowId = (int*)malloc(v0d2num * sizeof(int));
        sys->v0d2ColId = (int*)malloc(v0d2num * sizeof(int));
        sys->v0d2val = (double*)malloc(v0d2num * sizeof(double));
        sys->v0d2aRowId = (int*)malloc(v0d2anum * sizeof(int));
        sys->v0d2aColId = (int*)malloc(v0d2anum * sizeof(int));
        sys->v0d2aval = (double*)malloc(v0d2anum * sizeof(double));
        v0d2num = 0;
        v0d2anum = 0;
        leng_v0d2 = 0;
        leng_v0d2a = 0;
        for (i = 0; i < sys->numCdt; i++){
            if (sys->conductor[i].markPort != -1){
                ind = (int*)malloc(sys->cdtNumNode[i] * sizeof(int));
                indnum = 0;
                for (j = 0; j < sys->cdtNumNode[i]; j++){
                    ind[indnum] = (sys->conductor[i].node[j]);
                    indnum++;
                }
                status = nodeAdd(ind, indnum, sys->N_edge, sys, v0d2num, leng_v0d2, 2);    // used to generate v0d2
                if (status != 0)
                    return status;

                status = nodeAddAvg(ind, indnum, sys->N_edge, sys, v0d2anum, leng_v0d2a, 2);    // used to generate v0d2a
                if (status != 0)
                    return status;
                free(ind); ind = NULL;
            }
        }
        
    }
    cout << "leng_v0d2 is " << leng_v0d2 << endl;

    /* Construct V0c with row id, col id and its val */
    int leng_v0c = 0, v0cnum = 0;
    int leng_v0ca = 0, v0canum = 0;

    int numPortCdt = 0;
    int leng_Ac = 0;
    count = 0;


    sys->cindex = (int*)calloc((sys->numCdt + 1), (sizeof(int)));
    sys->acu_cnno = (int*)calloc((sys->numCdt + 1), sizeof(int));
    sys->cindex[0] = -1;    // the last index in the sparse form for each conductor in V0c, i denotes ith conductor (starting from 1)
    sys->acu_cnno[0] = 0;    // how many v0c are in each conductor, i denotes the ith conductor (starting from 1)
    

    for (i = 0; i < sys->numCdt; i++){
        n = sys->cdtNumNode[i] - 1;
        if (sys->conductor[i].markPort == -1)
            n = sys->cdtNumNode[i];

        for (j = 0; j < n; j++){    // there is one node needed to be removed
            if (sys->conductor[i].node[j] < sys->N_node_s || sys->conductor[i].node[j] >= sys->N_node - sys->N_node_s){
                continue;
            }
            for (k = 0; k < sys->nodeEdge[sys->conductor[i].node[j]].size(); k++){
                v0cnum++;
                
            }

            sys->acu_cnno[count + 1]++;
            
            leng_v0c++;

            inz = sys->conductor[i].node[j] / sys->N_node_s;
            inx = (sys->conductor[i].node[j] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
            iny = sys->conductor[i].node[j] % (sys->N_cell_y + 1);
            if (inz != 0 && inz != sys->N_cell_z){
                v0canum++;
                if (sys->markEdge[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_node_s >= sys->N_node_s && sys->conductor[i].node[j] - sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1]))
                    {
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
            }
            else if (inz == sys->N_cell_z){
                v0canum++;
                if (sys->markEdge[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_node_s >= sys->N_node_s && sys->conductor[i].node[j] - sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
            }

            if (iny != 0 && iny != sys->N_cell_y){
                v0canum++;
                v0canum++;
            }
            else if (iny == sys->N_cell_y){
                v0canum++;
            }
            else{
                v0canum++;
            }

            if (inx != 0 && inx != sys->N_cell_x){
                v0canum++;
                v0canum++;
            }
            else if (inx == sys->N_cell_x){
                v0canum++;
            }
            else{
                v0canum++;
            }

            if (inz != sys->N_cell_z && inz != 0){
                v0canum++;
            }
            else if (inz == 0){
                v0canum++;
            }

            leng_v0ca++;
            if (inx != 0 && inx != sys->N_cell_x){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_cell_y - 1 >= sys->N_node_s && sys->conductor[i].node[j] - sys->N_cell_y - 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_cell_y - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
                if (iny != 0 && iny != sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                }
                else if (iny == sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
                else{
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                }
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny] != 0 && sys->conductor[i].node[j] + sys->N_cell_y + 1 >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_cell_y + 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_cell_y + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
            }
            else if (inx == sys->N_cell_x){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_cell_y - 1 >= sys->N_node_s && sys->conductor[i].node[j] - sys->N_cell_y - 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_cell_y - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
                if (iny != 0 && iny != sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                }
                else if (iny == sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
                else{
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                }
            }
            else{
                if (iny != 0 && iny != sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                }
                else if (iny == sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
                else{
                    if (inz != 0 && inz != sys->N_cell_z){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else if (inz == 0){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    else{
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            leng_Ac++;
                            sys->cindex[count + 1]++;
                        }
                    }
                }
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny] != 0 && sys->conductor[i].node[j] + sys->N_cell_y + 1 >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_cell_y + 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_cell_y + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
            }

            if (inz != 0 && inz != sys->N_cell_z){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] + sys->N_node_s >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
            }
            else if (inz == 0){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] + sys->N_node_s >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        leng_Ac++;
                        sys->cindex[count + 1]++;
                    }
                }
            }

        }
        sys->cindex[count + 1] = sys->cindex[count] + sys->cindex[count + 1];
        sys->acu_cnno[count + 1] = sys->acu_cnno[count] + sys->acu_cnno[count + 1];
        count++;
    }
    

    map = (int*)calloc(sys->N_node, (sizeof(int)));
    sys->AcRowId = (int*)malloc(leng_Ac * sizeof(int));
    sys->AcColId = (int*)malloc(leng_Ac * sizeof(int));
    sys->Acval = (double*)malloc(leng_Ac * sizeof(double));
    sys->v0cRowId = (int*)malloc(v0cnum * sizeof(int));
    sys->v0cColId = (int*)malloc(v0cnum * sizeof(int));
    sys->v0cval = (double*)malloc(v0cnum * sizeof(double));
    sys->v0caRowId = (int*)malloc(v0canum * sizeof(int));
    sys->v0caColId = (int*)malloc(v0canum * sizeof(int));
    sys->v0caval = (double*)malloc(v0canum * sizeof(double));
    leng_Ac = 0;
    v0cnum = 0;
    leng_v0c = 0;
    v0canum = 0;
    leng_v0ca = 0;

    for (i = 0; i < sys->numCdt; i++){
        n = sys->cdtNumNode[i] - 1;
        if (sys->conductor[i].markPort == -1)
            n = sys->cdtNumNode[i];

        for (j = 0; j < n; j++){    // there is one node needed to be removed
            if (sys->conductor[i].node[j] < sys->N_node_s || sys->conductor[i].node[j] >= sys->N_node - sys->N_node_s){
                continue;
            }
            for (k = 0; k < sys->nodeEdge[sys->conductor[i].node[j]].size(); k++){
                sys->v0cColId[v0cnum] = (leng_v0c);
                sys->v0cRowId[v0cnum] = (sys->nodeEdge[sys->conductor[i].node[j]][k].first);
                sys->v0cval[v0cnum] = (sys->nodeEdge[sys->conductor[i].node[j]][k].second);
                v0cnum++;
            }
            map[sys->conductor[i].node[j]] = leng_v0c + 1;

            leng_v0c++;

            inz = sys->conductor[i].node[j] / sys->N_node_s;
            inx = (sys->conductor[i].node[j] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
            iny = sys->conductor[i].node[j] % (sys->N_cell_y + 1);
            if (inz != 0 && inz != sys->N_cell_z){
                sys->v0caRowId[v0canum] = ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                v0canum++;
                if (sys->markEdge[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_node_s >= sys->N_node_s && sys->conductor[i].node[j] - sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1]))
                    {
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - sys->N_node_s);
                        sys->Acval[leng_Ac] = (-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                }
            }
            else if (inz == sys->N_cell_z){
                sys->v0caRowId[v0canum] = ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                v0canum++;
                if (sys->markEdge[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_node_s >= sys->N_node_s &&  sys->conductor[i].node[j] - sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - sys->N_node_s);
                        sys->Acval[leng_Ac] = (-1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                }
            }

            if (iny != 0 && iny != sys->N_cell_y){
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                v0canum++;
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                v0canum++;
            }
            else if (iny == sys->N_cell_y){
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                v0canum++;
            }
            else{
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (1 / (sys->yn[iny + 1] - sys->yn[iny]));
                v0canum++;
            }

            if (inx != 0 && inx != sys->N_cell_x){
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                v0canum++;
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                v0canum++;
            }
            else if (inx == sys->N_cell_x){
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                v0canum++;
            }
            else{
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (1 / (sys->xn[inx + 1] - sys->xn[inx]));
                v0canum++;
            }

            if (inz != sys->N_cell_z && inz != 0){
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                v0canum++;
            }
            else if (inz == 0){
                sys->v0caRowId[v0canum] = (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0caColId[v0canum] = (leng_v0ca);
                sys->v0caval[v0canum] = (1 / (sys->zn[inz + 1] - sys->zn[inz]));
                v0canum++;
            }

            leng_v0ca++;
            if (inx != 0 && inx != sys->N_cell_x){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_cell_y - 1 >= sys->N_node_s && sys->conductor[i].node[j] - sys->N_cell_y - 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_cell_y - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - sys->N_cell_y - 1);
                        sys->Acval[leng_Ac] = (-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]);
                        leng_Ac++;
                    }
                }
                if (iny != 0 && iny != sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - 1);
                            sys->Acval[leng_Ac] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                            leng_Ac++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz])* sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + 1);
                            sys->Acval[leng_Ac] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                            leng_Ac++;
                        }
                    }
                }
                else if (iny == sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - 1);
                            sys->Acval[leng_Ac] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                            leng_Ac++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                }
                else{
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + 1);
                            sys->Acval[leng_Ac] = (-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                            leng_Ac++;
                        }
                    }
                }
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny] != 0 && sys->conductor[i].node[j] + sys->N_cell_y + 1 >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_cell_y + 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_cell_y + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + sys->N_cell_y + 1);
                        sys->Acval[leng_Ac] = (-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]);
                        leng_Ac++;
                    }
                }
            }
            else if (inx == sys->N_cell_x){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny] != 0 && sys->conductor[i].node[j] - sys->N_cell_y - 1 >= sys->N_node_s && sys->conductor[i].node[j] - sys->N_cell_y - 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - sys->N_cell_y - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - sys->N_cell_y - 1);
                        sys->Acval[leng_Ac] = (-1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]);
                        leng_Ac++;
                    }
                }
                if (iny != 0 && iny != sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - 1);
                            sys->Acval[leng_Ac] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                            leng_Ac++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + 1);
                            sys->Acval[leng_Ac] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                            leng_Ac++;
                        }
                    }
                }
                else if (iny == sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - 1);
                            sys->Acval[leng_Ac] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                            leng_Ac++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                }
                else{
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + 1);
                            sys->Acval[leng_Ac] = (-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                            leng_Ac++;
                        }
                    }
                }
            }
            else{
                if (iny != 0 && iny != sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - 1);
                            sys->Acval[leng_Ac] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                            leng_Ac++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + 1);
                            sys->Acval[leng_Ac] = (-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                            leng_Ac++;
                        }
                    }
                }
                else if (iny == sys->N_cell_y){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->conductor[i].node[j] - 1 >= sys->N_node_s && sys->conductor[i].node[j] - 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] - 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] - 1);
                            sys->Acval[leng_Ac] = (-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                            leng_Ac++;
                        }
                    }
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                }
                else{
                    if (inz != 0 && inz != sys->N_cell_z){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                            + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else if (inz == 0){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    else{
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->Acval[leng_Ac] = (1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                            + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                            + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->conductor[i].node[j] + 1 >= sys->N_node_s && sys->conductor[i].node[j] + 1 < sys->N_node - sys->N_node_s){
                        if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                            sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                            sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + 1);
                            sys->Acval[leng_Ac] = (-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                            leng_Ac++;
                        }
                    }
                }
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny] != 0 && sys->conductor[i].node[j] + sys->N_cell_y + 1 >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_cell_y + 1 < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_cell_y + 1 == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + sys->N_cell_y + 1);
                        sys->Acval[leng_Ac] = (-1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]);
                        leng_Ac++;
                    }
                }
            }

            if (inz != 0 && inz != sys->N_cell_z){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] + sys->N_node_s >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + sys->N_node_s);
                        sys->Acval[leng_Ac] = (-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                }
            }
            else if (inz == 0){
                if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny] != 0 && sys->conductor[i].node[j] + sys->N_node_s >= sys->N_node_s && sys->conductor[i].node[j] + sys->N_node_s < sys->N_node - sys->N_node_s){
                    if (!(n == sys->cdtNumNode[i] - 1 && sys->conductor[i].node[j] + sys->N_node_s == sys->conductor[i].node[sys->cdtNumNode[i] - 1])){
                        sys->AcRowId[leng_Ac] = (sys->conductor[i].node[j]);
                        sys->AcColId[leng_Ac] = (sys->conductor[i].node[j] + sys->N_node_s);
                        sys->Acval[leng_Ac] = (-1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny]);
                        leng_Ac++;
                    }
                }
            }

        }
    }


    cout << "The storage of Ac " << leng_Ac << endl;
    for (i = 0; i < leng_Ac; i++){
        sys->AcRowId[i] = map[sys->AcRowId[i]] - 1;
        sys->AcColId[i] = map[sys->AcColId[i]] - 1;

    }
    free(map); map = NULL;

    sys->v0cvalo = (double*)malloc(v0cnum * sizeof(double));
    for (i = 0; i < v0cnum; i++)
        sys->v0cvalo[i] = sys->v0cval[i];    // v0cvalo is the v0c values without D_sig
    for (i = 0; i < v0cnum; i++){
        sys->v0cval[i] = sys->v0cval[i] * sqrt(sys->sig[sys->v0cRowId[i]]);       // Compute the sparse form of D_sig*V0c
    }
    sys->v0cavalo = (double*)malloc(v0canum * sizeof(double));
    for (i = 0; i < v0canum; i++)
        sys->v0cavalo[i] = sys->v0caval[i];
    for (i = 0; i < v0canum; i++){
        sys->v0caval[i] = sys->v0caval[i] * sqrt(sys->sig[sys->v0caRowId[i]]);
    }



    //count = 0;
    //
    //block1_x = (sys->xn[sys->nx - 1] - sys->xn[0]) / 10;
    //block1_y = (sys->yn[sys->ny - 1] - sys->yn[0]) / 10;
    //block2_x = 0;// (sys->xn[sys->nx - 1] - sys->xn[0]) / 250;
    //block2_y = 0;// (sys->yn[sys->ny - 1] - sys->yn[0]) / 250;

    //status = merge_v0c(sys, block1_x, block1_y, block2_x, block2_y, v0cnum, leng_v0c, v0canum, leng_v0ca);


    //sys->v0cvalo = (double*)malloc(v0cnum * sizeof(double));
    //for (i = 0; i < v0cnum; i++)
    //    sys->v0cvalo[i] = sys->v0cval[i];    // v0cvalo is the v0c values without D_sig
    //for (i = 0; i < v0cnum; i++){
    //    sys->v0cval[i] = sys->v0cval[i] * sqrt(sys->sig[sys->v0cRowId[i]]);       // Compute the sparse form of D_sig*V0c
    //}
    //sys->v0cavalo = (double*)malloc(v0canum * sizeof(double));
    //for (i = 0; i < v0canum; i++)
    //    sys->v0cavalo[i] = sys->v0caval[i];
    //for (i = 0; i < v0canum; i++){
    //    sys->v0caval[i] = sys->v0caval[i] * sqrt(sys->sig[sys->v0caRowId[i]]);
    //}

    //status = matrixMulti(sys->v0caColId, sys->v0caRowId, sys->v0caval, v0canum, sys->v0cRowId, sys->v0cColId, sys->v0cval, v0cnum, sys, leng_Ac, 2);    // the last variable 1 is for v0d, 2 is for v0c
    //if (status != 0)
    //    return status;

    cout << "leng v0c " << leng_v0c << endl;
    cout << "The number of nonzeros in Ac is " << leng_Ac << endl;
    /*for (i = 0; i < sys->numCdt + 1; i++){
        cout << sys->acu_cnno[i] << " ";
    }
    cout << endl;
    for (i = 0; i < sys->numCdt + 1; i++){
        cout << sys->cindex[i] << " ";
    }
    cout << endl;
    */
    

    /*sys->AcRowId1 = (int*)malloc((leng_v0c + 1) * sizeof(int));
    status = COO2CSR_malloc(sys->AcRowId, sys->AcColId, sys->Acval, leng_Ac, leng_v0c, sys->AcRowId1);
    if (status != 0)
    return status;

    free(sys->AcRowId); sys->AcRowId = NULL;

    double *a = &(sys->Acval[0]);
    int *ia = &(sys->AcRowId1[0]);
    int *ja = &(sys->AcColId[0]);*/


    //outfile1.open("v0c.txt", std::ofstream::out | std::ofstream::trunc);
    //for (i = 0; i < v0cnum; i++){
    //outfile1 << sys->v0cRowId[i] + 1 << " " << sys->v0cColId[i] + 1 << " " << sys->v0cvalo[i] << endl;
    //}
    //outfile1.close();
    //outfile1.open("v0ca.txt", std::ofstream::out | std::ofstream::trunc);
    //for (i = 0; i < v0canum; i++){
    //outfile1 << sys->v0caRowId[i] + 1 << " " << sys->v0caColId[i] + 1 << " " << sys->v0cavalo[i] << endl;
    //}
    //outfile1.close();

    /* Compute the matrix V0c'*D_sig*V0c */
    /*status = matrixMulti(sys->v0caColId, sys->v0caRowId, sys->v0caval, sys->v0cRowId, sys->v0cColId, sys->v0cval, sys->AcRowId, sys->AcColId, sys->Acval);
    if (status != 0)
    return status;*/


    /*outfile1.open("v0d1.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d1num; i++){
        outfile1 << sys->v0d1RowId[i] + 1 << " " << sys->v0d1ColId[i] + 1 << " " << sys->v0d1val[i] << endl;
    }
    outfile1.close();
    outfile1.open("v0d1a.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d1anum; i++){
        outfile1 << sys->v0d1aRowId[i] + 1 << " " << sys->v0d1aColId[i] + 1 << " " << sys->v0d1aval[i] << endl;
    }
    outfile1.close();
    outfile1.open("v0d2.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d2num; i++){
        outfile1 << sys->v0d2RowId[i] + 1 << " " << sys->v0d2ColId[i] + 1 << " " << sys->v0d2val[i] << endl;
    }
    outfile1.close();
    outfile1.open("v0d2a.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d2anum; i++){
        outfile1 << sys->v0d2aRowId[i] + 1 << " " << sys->v0d2aColId[i] + 1 << " " << sys->v0d2aval[i] << endl;
    }
    outfile1.close();*/






    sys->AdRowId1 = (int*)malloc((leng_v0d1 + 1) * sizeof(int));

    status = COO2CSR_malloc(sys->AdRowId, sys->AdColId, sys->Adval, leng_Ad, leng_v0d1, sys->AdRowId1);
    if (status != 0)
        return status;
    free(sys->AdRowId); sys->AdRowId = NULL;
    

   /* double *d = &(sys->Adval[0]);
    int *id = &(sys->AdRowId1[0]);
    int *jd = &(sys->AdColId[0]);*/




    

    if (leng_v0d2 > 0){
        sys->v0d2valo = (double*)malloc(v0d2num * sizeof(double));
        for (i = 0; i < v0d2num; i++)
            sys->v0d2valo[i] = sys->v0d2val[i];
        for (i = 0; i < v0d2num; i++){    // compute sqrt(D_eps)*v0d2
            sys->v0d2val[i] = sys->v0d2val[i] * sqrt(sys->eps[sys->v0d2RowId[i]]);
        }
        sys->v0d2avalo = (double*)malloc(v0d2anum * sizeof(double));
        for (i = 0; i < v0d2anum; i++)
            sys->v0d2avalo[i] = sys->v0d2aval[i];
        for (i = 0; i < v0d2anum; i++){
            sys->v0d2aval[i] = sys->v0d2aval[i] * sqrt(sys->eps[sys->v0d2aRowId[i]]);
        }
        sys->v0d2ColIdo = (int*)malloc(v0d2num * sizeof(int));
        for (i = 0; i < v0d2num; i++)
            sys->v0d2ColIdo[i] = sys->v0d2ColId[i];
        //free(sys->v0d2ColId); 
        sys->v0d2ColId = (int*)malloc((leng_v0d2 + 1)*sizeof(int));
        status = COO2CSR_malloc(sys->v0d2ColIdo, sys->v0d2RowId, sys->v0d2val, v0d2num, leng_v0d2, sys->v0d2ColId);
        if (status != 0)
            return status;
        sys->v0d2aColIdo = (int*)malloc(v0d2anum * sizeof(int));
        for (i = 0; i < v0d2anum; i++)
            sys->v0d2aColIdo[i] = sys->v0d2aColId[i];
        free(sys->v0d2aColId); sys->v0d2aColId = (int*)malloc((leng_v0d2a + 1) * sizeof(int));
        status = COO2CSR_malloc(sys->v0d2aColIdo, sys->v0d2aRowId, sys->v0d2aval, v0d2anum, leng_v0d2a, sys->v0d2aColId);
        if (status != 0)
            return status;
    }
    int leng_v0d = leng_v0d1 + leng_v0d2;
    sys->v0d1ColIdo = (int*)malloc(v0d1num * sizeof(int));
    for (i = 0; i < v0d1num; i++)
        sys->v0d1ColIdo[i] = sys->v0d1ColId[i];
    free(sys->v0d1ColId); sys->v0d1ColId = (int*)malloc((leng_v0d1 + 1) * sizeof(int));
    status = COO2CSR_malloc(sys->v0d1ColIdo, sys->v0d1RowId, sys->v0d1val, v0d1num, leng_v0d1, sys->v0d1ColId);
    if (status != 0)
        return status;
    sys->v0d1aColIdo = (int*)malloc(v0d1anum * sizeof(int));
    for (i = 0; i < v0d1anum; i++)
        sys->v0d1aColIdo[i] = sys->v0d1aColId[i];
    free(sys->v0d1aColId); sys->v0d1aColId = (int*)malloc((leng_v0d1a + 1) * sizeof(int));
    status = COO2CSR_malloc(sys->v0d1aColIdo, sys->v0d1aRowId, sys->v0d1aval, v0d1anum, leng_v0d1a, sys->v0d1aColId);
    if (status != 0)
        return status;
    sys->v0cColIdo = (int*)malloc(v0cnum * sizeof(int));
    for (i = 0; i < v0cnum; i++)
        sys->v0cColIdo[i] = sys->v0cColId[i];
    free(sys->v0cColId); sys->v0cColId = (int*)malloc((leng_v0c + 1) * sizeof(int));
    status = COO2CSR_malloc(sys->v0cColIdo, sys->v0cRowId, sys->v0cval, v0cnum, leng_v0c, sys->v0cColId);
    if (status != 0)
        return status;
    sys->v0caColIdo = (int*)malloc(v0canum * sizeof(int));
    for (i = 0; i < v0canum; i++)
        sys->v0caColIdo[i] = sys->v0caColId[i];
    free(sys->v0caColId); sys->v0caColId = (int*)malloc((leng_v0ca + 1)*sizeof(int));
    status = COO2CSR_malloc(sys->v0caColIdo, sys->v0caRowId, sys->v0caval, v0canum, leng_v0ca, sys->v0caColId);
    if (status != 0)
        return status;
    /* Pick up a y0c2 cooresponding to one source port */
    ofstream outfile;
    complex<double> Zresult;
    double *bd1, *bd2;
    double *bdc1, *bdc2;
    double *temp;
    double mdone;
    int ione;
    double *solution_d2;
    double *v0d2, *b0d2;
    double *v0d2epsv0d2, *v0d2epsv0d1;
    int *ipiv;
    int info;
    double *workspace;
    double *xd2;
    double *temp1;
    int startCol;
    startCol = 0;
    sys->x = (complex<double>*) malloc(sys->numPorts * sys->numPorts * sizeof(complex<double>));
    for (i = 0; i < sys->numPorts*sys->numPorts; i++){
        sys->x[i] = complex<double>(0., 0.); // Complex double constructor from real and imaginary
    }
    

    double *b, *xc;    // the array of the right hand side
    xcol = 0;

    vector<int> v0csRowId, v0csColId;
    vector<double> v0csval;
    int leng_v0cs;
    double *epsy, *sigy;
    double *Jr, *Ji;
    int port, sourcePort;    // show which port it is using
    int node;
    double v0csedgeleng;
    complex<double> current;
    int port_N_node_s;
    double *crhs;
    double *yc_eps;
    sys->Y = (complex<double>*)calloc(sys->numPorts * sys->numPorts * sys->nfreq, sizeof(complex<double>));

    

    char transa;
    int m;
    double alpha = 1, beta = 0;
    char matdescra[6];
    matdescra[0] = 'G'; matdescra[3] = 'C';    // general matrix multi, 0-based indexing
    cout << "Begin!\n";
    int nrhs = leng_v0d2;
    v0d2epsv0d2 = (double*)calloc(nrhs*nrhs, sizeof(double));
    v0d2epsv0d1 = (double*)calloc(nrhs*nrhs, sizeof(double));
    solution_d2 = (double*)calloc(leng_v0d2*leng_v0d1, sizeof(double));
    v0d2 = (double*)calloc(leng_v0d2*sys->N_edge, sizeof(double));
    b0d2 = (double*)calloc(leng_v0d2*leng_v0d1, sizeof(double));
    if (leng_v0d2 != 0){
        transa = 'N';
        m = leng_v0d1;
        k = sys->N_edge;

        for (i = 0; i < v0d2num; i++){
            v0d2[sys->v0d2ColIdo[i] * sys->N_edge + sys->v0d2RowId[i]] = sys->v0d2val[i];
        }

        for (i = 1; i <= leng_v0d2; i++){
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1aval[0], &sys->v0d1aRowId[0], &sys->v0d1aColId[0], &sys->v0d1aColId[1], &v0d2[(i - 1)*sys->N_edge], &beta, &b0d2[(i - 1)*leng_v0d1]);
        }
        /*for (i = 0; i < leng_v0d1; i++){
            for (j = 0; j < leng_v0d2; j++){
                cout << b0d2[j * leng_v0d1 + i] << " ";
            }
            cout << endl;
        }*/
        
        for (i = 0; i < leng_v0d2; i++){
            //clock_t t1 = clock();
            /*status = interativeSolver(leng_v0d1, 1, &b0d2[i*leng_v0d1], &sys->v0d1ColId[0], &sys->v0d1RowId[0], &sys->v0d1val[0], &sys->v0d1aColId[0], &sys->v0d1aRowId[0], &sys->v0d1aval[0], &solution_d2[i*leng_v0d1], sys);
            if (status != 0)
                return status;*/
            status = pardisoSolve(sys, &b0d2[i*leng_v0d1], &solution_d2[i*leng_v0d1], leng_v0d1);
            //cout << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
        }
        /*for (i = 0; i < leng_v0d1; i++){
            for (j = 0; j < leng_v0d2; j++){
                cout << solution_d2[j*leng_v0d1 + i] << " ";
            }
            cout << endl;
        }*/
        cout << "The first iteration is done!\n";

        temp1 = (double*)malloc(sys->N_edge*nrhs*sizeof(double));
        m = nrhs;
        transa = 'N';
        for (i = 0; i < leng_v0d2; i++){
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2aval[0], &sys->v0d2aRowId[0], &sys->v0d2aColId[0], &sys->v0d2aColId[1], &v0d2[i*sys->N_edge], &beta, &v0d2epsv0d2[i*nrhs]);
        }
        transa = 'T';
        m = leng_v0d1;
        for (i = 0; i < nrhs; i++){
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1val[0], &sys->v0d1RowId[0], &sys->v0d1ColId[0], &sys->v0d1ColId[1], &solution_d2[i*leng_v0d1], &beta, &temp1[i*sys->N_edge]);
        }
        transa = 'N';
        m = nrhs;
        for (i = 0; i < nrhs; i++){
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2aval[0], &sys->v0d2aRowId[0], &sys->v0d2aColId[0], &sys->v0d2aColId[1], &temp1[i*sys->N_edge], &beta, &v0d2epsv0d1[i*nrhs]);
        }
        mdone = -1;
        ione = 1;
        cblas_daxpy(nrhs*nrhs, mdone, v0d2epsv0d1, ione, v0d2epsv0d2, ione);
        info = 0;
        ipiv = (int*)malloc(leng_v0d2*sizeof(int));
        workspace = (double*)malloc(leng_v0d2*sizeof(double));
        dgetrf(&leng_v0d2, &leng_v0d2, v0d2epsv0d2, &leng_v0d2, ipiv, &info);
        dgetri(&leng_v0d2, v0d2epsv0d2, &leng_v0d2, ipiv, workspace, &leng_v0d2, &info);
        /*for (i = 0; i < leng_v0d2; i++){
            for (j = 0; j < leng_v0d2; j++){
                cout << v0d2epsv0d2[j*leng_v0d2 + i] << " ";
            }
            cout << endl;
        }*/
    }

    double *ydcp;
    double *y0c;
    double *yc;
    double *yccp;
    double *dRhs;
    complex<double> *yd, *yd2;
    double *v0caJ;
    double leng;
    sourcePort = 0;
   
    while (sourcePort < sys->numPorts){

        sys->J = (double*)calloc(sys->N_edge, sizeof(double));
        for (i = 0; i < sys->portEdge[sourcePort].size(); i++){
            sys->J[sys->portEdge[sourcePort][i]] = sys->portCoor[sourcePort].portDirection;
        }

        dRhs = (double*)malloc(sys->N_edge*sizeof(double));
        for (i = 0; i < sys->N_edge; i++){
            dRhs[i] = -sys->J[i];
        }
        
        yd = (complex<double>*)malloc((sys->N_edge)*sizeof(complex<double>));
        status = solveV0dSystem(sys, dRhs, yd, v0d2epsv0d2, solution_d2, leng_v0d1, leng_v0d2);    // dRhs is the right hand side haven't multiply V0da^T
        
        /* Compute C right hand side */
        y0c = (double*)malloc(leng_v0c * sizeof(double));


        //cout << (clock() - t) * 1.0 / CLOCKS_PER_SEC << endl;
        v0caJ = (double*)calloc(leng_v0c, sizeof(double));
        transa = 'N';
        m = leng_v0ca;
        k = sys->N_edge;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0cavalo[0], &sys->v0caRowId[0], &sys->v0caColId[0], &sys->v0caColId[1], sys->J, &beta, v0caJ);
        for (i = 0; i < leng_v0c; i++){
            v0caJ[i] = -v0caJ[i];
        }
        crhs = (double*)malloc(leng_v0c*sizeof(double));
        ydcp = (double*)calloc(sys->N_edge, sizeof(double));
        for (i = 0; i < sys->N_edge; i++){
            ydcp[i] = yd[i].imag() * (2 * PI*sys->freqStart * sys->freqUnit)*sys->eps[i];
        }

        transa = 'N';
        m = leng_v0c;
        k = sys->N_edge;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0cavalo[0], &sys->v0caRowId[0], &sys->v0caColId[0], &sys->v0caColId[1], ydcp, &beta, crhs);

        for (i = 0; i < leng_v0c; i++){
            v0caJ[i] = (v0caJ[i] + crhs[i]);
        }

        cout << "h\n";
        /*solve c system block by block*/
        for (i = 1; i <= sys->numCdt; i++){
            pardisoSolve_c(sys, &v0caJ[sys->acu_cnno[i - 1]], &y0c[sys->acu_cnno[i - 1]], sys->acu_cnno[i - 1], sys->acu_cnno[i] - 1, sys->cindex[i - 1] + 1, sys->cindex[i]);
        }
        //pardisoSolve_c(sys, &v0caJ[0], &y0c[0], 0, leng_v0c - 1, 0, leng_Ac - 1);
        cout << "h\n";

        /* V0cy0c */
        yc = (double*)malloc(sys->N_edge * sizeof(double));
        yccp = (double*)malloc(sys->N_edge * sizeof(double));
        transa = 'T';
        k = sys->N_edge;
        m = leng_v0c;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0cvalo[0], &sys->v0cRowId[0], &sys->v0cColId[0], &sys->v0cColId[1], y0c, &beta, yc);

        for (i = 0; i < sys->N_edge; i++){
            yccp[i] = -yc[i] * 2 * PI*sys->freqStart*sys->freqUnit*sys->eps[i];

        }

        yd2 = (complex<double>*)malloc(sys->N_edge*sizeof(complex<double>));
        status = solveV0dSystem(sys, yccp, yd2, v0d2epsv0d2, solution_d2, leng_v0d1, leng_v0d2);    // dRhs is the right hand side haven't multiply V0da^T
        for (i = 0; i < sys->N_edge; i++){
            yd2[i] = -yd2[i].imag();
        }

        for (int id = 0; id < sys->N_edge; id++){
            yd[id] = yd[id].real() + yd2[id].real() + (1i)*(yd2[id].imag() + yd[id].imag());
        }


        sys->y = (complex<double>*)malloc(sys->N_edge*sizeof(complex<double>));
        for (i = 0; i < sys->N_edge; i++){
            sys->y[i] = yd[i].real() + yc[i] + (1i) * yd[i].imag();
        }
        for (i = 0; i < sys->numPorts; i++){
            for (j = 0; j < sys->portEdge[i].size(); j++){
                leng = pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3]), 2);
                leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3 + 1]), 2);
                leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3 + 2]), 2);
                leng = sqrt(leng);
                sys->x[i + sys->numPorts*xcol] = (sys->x[i + sys->numPorts*xcol].real() + sys->y[sys->portEdge[i][j]].real() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection))) + (1i)*(sys->y[sys->portEdge[i][j]].imag() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection)) + sys->x[i + sys->numPorts*xcol].imag());

            }
        }
        free(y0c); y0c = NULL;
        free(yd); yd = NULL;
        free(ydcp); ydcp = NULL;
        free(yc); yc = NULL;
        free(yccp); yccp = NULL;
        free(v0caJ); v0caJ = NULL;
        free(sys->y); sys->y = NULL;
        free(dRhs); dRhs = NULL;
        free(sys->J); sys->J = NULL;
        free(crhs); crhs = NULL;

        sourcePort++;
        xcol++;
        
    }
    for (int id = 0; id < sys->nfreq; id++){
        cout << "Z parameter at frequency " << (sys->freqStart + id * (sys->freqEnd - sys->freqStart) / (sys->nfreq - 1)) * sys->freqUnit << " is " << endl;
        for (i = 0; i < sys->numPorts; i++){
            for (j = 0; j < sys->numPorts; j++){
                Zresult = sys->x[j + i*sys->numPorts].real() + (1i) * sys->x[j + i*sys->numPorts].imag() * sys->freqStart / (sys->freqStart + id * (sys->freqEnd - sys->freqStart) / (sys->nfreq - 1));
                cout << Zresult << "\n";
            }
        }
    }
    free(v0d2epsv0d2); v0d2epsv0d2 = NULL;
    free(v0d2epsv0d1); v0d2epsv0d1 = NULL;
    free(v0d2); v0d2 = NULL;
    free(b0d2); b0d2 = NULL;
    free(solution_d2); solution_d2 = NULL;
    free(sys->cindex); sys->cindex = NULL;
    free(sys->acu_cnno); sys->acu_cnno = NULL;

    free(sys->AdColId); sys->AdColId = NULL;
    free(sys->Adval); sys->Adval = NULL;
    free(sys->AdRowId1); sys->AdRowId1 = NULL;
    free(sys->v0d1RowId); sys->v0d1RowId = NULL;
    free(sys->v0d1ColId); sys->v0d1ColId = NULL;
    free(sys->v0d1ColIdo); sys->v0d1ColIdo = NULL;
    free(sys->v0d1val); sys->v0d1val = NULL;
    free(sys->v0d1valo); sys->v0d1valo = NULL;
    free(sys->v0d1aRowId); sys->v0d1aRowId = NULL;
    free(sys->v0d1aColId); sys->v0d1aColId = NULL;
    free(sys->v0d1aColIdo); sys->v0d1aColIdo = NULL;
    free(sys->v0d1aval); sys->v0d1aval = NULL;
    free(sys->v0d1avalo); sys->v0d1avalo = NULL;
    if (leng_v0d2 > 0){
        free(sys->v0d2RowId);  sys->v0d2RowId = NULL;
        free(sys->v0d2ColId);  sys->v0d2ColId = NULL;
        free(sys->v0d2val); sys->v0d2val = NULL;
        free(sys->v0d2valo); sys->v0d2valo = NULL;
        free(sys->v0d2ColIdo);  sys->v0d2ColIdo = NULL;
        free(sys->v0d2aRowId); sys->v0d2aRowId = NULL;
        free(sys->v0d2aColId); sys->v0d2aColId = NULL;
        free(sys->v0d2aval); sys->v0d2aval = NULL;
    }
    free(sys->v0cRowId);  sys->v0cRowId = NULL;
    free(sys->v0cColId);  sys->v0cColId = NULL;
    free(sys->v0cColIdo); sys->v0cColIdo = NULL;
    free(sys->v0cval); sys->v0cval = NULL;
    free(sys->v0cvalo); sys->v0cvalo = NULL;
    free(sys->v0caRowId); sys->v0caRowId = NULL;
    free(sys->v0caColId); sys->v0caColId = NULL;
    free(sys->v0caColIdo); sys->v0caColIdo = NULL;
    free(sys->v0caval); sys->v0caval = NULL;
    free(sys->v0cavalo); sys->v0cavalo = NULL;
    /*free(sys->AcRowId1); sys->AcRowId1 = NULL;
    free(sys->AcColId); sys->AcColId = NULL;
    free(sys->Acval); sys->Acval = NULL;*/
    

    return 0;
}
int pardisoSolve_c(fdtdMesh *sys, double *rhs, double *solution, int nodestart, int nodeend, int indstart, int indend){
    /*solve Ac system block by block*/
    int leng = nodeend - nodestart + 1;
    int status;
    double *a = (double*)malloc((indend - indstart + 1) * sizeof(double));
    int *ia = (int*)malloc((indend - indstart + 1) * sizeof(int));
    int *ja = (int*)malloc((indend - indstart + 1) * sizeof(int));

    
    for (int i = 0; i < (indend - indstart + 1); i++){
        a[i] = sys->Acval[indstart + i];
        ia[i] = sys->AcRowId[indstart + i] - sys->AcRowId[indstart];
        ja[i] = sys->AcColId[indstart + i] - sys->AcColId[indstart];
        
    }


    int *ia1 = (int*)malloc((leng + 1) * sizeof(int));
   /* status = COO2CSR_malloc(ia, ja, a, indstart - indend + 1, leng, ia1);
    if (status != 0)
        return status;*/
    int count = 0;
    int i = 0;
    int k = 0;
    int start;
    ia1[k] = 0;
    k++;
    while (i < (indend - indstart + 1)){
        start = ia[i];
        while (i < (indend - indstart + 1) && ia[i] == start) {
            count++;
            i++;
        }
        ia1[k] = (count);
        k++;
    }


    
    

    void *pt[64];
    int mtype;
    int iparm[64];
    double dparm[64];
    int maxfct, mnum, phase, error, solver;
    int num_process;   //number of processors
    int v0csin;
    int perm;
    int nrhs = 1;
    int msglvl = 0;    //print statistical information

    mtype = 11;    // real and not symmetric
    solver = 0;
    error = 0;
    maxfct = 1;    //maximum number of numerical factorizations
    mnum = 1;    //which factorization to use
    phase = 13;    //analysis

    pardisoinit(pt, &mtype, iparm);
    nrhs = 1;
    iparm[38] = 1;
    iparm[34] = 1;    //0-based indexing
    
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &leng, a, ia1, ja, &perm, &nrhs, iparm, &msglvl, rhs, solution, &error);
    
    free(a); a = NULL;
    free(ia); ia = NULL;
    free(ja); ja = NULL;
    free(ia1); ia1 = NULL;
    

    return 0;
}

int pardisoSolve(fdtdMesh *sys, double *rhs, double *solution, int leng_v0d1){
    /* A\b1 */
    
    double *d = &(sys->Adval[0]);
    int *id = &(sys->AdRowId1[0]);
    int *jd = &(sys->AdColId[0]);

    void *ptd[64];
    int mtyped;
    int iparmd[64];
    double dparmd[64];
    int maxfctd, mnumd, phased, errord, solverd;
    int num_processd;   //number of processors
    int v0csin;
    int permd;
    int nrhs = 1;
    int msglvld = 0;    //print statistical information

    mtyped = 11;    // real and not symmetric
    solverd = 0;
    errord = 0;
    maxfctd = 1;    //maximum number of numerical factorizations
    mnumd = 1;    //which factorization to use
    phased = 13;    //analysis
   
    pardisoinit(ptd, &mtyped, iparmd);
    nrhs = 1;
    iparmd[38] = 1;
    iparmd[34] = 1;    //0-based indexing
    
    pardiso(ptd, &maxfctd, &mnumd, &mtyped, &phased, &leng_v0d1, d, id, jd, &permd, &nrhs, iparmd, &msglvld, rhs, solution, &errord);

    return 0;
}

int interativeSolver(int N, int nrhs, double *rhs, int *ia, int *ja, double *a, int *ib, int *jb, double *b, double *solution, fdtdMesh *sys){
    // ia, ja, a are CSR form with one-based indexing and it is only upper triangular elements (symmetric)
    double mdone = -1;
    int ione = 1;
    double euclidean_norm;
    int RCI_request, itercount;
    int *ipar;
    double *dpar;
    double *tmp;
    char U = 'U';
    double *temp;
    vector<char> transa;
    int m;    // number of rows
    int k;    // number of columns
    double alpha = 1, beta = 0;
    char matdescra[6];
    int *pntrb1, *pntre1;
    int *pntrb2, *pntre2;
    int i, j;
    double *y;
    int length = 128;
    ipar = (int*)malloc((length + 2 * nrhs) * sizeof(int));
    dpar = (double*)malloc((length + 2 * nrhs) * sizeof(double));
    tmp = (double*)malloc(N*(3 + nrhs) * sizeof(double));
    y = (double*)malloc(sys->N_edge * sizeof(double));
    temp = (double*)malloc(N * sizeof(double));
    transa.push_back('T'); transa.push_back('N');
    m = N;
    k = sys->N_edge;
    matdescra[0] = 'G'; matdescra[3] = 'C';
    pntrb1 = ia;
    pntre1 = &ia[1];
    pntrb2 = ib;
    pntre2 = &ib[1];
    alpha = 1; beta = 0;

    dfgmres_init(&N, solution, rhs, &RCI_request, ipar, dpar, tmp);
    ipar[14] = 2;
    //ipar[7] = 0;
    ipar[10] = 0;
    dpar[0] = 5e-3;
    ipar[4] = 1000;   // set the maximum iteration number
    dfgmres_check(&N, solution, rhs, &RCI_request, ipar, dpar, tmp);

    while (true){
        dfgmres(&N, solution, rhs, &RCI_request, ipar, dpar, tmp);
        //cout << RCI_request << " ";
        if (RCI_request == 0){
            break;
        }
        if (RCI_request == 1){
            mkl_dcsrmv(&transa[0], &m, &k, &alpha, matdescra, a, ja, pntrb1, pntre1, tmp, &beta, y);
            mkl_dcsrmv(&transa[1], &m, &k, &alpha, matdescra, b, jb, pntrb2, pntre2, y, &beta, &tmp[N]);
            continue;
        }
        if (RCI_request == 2){
            for (j = 0; j < nrhs; j++){
                mkl_dcsrmv(&transa[0], &m, &k, &alpha, matdescra, a, ja, pntrb1, pntre1, solution, &beta, y);
                mkl_dcsrmv(&transa[1], &m, &k, &alpha, matdescra, b, jb, pntrb2, pntre2, y, &beta, temp);
            }
            cblas_daxpy(N, mdone, rhs, ione, temp, ione);
            euclidean_norm = cblas_dnrm2(N, temp, ione) / cblas_dnrm2(N, rhs, ione);
            //cout << euclidean_norm << "\n";
            if (euclidean_norm > 1.e-3)
                continue;
            else
                break;
        }
        else{
            if (dpar[6] < 1e-12) break;
            else continue;
        }
    }

    dfgmres_get(&N, solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
    //cout << itercount << "\n";
    //cout << euclidean_norm << "\n";
    free(tmp);
    free(temp);
    free(ipar);
    free(dpar);
    free(y);

    return 0;
}



int matrixMulti(int *aRowId, int *aColId, double *aval, int anum, int *bRowId, int *bColId, double *bval, int bnum, fdtdMesh *sys, int &leng, int mark){
    //the first matrix is row by row, the second matrix is column by column
    // leng is used to record the number of elements in the derived matrix

    int i = 0, j = 0;
    int flaga, flagb, k;
    int starta;
    double sum;
    /*vector<int> cRowId;
    vector<int> cColId;
    vector<double> cval;*/
    

    flaga = aRowId[0];
    flagb = bColId[0];
    starta = 0;
    sum = 0;
    while (i < anum){
        while (j < bnum && bColId[j] == flagb && i < anum && aRowId[i] == flaga){
            if (aColId[i] == bRowId[j]){
                sum += aval[i] * bval[j];
                j++;
                i++;
            }
            else if (aColId[i] < bRowId[j]){
                i++;
            }
            else if (aColId[i] > bRowId[j]){
                j++;
            }
        }
        if (sum != 0){
            /*cRowId.push_back(flaga);
            cColId.push_back(flagb);
            cval.push_back(sum);*/
            sum = 0;
            leng++;
        }
        if (i == anum){
            if (j == bnum)
                break;
            else{
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bnum){
                        while (i < anum && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == anum)    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
                continue;
            }
        }
        if (j == bnum){
            while (i < anum && aRowId[i] == flaga){
                i++;
            }
            starta = i;
            if (i == anum)    //run all of the datas
                break;
            flaga = aRowId[i];
            j = 0;
        }
        else{
            if (bColId[j] != flagb && aRowId[i] != flaga){
                flagb = bColId[j];
                i = starta;
            }
            else if (bColId[j] != flagb){
                flagb = bColId[j];
                i = starta;
            }
            else if (aRowId[i] != flaga){
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bnum){
                        while (i < anum && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == anum)    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
            }

        }
    }
    if (mark == 1){
        sys->AdRowId = (int*)malloc(leng * sizeof(int));
        sys->AdColId = (int*)malloc(leng * sizeof(int));
        sys->Adval = (double*)malloc(leng * sizeof(double));
    }
    else if (mark == 2) {
        sys->AcRowId = (int*)malloc(leng * sizeof(int));
        sys->AcColId = (int*)malloc(leng * sizeof(int));
        sys->Acval = (double*)malloc(leng * sizeof(double));
    }
    leng = 0;

    flaga = aRowId[0];
    flagb = bColId[0];
    starta = 0;
    sum = 0;
    i = 0;
    j = 0;
    int count = 1;
    sys->cindex[count] = sys->cindex[count - 1];
    while (i < anum){
        while (j < bnum && bColId[j] == flagb && i < anum && aRowId[i] == flaga){
            if (aColId[i] == bRowId[j]){
                sum += aval[i] * bval[j];
                j++;
                i++;
            }
            else if (aColId[i] < bRowId[j]){
                i++;
            }
            else if (aColId[i] > bRowId[j]){
                j++;
            }
        }
        if (sum != 0){
            if (mark == 1){
                sys->AdRowId[leng] = (flaga);
                sys->AdColId[leng] = (flagb);
                sys->Adval[leng] = (sum);
            }
            else if (mark == 2){
                sys->AcRowId[leng] = (flaga);
                sys->AcColId[leng] = (flagb);
                sys->Acval[leng] = (sum);
                if (flaga < sys->acu_cnno[count]){
                    sys->cindex[count]++;
                }
                else{
                    count++;
                    sys->cindex[count] = sys->cindex[count - 1];
                    sys->cindex[count]++;
                }
            }
            sum = 0;
            leng++;
        }
        if (i == anum){
            if (j == bnum)
                break;
            else{
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bnum){
                        while (i < anum && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == anum)    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
                continue;
            }
        }
        if (j == bnum){
            while (i < anum && aRowId[i] == flaga){
                i++;
            }
            starta = i;
            if (i == anum)    //run all of the datas
                break;
            flaga = aRowId[i];
            j = 0;
        }
        else{
            if (bColId[j] != flagb && aRowId[i] != flaga){
                flagb = bColId[j];
                i = starta;
            }
            else if (bColId[j] != flagb){
                flagb = bColId[j];
                i = starta;
            }
            else if (aRowId[i] != flaga){
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bnum){
                        while (i < anum && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == anum)    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
            }

        }
    }
    
    return 0;
}

int matrixMul(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> bRowId, vector<int> bColId, vector<double> bval, vector<int> &cRowId, vector<int> &cColId, vector<double> &cval){
    //the first matrix is row by row, the second matrix is column by column

    int i = 0, j = 0;
    int flaga, flagb, k;
    int starta;
    double sum;

    flaga = aRowId[0];
    flagb = bColId[0];
    starta = 0;
    sum = 0;
    while (i < aRowId.size()){
        while (j < bColId.size() && bColId[j] == flagb && i < aRowId.size() && aRowId[i] == flaga){
            if (aColId[i] == bRowId[j]){
                sum += aval[i] * bval[j];
                j++;
                i++;
            }
            else if (aColId[i] < bRowId[j]){
                i++;
            }
            else if (aColId[i] > bRowId[j]){
                j++;
            }
        }
        if (sum != 0){
            cRowId.push_back(flaga);
            cColId.push_back(flagb);
            cval.push_back(sum);
            sum = 0;
        }
        if (i == aRowId.size()){
            if (j == bColId.size())
                break;
            else{
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bColId.size()){
                        while (i < aRowId.size() && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == aRowId.size())    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
                continue;
            }
        }
        if (j == bColId.size()){
            while (i < aRowId.size() && aRowId[i] == flaga){
                i++;
            }
            starta = i;
            if (i == aRowId.size())    //run all of the datas
                break;
            flaga = aRowId[i];
            j = 0;
        }
        else{
            if (bColId[j] != flagb && aRowId[i] != flaga){
                flagb = bColId[j];
                i = starta;
            }
            else if (bColId[j] != flagb){
                flagb = bColId[j];
                i = starta;
            }
            else if (aRowId[i] != flaga){
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bColId.size()){
                        while (i < aRowId.size() && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == aRowId.size())    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
            }

        }
    }

    return 0;
}


int COO2CSR(vector<int> &rowId, vector<int> &ColId, vector<double> &val){
    int i;
    vector<int> rowId2;
    int count, start;

    rowId2.push_back(0);
    count = 0;
    i = 0;
    while (i < rowId.size()){
        start = rowId[i];
        while (i < rowId.size() && rowId[i] == start) {
            count++;
            i++;
        }
        rowId2.push_back(count);
    }

    rowId.clear();
    rowId = rowId2;
    return 0;
}

int COO2CSR_malloc(int *rowId, int *ColId, double *val, int totalnum, int leng, int *rowId1){    // totalnum is the total number of entries, leng is the row number
    int i;
    int *rowId2;
    int count, start;
    int k;

    rowId2 = (int*)malloc((leng + 1) * sizeof(int));
    count = 0;
    i = 0;
    k = 0;
    rowId2[k] = 0;
    k++;
    while (i < totalnum){
        start = rowId[i];
        while (i < totalnum && rowId[i] == start) {
            count++;
            i++;
        }
        rowId2[k] = (count);
        k++;
    }

    for (i = 0; i <= leng; i++){
        rowId1[i] = rowId2[i];
    }

    free(rowId2); rowId2 = NULL;
    return 0;
}

int mvMulti(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> &bRowId, vector<int> &bColId, vector<double> &bval, double *index_val, int size){
    //the same sequence in aColId and index
    double *v;
    int i;

    i = 0;
    v = (double*)calloc(size, sizeof(double));
    while (i < aColId.size()){
        v[aRowId[i]] += index_val[aColId[i]] * aval[i];
        i++;
    }
    for (i = 0; i < size; i++){
        if (abs(v[i]) > 1.e-1){
            bRowId.push_back(i);
            bColId.push_back(0);
            bval.push_back(v[i]);
        }
    }

    return 0;
}

int nodeAdd_count(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng){    // size is the size of index, total_size is the size of the vector

    int i, j;
    double *v;

    v = (double*)calloc(total_size, sizeof(double));
    for (i = 0; i < size; i++){
        for (j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            v[sys->nodeEdge[index[i]][j].first] += sys->nodeEdge[index[i]][j].second;
        }
    }
    for (i = 0; i < total_size; i++){
        if (abs(v[i]) > 1e-5) {
            num++;
        }
    }
    leng++;

    free(v);
    v = NULL;
    return 0;

}

int nodeAdd(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int mark){    // size is the size of index, total_size is the size of the vector

    int i, j;
    double *v;

    v = (double*)calloc(total_size, sizeof(double));
    for (i = 0; i < size; i++){
        for (j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            v[sys->nodeEdge[index[i]][j].first] += sys->nodeEdge[index[i]][j].second;
        }
    }
    if (mark == 1){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d1RowId[num] = (i);
                sys->v0d1ColId[num] = (leng);
                sys->v0d1val[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 2){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d2RowId[num] = (i);
                sys->v0d2ColId[num] = (leng);
                sys->v0d2val[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 3){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0cRowId[num] = (i);
                sys->v0cColId[num] = (leng);
                sys->v0cval[num] = (v[i]);
                
                num++;
            }
        }
        leng++;
    }

    free(v);
    v = NULL;
    return 0;

}

int nodeAddLarger(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int *RowId, int *ColId, double *Val){    // size is the size of index, total_size is the size of the vector

    int i, j;
    double *v;

    v = (double*)calloc(total_size, sizeof(double));
    for (i = 0; i < size; i++){
        for (j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            v[sys->nodeEdge[index[i]][j].first] += sys->nodeEdge[index[i]][j].second;
        }
    }

    for (i = 0; i < total_size; i++){
        if (abs(v[i]) > 1e-5) {
            RowId[num] = (i);
            ColId[num] = (leng);
            Val[num] = (v[i]);
            num++;
        }
    }
    leng++;

    free(v);
    v = NULL;
    return 0;

}

int nodeAddAvg_count(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng){    // Get the average V0d2 (around the conductor)
    int i, j;
    double *v;
    int inx, iny, inz;
    vector<int> rowId, colId;
    vector<double> val;

    int *nodeset = (int*)calloc(sys->N_node, sizeof(int));
    for (i = 0; i < size; i++){
        nodeset[index[i]] = 1;
    }


    v = (double*)calloc(total_size, sizeof(double));
    inz = index[0] / sys->N_node_s;
    inx = (index[0] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    iny = index[0] % (sys->N_cell_y + 1);
    if (iny == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    }
    else if (iny == sys->N_cell_y){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    }

    if (inx == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    }
    else if (inx == sys->N_cell_x){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
        colId.push_back(1);
        val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    }

    if (inz == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    }
    else if (inz == sys->N_cell_z){
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    }

    for (i = 0; i < val.size(); i++){
        v[rowId[i]] = val[i];
    }

    while (!rowId.empty()){
        rowId.pop_back();
        colId.pop_back();
        val.pop_back();
    }

    int *visited;
    stack<int> st;
    double ratio;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    int mark;
    int count;
    st.push(index[0]);
    visited[index[0]] = 1;
    while (!st.empty()){
        mark = 0;
        for (j = 0; j < sys->nodeEdge[st.top()].size(); j++){
            if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 0) && nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 1){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] % (sys->N_cell_y + 1);
                /*cout << "h\n";*/
                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] = v[rowId[i]] + ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }

                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]);
                mark = 1;

                break;
            }
            else if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 0) && nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 1){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] % (sys->N_cell_y + 1);

                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] += ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }
                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]);
                mark = 1;
                break;
            }
        }
        if (mark == 0){
            st.pop();
        }
    }

    for (i = 0; i < total_size; i++){
        if (abs(v[i]) > 1e-5) {
            num++;
        }
    }
    leng++;
    free(visited);
    visited = NULL;
    free(v);
    v = NULL;
    free(nodeset);
    nodeset = NULL;

    return 0;
}

int nodeAddAvg(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int mark){    // Get the average V0d2 (around the conductor)
    int i, j;
    double *v;
    int inx, iny, inz;
    vector<int> rowId, colId;
    vector<double> val;

    int *nodeset = (int*)calloc(sys->N_node, sizeof(int));
    for (i = 0; i < size; i++){
        nodeset[index[i]] = 1;
    }

    v = (double*)calloc(total_size, sizeof(double));
    inz = index[0] / sys->N_node_s;
    inx = (index[0] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    iny = index[0] % (sys->N_cell_y + 1);
    if (iny == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    }
    else if (iny == sys->N_cell_y){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    }

    if (inx == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    }
    else if (inx == sys->N_cell_x){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
        colId.push_back(1);
        val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    }

    if (inz == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    }
    else if (inz == sys->N_cell_z){
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    }

    for (i = 0; i < val.size(); i++){
        v[rowId[i]] = val[i];
    }
    
    while (!rowId.empty()){
        rowId.pop_back();
        colId.pop_back();
        val.pop_back();
    }

    int *visited;
    stack<int> st;
    double ratio;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    int record;
    int count;
    st.push(index[0]);
    visited[index[0]] = 1;
    while (!st.empty()){
        record = 0;
        for (j = 0; j < sys->nodeEdge[st.top()].size(); j++){
            if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] % (sys->N_cell_y + 1);
                /*cout << "h\n";*/
                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] = v[rowId[i]] + ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }

                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]);
                record = 1;

                break;
            }
            else if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] % (sys->N_cell_y + 1);

                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] += ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }
                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]);
                record = 1;
                break;
            }
        }
        if (record == 0){
            st.pop();
        }
    }

    // mark = 1 is for v0d1a, mark = 2 is for v0d2a, mark = 3 is for v0c
    if (mark == 1){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d1aRowId[num] = (i);
                sys->v0d1aColId[num] = (leng);
                sys->v0d1aval[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 2){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d2aRowId[num] = (i);
                sys->v0d2aColId[num] = (leng);
                sys->v0d2aval[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 3){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0caRowId[num] = (i);
                sys->v0caColId[num] = (leng);
                sys->v0caval[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }

    free(visited);
    visited = NULL;
    free(v);
    v = NULL;
    free(nodeset);
    nodeset = NULL;

    return 0;
}

int nodeAddAvgLarger(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int* RowId, int *ColId, double *Val){    // Get the average V0d2 (around the conductor)
    int i, j;
    double *v;
    int inx, iny, inz;
    vector<int> rowId, colId;
    vector<double> val;

    int *nodeset = (int*)calloc(sys->N_node, sizeof(int));
    for (i = 0; i < size; i++){
        nodeset[index[i]] = 1;
    }

    v = (double*)calloc(total_size, sizeof(double));
    inz = index[0] / sys->N_node_s;
    inx = (index[0] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    iny = index[0] % (sys->N_cell_y + 1);
    if (iny == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    }
    else if (iny == sys->N_cell_y){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    }

    if (inx == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    }
    else if (inx == sys->N_cell_x){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
        colId.push_back(1);
        val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    }

    if (inz == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    }
    else if (inz == sys->N_cell_z){
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    }

    for (i = 0; i < val.size(); i++){
        v[rowId[i]] = val[i];
    }

    while (!rowId.empty()){
        rowId.pop_back();
        colId.pop_back();
        val.pop_back();
    }

    int *visited;
    stack<int> st;
    double ratio;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    int record;
    int count;
    st.push(index[0]);
    visited[index[0]] = 1;
    while (!st.empty()){
        record = 0;
        for (j = 0; j < sys->nodeEdge[st.top()].size(); j++){
            if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] % (sys->N_cell_y + 1);
                /*cout << "h\n";*/
                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] = v[rowId[i]] + ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }

                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]);
                record = 1;

                break;
            }
            else if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] % (sys->N_cell_y + 1);

                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] += ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }
                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]);
                record = 1;
                break;
            }
        }
        if (record == 0){
            st.pop();
        }
    }


    for (i = 0; i < total_size; i++){
        if (abs(v[i]) > 1e-5) {
            RowId[num] = (i);
            ColId[num] = (leng);
            Val[num] = (v[i]);
            num++;
        }
    }
    leng++;



    free(visited);
    visited = NULL;
    free(v);
    v = NULL;
    free(nodeset);
    nodeset = NULL;

    return 0;
}


int solveV0dSystem(fdtdMesh *sys, double *dRhs, complex<double> *y0d, double *v0d2epsv0d2, double *solution_d2, int leng_v0d1, int leng_v0d2){
    double *mV0da1J = (double*)calloc(leng_v0d1, sizeof(double));
    double *mV0da2J = (double*)calloc(leng_v0d2, sizeof(double));
    double *y0d2 = (double*)malloc(leng_v0d2*sizeof(double));
    char transa;
    double alpha = 1, beta = 0;
    char matdescra[6];
    ofstream outfile;
    matdescra[0] = 'G'; matdescra[3] = 'C';    // general matrix multi, 0-based indexing

    
    transa = 'N';
    int k = sys->N_edge;
    int m = leng_v0d1;
    mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1avalo[0], &sys->v0d1aRowId[0], &sys->v0d1aColId[0], &sys->v0d1aColId[1], dRhs, &beta, mV0da1J);
    m = leng_v0d2;
    mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2avalo[0], &sys->v0d2aRowId[0], &sys->v0d2aColId[0], &sys->v0d2aColId[1], dRhs, &beta, mV0da2J);
    
    
    /* A\b1 */
    double *d = &(sys->Adval[0]);
    int *id = &(sys->AdRowId1[0]);
    int *jd = &(sys->AdColId[0]);
    double *Ainvb1 = (double*)malloc(leng_v0d1*sizeof(double));

    void *ptd[64];
    int mtyped;
    int iparmd[64];
    double dparmd[64];
    int maxfctd, mnumd, phased, errord, solverd;
    int num_processd;   //number of processors
    int v0csin;
    int permd;
    int nrhs = 1;
    int msglvld = 0;    //print statistical information

    mtyped = 11;    // real and not symmetric
    solverd = 0;
    errord = 0;
    maxfctd = 1;    //maximum number of numerical factorizations
    mnumd = 1;    //which factorization to use
    phased = 13;    //analysis

    pardisoinit(ptd, &mtyped, iparmd);
    nrhs = 1;
    iparmd[38] = 1;
    iparmd[34] = 1;    //0-based indexing
    pardiso(ptd, &maxfctd, &mnumd, &mtyped, &phased, &leng_v0d1, d, id, jd, &permd, &nrhs, iparmd, &msglvld, mV0da1J, Ainvb1, &errord);
    
    /* b2-CA\b1 */

    transa = 'T';
    k = sys->N_edge;
    m = leng_v0d1;
    double *temp = (double*)malloc(sys->N_edge*sizeof(double));
    mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1val[0], &sys->v0d1RowId[0], &sys->v0d1ColId[0], &sys->v0d1ColId[1], Ainvb1, &beta, temp);
    transa = 'N';
    m = leng_v0d2;
    double *cAinvb1 = (double*)malloc(leng_v0d2*sizeof(double));
    mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2aval[0], &sys->v0d2aRowId[0], &sys->v0d2aColId[0], &sys->v0d2aColId[1], temp, &beta, cAinvb1);
    double mdone = -1;
    int ione = 1;
    for (int id = 0; id < leng_v0d2; id++){
        mV0da2J[id] = mV0da2J[id] - cAinvb1[id];
    }
    
    /*cblas_daxpy(nrhs, mdone, cAinvb1, ione, mV0da2J, ione);*/

    free(temp);
    temp = NULL;
    /* y0d2 */
    y0d2 = (double*)calloc(leng_v0d2, sizeof(double));
    for (int id = 0; id < leng_v0d2; id++){
        for (int jd = 0; jd < leng_v0d2; jd++){
            y0d2[id] = y0d2[id] + v0d2epsv0d2[id + jd*leng_v0d2] * mV0da2J[jd];
        }
    }

    /* y0d1 */
    double *y0d1 = (double*)calloc(leng_v0d1, sizeof(double));
    for (int id = 0; id < leng_v0d2; id++){
        for (int jd = 0; jd < leng_v0d1; jd++){
            y0d1[jd] = y0d1[jd] - solution_d2[jd + id*leng_v0d1] * y0d2[id];
        }
    }
    mdone = 1;
    cblas_daxpy(leng_v0d1, mdone, Ainvb1, ione, y0d1, ione);
    
    /* final y0d */
    double *y0dd = (double*)malloc(sys->N_edge * sizeof(double));
    temp = (double*)malloc(sys->N_edge * sizeof(double));
    transa = 'T';
    m = leng_v0d1;
    k = sys->N_edge;
    mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1valo[0], &sys->v0d1RowId[0], &sys->v0d1ColId[0], &sys->v0d1ColId[1], y0d1, &beta, temp);
    /*outfile.open("y0d1.txt", std::ofstream::out | std::ofstream::trunc);
    for (int i = 0; i < leng_v0d1; i++){
        outfile << y0d1[i] << endl;
    }
    outfile.close();*/

    m = leng_v0d2;
    mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2valo[0], &sys->v0d2RowId[0], &sys->v0d2ColId[0], &sys->v0d2ColId[1], y0d2, &beta, y0dd);
    mdone = 1;
    cblas_daxpy(sys->N_edge, mdone, temp, ione, y0dd, ione);
    
    complex<double> ima(0.0, 1.0);
    for (int id = 0; id < sys->N_edge; id++){
        y0d[id] = -y0dd[id] / (2 * PI*sys->freqStart * sys->freqUnit) * ima;
    }
    
    free(y0dd); y0dd = NULL;
    free(temp); temp = NULL;

    return 1;

}

int merge_v0d1(fdtdMesh *sys, double block1_x, double block1_y, double block2_x, double block2_y, int &v0d1num, int &leng_v0d1, int &v0d1anum, int &leng_v0d1a){
    int *visited;

    double startx, starty;    // the start coordinates of each block
    vector<int> st;    // dfs stack
    vector<int> ind;
    int indsize;
    int indx, indy;
    int mark;
    int status;
    int indnum;

    int *markLayerNode = (int*)calloc(sys->N_node_s, sizeof(int));
    for (int i = 0; i < sys->numPorts; i++){
        for (int j = 0; j < sys->cdtNumNode[sys->portCoor[i].portCnd-1]; j++){
            markLayerNode[sys->conductor[sys->portCoor[i].portCnd - 1].node[j] % (sys->N_node_s)] = 1;
        }
    }

    leng_v0d1 = 0;
    leng_v0d1a = 0;
    v0d1num = 0;
    v0d1anum = 0;
    /* first asign a larger number of storage, don't need to calculate the entries twice */
    int *v0d1RowId = (int*)malloc(sys->N_edge * 2 * sizeof(int));
    int *v0d1ColId = (int*)malloc(sys->N_edge * 2 * sizeof(int));
    double *v0d1val = (double*)malloc(sys->N_edge * 2 * sizeof(double));
    int *v0d1aRowId = (int*)malloc(sys->N_edge * 2 * sizeof(int));
    int *v0d1aColId = (int*)malloc(sys->N_edge * 2 * sizeof(int));
    double *v0d1aval = (double*)malloc(sys->N_edge * 2 * sizeof(double));

    for (int iz = 1; iz < sys->nz - 1; iz++){    // merge on each layer, not in the conductor
        visited = (int*)calloc(sys->nx * sys->ny, sizeof(int));
        for (int ix = 0; ix < sys->nx; ix++){
            for (int iy = 0; iy < sys->ny; iy++){
                if (visited[ix * (sys->N_cell_y + 1) + iy] == 0 && sys->markNode[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] == 0 && markLayerNode[ix * (sys->N_cell_y + 1) + iy] == 0){ //&& sys->exciteCdtLayer[iz] == 0) {    // this point is not visited and it is outside the conductor, not in the projection of the excited conductor
                    if (!ind.empty())
                        ind.clear();
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    ind.push_back(iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy);
                    st.push_back(ix * (sys->N_cell_y + 1) + iy);
                    visited[ix * (sys->N_cell_y + 1) + iy] = 1;
                    while (!st.empty()){
                        mark = 0;
                        indx = (st.back()) / (sys->N_cell_y + 1);
                        indy = st.back() % (sys->N_cell_y + 1);
                        if (indx != sys->nx-1){    // it must have a right x edge, thus right x node
                            if (sys->markNode[iz * sys->N_node_s + st.back() + sys->N_cell_y + 1] == 0 && visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block1_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block1_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markNode[iz * sys->N_node_s + st.back() - sys->N_cell_y - 1] == 0 && visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block1_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block1_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markNode[iz * sys->N_node_s + st.back() + 1] == 0 && visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 0){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block1_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block1_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markNode[iz * sys->N_node_s + st.back() - 1] == 0 && visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block1_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block1_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[(indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (mark == 0){
                            st.pop_back();
                        }
                    }
                    indnum = ind.size();
                    
                    if (indnum != 0){
                        status = nodeAddLarger(&ind[0], indnum, sys->N_edge, sys, v0d1num, leng_v0d1, v0d1RowId, v0d1ColId, v0d1val);    // used to find v0d1num and leng_v0d1
                        if (status != 0)
                            return status;

                        status = nodeAddAvgLarger(&ind[0], indnum, sys->N_edge, sys, v0d1anum, leng_v0d1a, v0d1aRowId, v0d1aColId, v0d1aval);    // used to find v0d1anum and leng_v0d1a
                        if (status != 0)
                            return status;
                    }

                }

                else if (visited[ix * (sys->N_cell_y + 1) + iy] == 0 && sys->markNode[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] == 0 && markLayerNode[ix * (sys->N_cell_y + 1) + iy] == 1){//&& sys->exciteCdtLayer[iz] == 1){    // this point is not visited and it is outside the conductor, in the projection of the excited conductor
                    while (!ind.empty())
                        ind.clear();
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    ind.push_back(iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy);
                    st.push_back(ix * (sys->N_cell_y + 1) + iy);
                    visited[ix * (sys->N_cell_y + 1) + iy] = 1;
                    while (!st.empty()){
                        mark = 0;
                        indx = (st.back()) / (sys->N_cell_y + 1);
                        indy = st.back() % (sys->N_cell_y + 1);
                        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                            if (sys->markNode[iz * sys->N_node_s + st.back() + sys->N_cell_y + 1] == 0 && visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 1){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markNode[iz * sys->N_node_s + st.back() - sys->N_cell_y - 1] == 0 && visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 1){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markNode[iz * sys->N_node_s + st.back() + 1] == 0 && visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 1){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markNode[iz * sys->N_node_s + st.back() - 1] == 0 && visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] == 1){    // this node is in dielectric and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[(indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (mark == 0){
                            st.pop_back();
                        }
                    }
                    indnum = ind.size();
                    if (indnum != 0){
                        status = nodeAddLarger(&ind[0], indnum, sys->N_edge, sys, v0d1num, leng_v0d1, v0d1RowId, v0d1ColId, v0d1val);    // generate v0d1RowId, v0d1ColId, v0d1val with a larger storage
                        if (status != 0)
                            return status;

                        status = nodeAddAvgLarger(&ind[0], indnum, sys->N_edge, sys, v0d1anum, leng_v0d1a, v0d1aRowId, v0d1aColId, v0d1aval);    // generate v0d1aRowId, v0d1aColId, v0d1aval with a larger storage
                        if (status != 0)
                            return status;
                    }
                }

            }
        }
        free(visited); visited = NULL;
    }

    
    sys->v0d1RowId = (int*)malloc(v0d1num * sizeof(int));
    sys->v0d1ColId = (int*)malloc(v0d1num * sizeof(int));
    sys->v0d1val = (double*)malloc(v0d1num * sizeof(double));
    sys->v0d1aRowId = (int*)malloc(v0d1anum * sizeof(int));
    sys->v0d1aColId = (int*)malloc(v0d1anum * sizeof(int));
    sys->v0d1aval = (double*)malloc(v0d1anum * sizeof(double));
    

    for (int i = 0; i < v0d1num; i++){
        sys->v0d1RowId[i] = v0d1RowId[i];
        sys->v0d1ColId[i] = v0d1ColId[i];
        sys->v0d1val[i] = v0d1val[i];
    }
    for (int i = 0; i < v0d1anum; i++){
        sys->v0d1aRowId[i] = v0d1aRowId[i];
        sys->v0d1aColId[i] = v0d1aColId[i];
        sys->v0d1aval[i] = v0d1aval[i];
    }

    free(v0d1RowId); v0d1RowId = NULL;
    free(v0d1ColId); v0d1ColId = NULL;
    free(v0d1val); v0d1val = NULL;
    free(v0d1aRowId); v0d1aRowId = NULL;
    free(v0d1aColId); v0d1aColId = NULL;
    free(v0d1aval); v0d1aval = NULL;

    ind.clear();

    return 1;
}

int setsideLen(int node, int iz, double sideLen, int *markLayerNode, fdtdMesh *sys){
    queue<int> q;
    int *visited;
    int indx, indy;
    int mark;
    double startx, starty;

    if (sideLen == 0)
        return 0;

    visited = (int*)calloc(sys->N_node_s, sizeof(int));
    q.push(node);
    visited[node] = 1;
    startx = sys->xn[node / (sys->N_cell_y + 1)];
    starty = sys->yn[node % (sys->N_cell_y + 1)];

    while (!q.empty()){
        indx = q.front() / (sys->N_cell_y + 1);
        indy = q.front() % (sys->N_cell_y + 1);

        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node

            if (sys->markNode[iz * sys->N_node_s + q.front() + sys->N_cell_y + 1] == 0 && visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx + 1] - startx), 2) + pow((sys->yn[indy] - starty), 2)) <= sideLen){    // this node is within the block area
                    q.push((indx + 1)*(sys->N_cell_y + 1) + indy);
                    visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                    markLayerNode[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                }
            }
        }

        if (indx != 0){    // it must have a left x edge, thus left x node
            if (sys->markNode[iz * sys->N_node_s + q.front() - sys->N_cell_y - 1] == 0 && visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx - 1] - startx), 2) + pow((sys->yn[indy] - starty), 2)) <= sideLen){    // this node is within the block area
                    q.push((indx - 1)*(sys->N_cell_y + 1) + indy);
                    visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                    markLayerNode[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                }
            }
        }

        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
            if (sys->markNode[iz * sys->N_node_s + q.front() + 1] == 0 && visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 1){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx] - startx), 2) + pow((sys->yn[indy + 1] - starty), 2)) <= sideLen){    // this node is within the block area
                    q.push((indx)*(sys->N_cell_y + 1) + indy + 1);
                    visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                    markLayerNode[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                }
            }
        }

        if (indy != 0){    // it must have a closer y edge, thus closer y node
            if (sys->markNode[iz * sys->N_node_s + q.front() - 1] == 0 && visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] == 1){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx] - startx), 2) + pow((sys->yn[indy - 1] - starty), 2)) <= sideLen){    // this node is within the block area
                    q.push((indx)*(sys->N_cell_y + 1) + indy - 1);
                    visited[(indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                    markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] = 1;
                }
            }
        }

        q.pop();
    }

    free(visited); visited = NULL;

    return 0;
}

int merge_v0c(fdtdMesh *sys, double block_x, double block_y, double block2_x, double block2_y, int &v0cnum, int &leng_v0c, int &v0canum, int &leng_v0ca){
    int *visited;

    double startx, starty;    // the start coordinates of each block
    vector<int> st;    // dfs stack
    vector<int> ind;
    int indsize;
    int indx, indy;
    int mark;
    int status;
    int markcond;
    int count = 0;

    leng_v0c = 0;
    leng_v0ca = 0;
    v0cnum = 0;
    v0canum = 0;
    int indnum;
    int ix, iy, iz;
    int n;

    for (int ic = 0; ic < sys->numCdt; ic++){
        cout << ic << endl;
        if (sys->conductor[ic].markPort <= 0){
            markcond = sys->markNode[sys->conductor[ic].node[0]];
            visited = (int*)calloc(sys->N_node, sizeof(int));
            n = sys->cdtNumNode[ic] - 1;
            if (sys->conductor[ic].markPort == -1)
                n = sys->cdtNumNode[ic];
            for (int jc = 0; jc < n; jc++){
                if (visited[sys->conductor[ic].node[jc]] == 0 && sys->conductor[ic].node[jc] >= sys->N_node_s && sys->conductor[ic].node[jc] < sys->N_node - sys->N_node_s){

                    if (!ind.empty())
                        ind.clear();
                    iz = sys->conductor[ic].node[jc] / (sys->N_node_s);
                    ix = (sys->conductor[ic].node[jc] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iy = sys->conductor[ic].node[jc] % (sys->N_cell_y + 1);
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    ind.push_back(sys->conductor[ic].node[jc]);
                    st.push_back(ix * (sys->N_cell_y + 1) + iy);
                    visited[sys->conductor[ic].node[jc]] = 1;

                    while (!st.empty()){

                        mark = 0;
                        indx = (st.back()) / (sys->N_cell_y + 1);
                        indy = st.back() % (sys->N_cell_y + 1);

                        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * indx + indy] == markcond && visited[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy] == markcond && visited[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy - 1] == markcond && visited[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy - 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (mark == 0){

                            st.pop_back();

                        }
                    }

                    indnum = ind.size();


                    if (indnum != 0){

                        status = nodeAdd_count(&ind[0], indnum, sys->N_edge, sys, v0cnum, leng_v0c);    // used to find v0d1num and leng_v0d1
                        if (status != 0)
                            return status;

                        status = nodeAddAvg_count(&ind[0], indnum, sys->N_edge, sys, v0canum, leng_v0ca);    // used to find v0d1anum and leng_v0d1a
                        if (status != 0)
                            return status;

                    }
                }
            }

            sys->acu_cnno[count + 1] = leng_v0c;
            count++;
            free(visited); visited = NULL;
        }
        else{
            //for (int j = 0; j < sys->cdtNumNode[ic] - 1; j++){    // there is one node needed to be removed
            //    
            //    for (int k = 0; k < sys->nodeEdge[sys->conductor[ic].node[j]].size(); k++){
            //        v0cnum++;
            //    }

            //    leng_v0c++;

            //    iz = sys->conductor[ic].node[j] / sys->N_node_s;
            //    ix = (sys->conductor[ic].node[j] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
            //    iy = sys->conductor[ic].node[j] % (sys->N_cell_y + 1);
            //    if (iz != 0 && iz != sys->N_cell_z){
            //        v0canum++;
            //    }
            //    else if (iz == sys->N_cell_z){
            //        v0canum++;
            //    }

            //    if (iy != 0 && iy != sys->N_cell_y){
            //        v0canum++;
            //        v0canum++;
            //    }
            //    else if (iy == sys->N_cell_y){
            //        v0canum++;
            //    }
            //    else{
            //        v0canum++;
            //    }

            //    if (ix != 0 && ix != sys->N_cell_x){
            //        v0canum++;
            //        v0canum++;
            //    }
            //    else if (ix == sys->N_cell_x){
            //        v0canum++;
            //    }
            //    else{
            //        v0canum++;
            //    }

            //    if (iz != sys->N_cell_z && iz != 0){
            //        v0canum++;
            //    }
            //    else if (iz == 0){
            //        v0canum++;
            //    }

            //    leng_v0ca++;
            //}

            markcond = sys->markNode[sys->conductor[ic].node[0]];
            visited = (int*)calloc(sys->N_node, sizeof(int));
            n = sys->cdtNumNode[ic] - 1;
            for (int jc = 0; jc < n; jc++){
                if (visited[sys->conductor[ic].node[jc]] == 0 && sys->conductor[ic].node[jc] >= sys->N_node_s && sys->conductor[ic].node[jc] < sys->N_node - sys->N_node_s){

                    if (!ind.empty())
                        ind.clear();
                    iz = sys->conductor[ic].node[jc] / (sys->N_node_s);
                    ix = (sys->conductor[ic].node[jc] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iy = sys->conductor[ic].node[jc] % (sys->N_cell_y + 1);
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    ind.push_back(sys->conductor[ic].node[jc]);
                    st.push_back(ix * (sys->N_cell_y + 1) + iy);
                    visited[sys->conductor[ic].node[jc]] = 1;

                    while (!st.empty()){

                        mark = 0;
                        indx = (st.back()) / (sys->N_cell_y + 1);
                        indy = st.back() % (sys->N_cell_y + 1);

                        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * indx + indy] == markcond && visited[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy] == markcond && visited[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block2_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy - 1] == markcond && visited[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy - 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block2_y){    // this node is within the block area

                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    mark = 1;

                                    continue;
                                }
                            }
                        }
                        if (mark == 0){

                            st.pop_back();

                        }
                    }

                    indnum = ind.size();


                    if (indnum != 0){

                        status = nodeAdd_count(&ind[0], indnum, sys->N_edge, sys, v0cnum, leng_v0c);    // used to find v0d1num and leng_v0d1
                        if (status != 0)
                            return status;

                        status = nodeAddAvg_count(&ind[0], indnum, sys->N_edge, sys, v0canum, leng_v0ca);    // used to find v0d1anum and leng_v0d1a
                        if (status != 0)
                            return status;

                    }
                }
            }

            sys->acu_cnno[count + 1] = leng_v0c;
            count++;
            free(visited); visited = NULL;

        }
    }


    sys->v0cRowId = (int*)malloc(v0cnum * sizeof(int));
    sys->v0cColId = (int*)malloc(v0cnum * sizeof(int));
    sys->v0cval = (double*)malloc(v0cnum * sizeof(double));
    sys->v0caRowId = (int*)malloc(v0canum * sizeof(int));
    sys->v0caColId = (int*)malloc(v0canum * sizeof(int));
    sys->v0caval = (double*)malloc(v0canum * sizeof(double));
    leng_v0c = 0;
    leng_v0ca = 0;
    v0cnum = 0;
    v0canum = 0;

    for (int ic = 0; ic < sys->numCdt; ic++){
        cout << ic << endl;
        if (sys->conductor[ic].markPort <= 0){
            markcond = sys->markNode[sys->conductor[ic].node[0]];
            visited = (int*)calloc(sys->N_node, sizeof(int));
            n = sys->cdtNumNode[ic] - 1;
            if (sys->conductor[ic].markPort == -1)
                n = sys->cdtNumNode[ic];
            for (int jc = 0; jc < n; jc++){
                if (visited[sys->conductor[ic].node[jc]] == 0 && sys->conductor[ic].node[jc] >= sys->N_node_s && sys->conductor[ic].node[jc] < sys->N_node - sys->N_node_s){
                    if (!ind.empty())
                        ind.clear();
                    iz = sys->conductor[ic].node[jc] / (sys->N_node_s);
                    ix = (sys->conductor[ic].node[jc] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iy = sys->conductor[ic].node[jc] % (sys->N_cell_y + 1);
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    ind.push_back(sys->conductor[ic].node[jc]);
                    st.push_back(ix * (sys->N_cell_y + 1) + iy);
                    visited[sys->conductor[ic].node[jc]] = 1;
                    while (!st.empty()){
                        mark = 0;
                        indx = (st.back()) / (sys->N_cell_y + 1);
                        indy = st.back() % (sys->N_cell_y + 1);
                        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * indx + indy] == markcond && visited[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy] == markcond && visited[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy - 1] == markcond && visited[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy - 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (mark == 0){
                            st.pop_back();
                        }
                    }
                    indnum = ind.size();

                    if (indnum != 0){
                        status = nodeAdd(&ind[0], indnum, sys->N_edge, sys, v0cnum, leng_v0c, 3);    // used to generate v0c
                        if (status != 0)
                            return status;

                        status = nodeAddAvg(&ind[0], indnum, sys->N_edge, sys, v0canum, leng_v0ca, 3);    // used to generate v0dc
                        if (status != 0)
                            return status;
                    }
                }
            }
            free(visited); visited = NULL;
        }
        else{
            //for (int j = 0; j < sys->cdtNumNode[ic] - 1; j++){    // there is one node needed to be removed
            //
            //    for (int k = 0; k < sys->nodeEdge[sys->conductor[ic].node[j]].size(); k++){
            //        sys->v0cColId[v0cnum] = (leng_v0c);
            //        sys->v0cRowId[v0cnum] = (sys->nodeEdge[sys->conductor[ic].node[j]][k].first);
            //        sys->v0cval[v0cnum] = (sys->nodeEdge[sys->conductor[ic].node[j]][k].second);
            //        v0cnum++;
            //    }

            //    leng_v0c++;

            //    iz = sys->conductor[ic].node[j] / sys->N_node_s;
            //    ix = (sys->conductor[ic].node[j] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
            //    iy = sys->conductor[ic].node[j] % (sys->N_cell_y + 1);
            //    if (iz != 0 && iz != sys->N_cell_z){
            //        sys->v0caRowId[v0canum] = ((iz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix*(sys->N_cell_y + 1) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (-2 / (sys->zn[iz + 1] - sys->zn[iz - 1]));
            //        v0canum++;
            //        
            //    }
            //    else if (iz == sys->N_cell_z){
            //        sys->v0caRowId[v0canum] = ((iz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix*(sys->N_cell_y + 1) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (-1 / (sys->zn[iz] - sys->zn[iz - 1]));
            //        v0canum++;
            //        
            //    }

            //    if (iy != 0 && iy != sys->N_cell_y){
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + ix*sys->N_cell_y + iy - 1);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (-2 / (sys->yn[iy + 1] - sys->yn[iy - 1]));
            //        v0canum++;
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + ix*sys->N_cell_y + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (2 / (sys->yn[iy + 1] - sys->yn[iy - 1]));
            //        v0canum++;
            //    }
            //    else if (iy == sys->N_cell_y){
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + ix*sys->N_cell_y + iy - 1);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (-1 / (sys->yn[iy] - sys->yn[iy - 1]));
            //        v0canum++;
            //    }
            //    else{
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + ix*(sys->N_cell_y) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (1 / (sys->yn[iy + 1] - sys->yn[iy]));
            //        v0canum++;
            //    }

            //    if (ix != 0 && ix != sys->N_cell_x){
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(ix - 1) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (-2 / (sys->xn[ix + 1] - sys->xn[ix - 1]));
            //        v0canum++;
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(ix)+iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (2 / (sys->xn[ix + 1] - sys->xn[ix - 1]));
            //        v0canum++;
            //    }
            //    else if (ix == sys->N_cell_x){
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(ix - 1) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (-1 / (sys->xn[ix] - sys->xn[ix - 1]));
            //        v0canum++;
            //    }
            //    else{
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (1 / (sys->xn[ix + 1] - sys->xn[ix]));
            //        v0canum++;
            //    }

            //    if (iz != sys->N_cell_z && iz != 0){
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix*(sys->N_cell_y + 1) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (2 / (sys->zn[iz + 1] - sys->zn[iz - 1]));
            //        v0canum++;
            //    }
            //    else if (iz == 0){
            //        sys->v0caRowId[v0canum] = (iz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + ix*(sys->N_cell_y + 1) + iy);
            //        sys->v0caColId[v0canum] = (leng_v0ca);
            //        sys->v0caval[v0canum] = (1 / (sys->zn[iz + 1] - sys->zn[iz]));
            //        v0canum++;
            //    }

            //    leng_v0ca++;

            //}

            markcond = sys->markNode[sys->conductor[ic].node[0]];
            visited = (int*)calloc(sys->N_node, sizeof(int));
            n = sys->cdtNumNode[ic] - 1;

            for (int jc = 0; jc < n; jc++){
                if (visited[sys->conductor[ic].node[jc]] == 0 && sys->conductor[ic].node[jc] >= sys->N_node_s && sys->conductor[ic].node[jc] < sys->N_node - sys->N_node_s){
                    if (!ind.empty())
                        ind.clear();
                    iz = sys->conductor[ic].node[jc] / (sys->N_node_s);
                    ix = (sys->conductor[ic].node[jc] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iy = sys->conductor[ic].node[jc] % (sys->N_cell_y + 1);
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    ind.push_back(sys->conductor[ic].node[jc]);
                    st.push_back(ix * (sys->N_cell_y + 1) + iy);
                    visited[sys->conductor[ic].node[jc]] = 1;
                    while (!st.empty()){
                        mark = 0;
                        indx = (st.back()) / (sys->N_cell_y + 1);
                        indy = st.back() % (sys->N_cell_y + 1);
                        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * indx + indy] == markcond && visited[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy] == markcond && visited[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy - 1] == markcond && visited[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy - 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block2_y){    // this node is within the block area
                                    ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    mark = 1;
                                    continue;
                                }
                            }
                        }
                        if (mark == 0){
                            st.pop_back();
                        }
                    }
                    indnum = ind.size();

                    if (indnum != 0){
                        status = nodeAdd(&ind[0], indnum, sys->N_edge, sys, v0cnum, leng_v0c, 3);    // used to generate v0c
                        if (status != 0)
                            return status;

                        status = nodeAddAvg(&ind[0], indnum, sys->N_edge, sys, v0canum, leng_v0ca, 3);    // used to generate v0dc
                        if (status != 0)
                            return status;
                    }
                }
            }
            free(visited); visited = NULL;
        }
    }
    

    ind.clear();

    return 1;
}