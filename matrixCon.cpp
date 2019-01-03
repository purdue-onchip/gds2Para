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
    vector<int> rowId, rowId1;
    vector<int> colId, colId1;
    vector<double> val, val1;
    vector<int> temp2;

    //cout << "Number of conductors: " << sys->numCdt << endl;
    cout << "Number of edges: " << sys->N_edge << endl;
    cout << "Number of conductors: " << sys->numCdt << endl;
    cout << "Number of ports: " << sys->numPorts << endl;

    /* Merged conductor node */
    /*double dcx = 5e-7, dcy = 5e-7, dcz = 5e-7;
    int ncx = ceil(sys->xlim2 / dcx), ncy = ceil(sys->ylim2 / dcy), ncz = ceil(sys->zlim2 / dcz);
    int mcnno;    // merged node number
    unordered_map<int, vector<int> > groupc;
    int leng_v0c = 0;
    int *cindex = (int*)calloc((sys->numCdt + 1), (sizeof(int)));
    int *acu_cnno = (int*)calloc((sys->numCdt), sizeof(int));
    cindex[0] = -1;

    for (i = 0; i < sys->numCdt; i++){
        for (j = 0; j < sys->cdtNumNode[i] - 1; j++){    // No port, need one node removed
            mcnno = 0;
            if (sys->nodepos[sys->conductor[i].node[j] * 3 + 2] == 0){
                mcnno = mcnno + 0;
            }
            else{
                mcnno = mcnno + (ncx * ncy) * (ceil((sys->nodepos[sys->conductor[i].node[j] * 3 + 2] / dcz)) - 1);
            }
            if (sys->nodepos[sys->conductor[i].node[j] * 3] == 0){
                mcnno = mcnno + 0;
            }
            else{
                mcnno = mcnno + ncy * (ceil((sys->nodepos[sys->conductor[i].node[j] * 3] / dcx)) - 1);
            }
            if (sys->nodepos[sys->conductor[i].node[j] * 3 + 1] == 0){
                mcnno = mcnno + 0;
            }
            else{
                mcnno = mcnno + ceil(sys->nodepos[sys->conductor[i].node[j] * 3 + 1] / dcy) - 1;
            }
            groupc[mcnno].push_back(sys->conductor[i].node[j]);
        }
        for (auto gc : groupc){
            status = nodeAdd(rowId, colId, val, gc.second, sys->N_edge, sys);
            if (status != 0)
                return status;
            for (j = 0; j < rowId.size(); j++){
                sys->v0cRowId.push_back(rowId[j]);
                sys->v0cColId.push_back(leng_v0c);
                sys->v0cval.push_back(val[j]);
                cindex[i + 1]++;
            }
            leng_v0c++;
        }
        acu_cnno[i] = leng_v0c;
        groupc.clear();
        cindex[i + 1] = cindex[i] + cindex[i + 1];
    }
    sys->v0cvalo = sys->v0cval;    // v0cvalo is the v0c values without D_sig
    for (i = 0; i < sys->v0cColId.size(); i++){
        sys->v0cval[i] = sys->v0cval[i] * sqrt(sys->sig[sys->v0cRowId[i]]);       // Compute the sparse form of D_sig*V0c
    }

    /* Construct V0c with row id, col id and its val */
    int leng_v0c = 0;
    int leng_v0ca = 0;
    int inx, iny, inz;
    int numPortCdt = 0;

    for (i = 0; i < sys->numCdt; i++){
        n = sys->cdtNumNode[i];
        sys->conductor[i].markPort = 0;
        for (j = 0; j < n; j++){    // No port, need one node removed
            if (sys->portNno.find(sys->conductor[i].node[j]) != sys->portNno.end()){
                sys->conductor[i].markPort = 1;
                numPortCdt++;
                break;
            }
        }
    }
    int *cindex = (int*)calloc(numPortCdt + 1, (sizeof(int)));
    int *acu_cnno = (int*)calloc((numPortCdt + 1), sizeof(int));
    cindex[0] = -1;    // the last index in the sparse form for each conductor in V0c
    acu_cnno[0] = 0;
    count = 0;
    int *map = (int*)calloc(sys->N_node, (sizeof(int)));

    for (i = 0; i < sys->numCdt; i++){
        if (sys->conductor[i].markPort == 1){
            acu_cnno[count + 1] = acu_cnno[count];
            n = sys->cdtNumNode[i];
            for (j = 0; j < n; j++){
                if (sys->portNno.find(sys->conductor[i].node[j]) != sys->portNno.end()){
                    continue;
                }
                for (k = 0; k < sys->nodeEdge[sys->conductor[i].node[j]].size(); k++){
                    sys->v0cColId.push_back(leng_v0c);
                    sys->v0cRowId.push_back(sys->nodeEdge[sys->conductor[i].node[j]][k].first);
                    sys->v0cval.push_back(sys->nodeEdge[sys->conductor[i].node[j]][k].second);
                    cindex[count + 1]++;
                }
                map[sys->conductor[i].node[j]] = leng_v0c + 1;

                acu_cnno[count + 1]++;
                leng_v0c++;

                inz = sys->conductor[i].node[j] / sys->N_node_s;
                inx = (sys->conductor[i].node[j] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->conductor[i].node[j] % (sys->N_cell_y + 1);
                if (inz != 0 && inz != sys->N_cell_z){
                    sys->v0caRowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (sys->markEdge[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] - sys->N_node_s) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] - sys->N_node_s);
                        sys->Acval.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                    }
                }
                else if (inz == sys->N_cell_z){
                    sys->v0caRowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if (sys->markEdge[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] - sys->N_node_s) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] - sys->N_node_s);
                        sys->Acval.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                    } 
                }

                if (iny != 0 && iny != sys->N_cell_y){
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                }
                else if (iny == sys->N_cell_y){
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                }
                else{
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                }

                if (inx != 0 && inx != sys->N_cell_x){
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                }
                else if (inx == sys->N_cell_x){
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                }
                else{
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                }

                if (inz != sys->N_cell_z && inz != 0){
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                }
                else if (inz == 0){
                    sys->v0caRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    sys->v0caColId.push_back(leng_v0ca);
                    sys->v0caval.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                }

                leng_v0ca++;

                if (inx != 0 && inx != sys->N_cell_x){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] - sys->N_cell_y - 1) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] - sys->N_cell_y - 1);
                        sys->Acval.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]);
                    }
                    if (iny != 0 && iny != sys->N_cell_y){
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->portNno.find(sys->conductor[i].node[j] - 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] - 1);
                            sys->Acval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                        }
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz])* sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] + 1);
                            sys->Acval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                        }
                    }
                    else if (iny == sys->N_cell_y){
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->portNno.find(sys->conductor[i].node[j] - 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] - 1);
                            sys->Acval.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                        }
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                    }
                    else{
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] + 1);
                            sys->Acval.push_back(-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                        }
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + sys->N_cell_y + 1) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] + sys->N_cell_y + 1);
                        sys->Acval.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]);
                    }
                }
                else if (inx == sys->N_cell_x){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] - sys->N_cell_y - 1) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] - sys->N_cell_y - 1);
                        sys->Acval.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]);
                    }
                    if (iny != 0 && iny != sys->N_cell_y){
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->portNno.find(sys->conductor[i].node[j] - 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] - 1);
                            sys->Acval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                        }
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] + 1);
                            sys->Acval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                        }
                    }
                    else if (iny == sys->N_cell_y){
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->portNno.find(sys->conductor[i].node[j] - 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] - 1);
                            sys->Acval.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                        }
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                    }
                    else{
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / (sys->xn[inx] - sys->xn[inx - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] + 1);
                            sys->Acval.push_back(-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                        }
                    }
                }
                else{
                    if (iny != 0 && iny != sys->N_cell_y){
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->portNno.find(sys->conductor[i].node[j] - 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] - 1);
                            sys->Acval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                        }
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] + 1);
                            sys->Acval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                        }
                    }
                    else if (iny == sys->N_cell_y){
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1] != 0 && sys->portNno.find(sys->conductor[i].node[j] - 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] - 1);
                            sys->Acval.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]);
                        }
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / (sys->yn[iny] - sys->yn[iny - 1]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                    }
                    else{
                        if (inz != 0 && inz != sys->N_cell_z){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]
                                + 2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else if (inz == 0){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        else{
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j]);
                            sys->Acval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]
                                + 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]
                                + 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / (sys->zn[inz] - sys->zn[inz - 1]) * sys->sig[(inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny]);
                        }
                        if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + 1) == sys->portNno.end()){
                            sys->AcRowId.push_back(sys->conductor[i].node[j]);
                            sys->AcColId.push_back(sys->conductor[i].node[j] + 1);
                            sys->Acval.push_back(-1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / (sys->yn[iny + 1] - sys->yn[iny]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny]);
                        }
                    }
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + sys->N_cell_y + 1) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] + sys->N_cell_y + 1);
                        sys->Acval.push_back(-1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / (sys->xn[inx + 1] - sys->xn[inx]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*inx + iny]);
                    }
                }

                if (inz != 0 && inz != sys->N_cell_z){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + sys->N_node_s) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] + sys->N_node_s);
                        sys->Acval.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny]);
                    }
                }
                else if (inz == 0){
                    if (sys->markEdge[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny] != 0 && sys->portNno.find(sys->conductor[i].node[j] + sys->N_node_s) == sys->portNno.end()){
                        sys->AcRowId.push_back(sys->conductor[i].node[j]);
                        sys->AcColId.push_back(sys->conductor[i].node[j] + sys->N_node_s);
                        sys->Acval.push_back(-1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / (sys->zn[inz + 1] - sys->zn[inz]) * sys->sig[inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx * (sys->N_cell_y + 1) + iny]);
                    }
                }
 
            }
            cindex[count + 1] = cindex[count] + cindex[count + 1];
            count++;
        }
    }

    for (i = 0; i < sys->AcRowId.size(); i++){
        sys->AcRowId[i] = map[sys->AcRowId[i]] - 1;
        sys->AcColId[i] = map[sys->AcColId[i]] - 1;
    }
 
    int leng_v0c2 = 0;
    for (i = 0; i < sys->numPorts; i++){
        temp2.clear();
        for (j = 0; j < sys->portCoor[i].nodenum; j++){
            temp2.push_back(sys->portCoor[i].node[j]);
        }
        rowId.clear();
        colId.clear();
        val.clear();
        status = nodeAdd(rowId, colId, val, temp2, sys->N_edge, sys);
        if (status != 0)
            return status;
        for (j = 0; j < rowId.size(); j++){
            sys->v0c2ColId.push_back(leng_v0c2);
            sys->v0c2RowId.push_back(rowId[j]);
            sys->v0c2val.push_back(val[j]);
        }
        leng_v0c2++;
    }
    /*for (i = 0; i < sys->v0c2RowId.size(); i++){
        cout << sys->v0c2RowId[i] << " " << sys->v0c2ColId[i] << " " << sys->v0c2val[i] << endl;
    }*/

    sys->v0cvalo = sys->v0cval;    // v0cvalo is the v0c values without D_sig
    for (i = 0; i < sys->v0cColId.size(); i++){
        sys->v0cval[i] = sys->v0cval[i] * sqrt(sys->sig[sys->v0cRowId[i]]);       // Compute the sparse form of D_sig*V0c
    }
    sys->v0cavalo = sys->v0caval;
    for (i = 0; i < sys->v0caColId.size(); i++){
        sys->v0caval[i] = sys->v0caval[i] * sqrt(sys->sig[sys->v0caRowId[i]]);
    }

    int portCdt = 0;    // how many conductors contain ports
    for (i = 0; i < sys->numCdt; i++){
        portCdt = portCdt + sys->conductor[i].markPort;
    }

    /*outfile1.open("v0c.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->v0cRowId.size(); i++){
        outfile1 << sys->v0cRowId[i] + 1 << " " << sys->v0cColId[i] + 1 << " " << sys->v0cvalo[i] << endl;
    }
    outfile1.close();
    outfile1.open("v0ca.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->v0caval.size(); i++){
        outfile1 << sys->v0caRowId[i] + 1 << " " << sys->v0caColId[i] + 1 << " " << sys->v0cavalo[i] << endl;
    }
    outfile1.close();*/

    /* Compute the matrix V0c'*D_sig*V0c */
    /*status = matrixMulti(sys->v0caColId, sys->v0caRowId, sys->v0caval, sys->v0cRowId, sys->v0cColId, sys->v0cval, sys->AcRowId, sys->AcColId, sys->Acval);
    if (status != 0)
        return status;*/
    status = COO2CSR(sys->AcRowId, sys->AcColId, sys->Acval);
    if (status != 0)
        return status;

    double *a = &(sys->Acval[0]);
    int *ia = &(sys->AcRowId[0]);
    int *ja = &(sys->AcColId[0]);

    /* Construct V0d with row id, col id and its val */
    int leng_v0d1 = 0;    // store the num of v0d1 vectors, which are nodes outside the conductors
    int leng_v0d1a = 0;

    for (i = 0; i < sys->N_node; i++){
        if (sys->markNode[i] == 0){
            for (j = 0; j < sys->nodeEdge[i].size(); j++){
                sys->v0d1RowId.push_back(sys->nodeEdge[i][j].first);
                sys->v0d1ColId.push_back(leng_v0d1);
                sys->v0d1val.push_back(sys->nodeEdge[i][j].second);
            }
            leng_v0d1++;

            inz = i / sys->N_node_s;
            inx = (i - inz * sys->N_node_s) / (sys->N_cell_y + 1);
            iny = i % (sys->N_cell_y + 1);
            if (inz != 0 && inz != sys->N_cell_z){
                sys->v0d1aRowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
            }
            else if (inz == sys->N_cell_z){
                sys->v0d1aRowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
            }

            if (iny != 0 && iny != sys->N_cell_y){
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
            }
            else if (iny == sys->N_cell_y){
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
            }
            else{
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
            }

            if (inx != 0 && inx != sys->N_cell_x){
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
            }
            else if (inx == sys->N_cell_x){
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
            }
            else{
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
            }

            if (inz != sys->N_cell_z && inz != 0){
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
            }
            else if (inz == 0){
                sys->v0d1aRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                sys->v0d1aColId.push_back(leng_v0d1a);
                sys->v0d1aval.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
            }

            leng_v0d1a++;

        }
    }

    mark = sys->v0d1ColId.back();
    while (sys->v0d1ColId.back() == mark){
        sys->v0d1RowId.pop_back();   // one node removed to make V0 full rank
        sys->v0d1ColId.pop_back();
        sys->v0d1val.pop_back();
    }
    leng_v0d1--;
    mark = sys->v0d1aColId.back();
    while (sys->v0d1aColId.back() == mark){
        sys->v0d1aRowId.pop_back();
        sys->v0d1aColId.pop_back();
        sys->v0d1aval.pop_back();
    }
    leng_v0d1a--;

    int leng_v0d2 = 0;
    int leng_v0d2a = 0;
    /*if (sys->numCdt - portCdt > 0){
        vector<vector<int> > ind(sys->numCdt);
        k = 0;
        for (i = 0; i < sys->numCdt; i++){
            for (j = 0; j < sys->cdtNumNode[i]; j++){
                ind[k].push_back(sys->conductor[i].node[j]);
            }
            k++;
        }
        for (i = 0; i < sys->numCdt; i++){
            if (sys->conductor[i].markPort == 0){    // if this conductor doesn't contain any port
                rowId.clear();
                colId.clear();
                val.clear();
                status = nodeAdd(rowId, colId, val, ind[i], sys->N_edge, sys);
                if (status != 0)
                    return status;
                for (j = 0; j < rowId.size(); j++){
                    sys->v0d2RowId.push_back(rowId[j]);
                    sys->v0d2ColId.push_back(leng_v0d2);
                    sys->v0d2val.push_back(val[j]);
                }
                leng_v0d2++;

                rowId.clear();
                colId.clear();
                val.clear();
                status = nodeAddAvg(rowId, colId, val, ind[i][0], sys->N_edge, sys);
                if (status != 0)
                    return status;
                for (j = 0; j < rowId.size(); j++){
                    sys->v0d2aRowId.push_back(rowId[j]);
                    sys->v0d2aColId.push_back(leng_v0d2a);
                    sys->v0d2aval.push_back(val[j]);
                }
                leng_v0d2a++;
            }
        }
        ind.clear();
    }
    outfile1.open("v0d2a.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->v0d2aRowId.size(); i++){
        outfile1 << sys->v0d2aRowId[i] + 1 << " " << sys->v0d2aColId[i] + 1 << " " << sys->v0d2aval[i] << endl;
    }
    outfile1.close();
    cout << leng_v0d1a << endl;

    outfile1.open("v0d1.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->v0d1RowId.size(); i++){
        outfile1 << sys->v0d1RowId[i] + 1 << " " << sys->v0d1ColId[i] + 1 << " " << sys->v0d1val[i] << endl;
    }
    outfile1.close();
    cout << leng_v0d1 << endl;*/

    /* Generate the merged V0d1 */
    /*double dx = 3e-6, dy = 1e-6; // merge layer by layer
    int nx = ceil((sys->xlim2 - sys->xlim1) / dx), ny = ceil((sys->ylim2 - sys->ylim1) / dy);
    int mnno = 0;    // merged node number
    unordered_map<int, vector<int> > group;
    int leng_v0d1 = 0;
    double sta;

    group.clear();
    for (i = 0; i < sys->numStack + 1; i++){
        for (j = 0; j < (sys->N_cell_y + 1); j++){
            sta = sys->xlim1;
            k = 0;
            while (k < (sys->N_cell_x + 1)){
                while (sys->markNode[i * sys->N_node_s + k * (sys->N_cell_y + 1) + j] == 0 && abs(sys->nodepos[(i * sys->N_node_s + k * (sys->N_cell_y + 1) + j) * 3] - sta) < dx){
                    group[mnno].push_back(i * sys->N_node_s + k * (sys->N_cell_y + 1) + j);
                    k++;
                }
                mnno++;
                while (sys->markNode[i * sys->N_node_s + k * (sys->N_cell_y + 1) + j] == 1 && k < (sys->N_cell_x + 1)){
                    k++;
                }
                sta = sys->nodepos[(i * sys->N_node_s + k * (sys->N_cell_y + 1) + j) * 3];
            }
        }
    }
    double *v;
    for (auto g : group){
        rowId.clear();
        colId.clear();
        val.clear();
        //status = nodeAdd(rowId, colId, val, g.second, sys->N_edge, sys);

        v = (double*)calloc(sys->N_edge, sizeof(double));
        for (l = 0; l < g.second.size(); l++){
            for (n = 0; n < sys->nodeEdge[g.second[l]].size(); n++){
                v[sys->nodeEdge[g.second[l]][n].first] = v[sys->nodeEdge[g.second[l]][n].first] + sys->nodeEdge[g.second[l]][n].second;
            }
        }
        for (n = 0; n < sys->N_edge; n++){
            if (abs(v[n]) > 1e-5) {
                rowId.push_back(n);
                colId.push_back(0);
                val.push_back(v[n]);
            }
        }
        free(v);
        v = NULL;

        if (status != 0)
            return status;
        for (j = 0; j < rowId.size(); j++){
            sys->v0d1RowId.push_back(rowId[j]);
            sys->v0d1ColId.push_back(leng_v0d1);
            sys->v0d1val.push_back(val[j]);
        }
        leng_v0d1++;
    }
    group.clear();
    cout << leng_v0d1 << endl;*/

    /*outfile1.open("v0d1.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->v0d1RowId.size(); i++){
        outfile1 << sys->v0d1RowId[i]+1 << " " << sys->v0d1ColId[i]+1 << " " << sys->v0d1val[i] << endl;
    }
    outfile1.close();*/

    sys->v0d1valo = sys->v0d1val;
    for (i = 0; i < sys->v0d1RowId.size(); i++){    // compute sqrt(D_eps)*V0d1
        sys->v0d1val[i] = sys->v0d1val[i] * sqrt(sys->eps[sys->v0d1RowId[i]]);
    }
    sys->v0d1avalo = sys->v0d1aval;
    for (i = 0; i < sys->v0d1aRowId.size(); i++){
        sys->v0d1aval[i] = sys->v0d1aval[i] * sqrt(sys->eps[sys->v0d1aRowId[i]]);
    }

    if (leng_v0d2 > 0){
        sys->v0d2valo = sys->v0d2val;
        for (i = 0; i < sys->v0d2RowId.size(); i++){    // compute sqrt(D_eps)*v0d2
            sys->v0d2val[i] = sys->v0d2val[i] * sqrt(sys->eps[sys->v0d2RowId[i]]);
        }
        sys->v0d2avalo = sys->v0d2aval;
        for (i = 0; i < sys->v0d2aRowId.size(); i++){
            sys->v0d2aval[i] = sys->v0d2aval[i] * sqrt(sys->eps[sys->v0d2aRowId[i]]);
        }
        sys->v0d2ColIdo = sys->v0d2ColId;
        status = COO2CSR(sys->v0d2ColId, sys->v0d2RowId, sys->v0d2val);
        if (status != 0)
            return status;
        sys->v0d2aColIdo = sys->v0d2aColId;
        status = COO2CSR(sys->v0d2aColId, sys->v0d2aRowId, sys->v0d2aval);
        if (status != 0)
            return status;
    }

    int leng_v0d = leng_v0d1 + leng_v0d2;
    status = COO2CSR(sys->v0d1ColId, sys->v0d1RowId, sys->v0d1val);
    if (status != 0)
        return status;
    status = COO2CSR(sys->v0d1aColId, sys->v0d1aRowId, sys->v0d1aval);
    if (status != 0)
        return status;

    sys->v0cColIdo = sys->v0cColId;
    status = COO2CSR(sys->v0cColId, sys->v0cRowId, sys->v0cval);
    if (status != 0)
        return status;
    sys->v0caColIdo = sys->v0caColId;
    status = COO2CSR(sys->v0caColId, sys->v0caRowId, sys->v0caval);
    if (status != 0)
        return status;

    sys->v0c2ColIdo = sys->v0c2ColId;
    status = COO2CSR(sys->v0c2ColId, sys->v0c2RowId, sys->v0c2val);
    if (status != 0)
        return status;

    /* Pick up a y0c2 cooresponding to one source port */
    ofstream outfile;
    double freq;
    double *bd1, *bd2;
    double *bdc1, *bdc2;
    double *solution_b1;
    double *temp;
    double mdone;
    int ione;
    double *solution_d2;
    double *v0d2, *b0d2;
    double *v0d2epsv0d2, *v0d2epsv0d1;
    double *v0d2epsv0d1b1;
    int *ipiv;
    int info;
    double *workspace;
    double *xd2;
    double *xd1rhs_re, *xd1rhs_im;
    double *xd1;
    double *temp1;
    int startCol;
    startCol = 0;
    /*sys->x = (complex<double>*) malloc(sys->numPorts * sys->numPorts * sizeof(complex<double>));
    for (i = 0; i < sys->numPorts*sys->numPorts; i++){
        sys->x[i] = complex<double> (0., 0.); // Complex double constructor from real and imaginary
    }*/
    void *pt[64];
    int mtype;
    int iparm[64];
    double dparm[64];
    int maxfct, mnum, phase, error, solver;
    int num_proces;   //number of processors
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

    double *b, *xc;    // the array of the right hand side
    xcol = 0;

    vector<int> v0csRowId, v0csColId;
    vector<double> v0csval;
    double *epsv0dy0d, *epsv0cy0c, *sigv0cy0c;
    double *v0csy0d, *v0csy0c, *v0cssigy0c;
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
    sys->Y = (complex<double>*)calloc(sys->numPorts * sys->numPorts, sizeof(complex<double>));


    char transa;
    int m;
    double alpha = 1, beta = 0;
    char matdescra[6];
    matdescra[0] = 'G'; matdescra[3] = 'C';    // general matrix multi, 0-based indexing
    cout << "Begin!\n";
    nrhs = leng_v0d2;
    v0d2epsv0d2 = (double*)calloc(nrhs*nrhs, sizeof(double));
    v0d2epsv0d1 = (double*)calloc(nrhs*nrhs, sizeof(double));
    solution_d2 = (double*)calloc(leng_v0d2*leng_v0d1, sizeof(double));
    v0d2 = (double*)calloc(leng_v0d2*sys->N_edge, sizeof(double));
    b0d2 = (double*)calloc(leng_v0d2*leng_v0d1, sizeof(double));
    if (leng_v0d2 != 0){
        transa = 'N';
        m = leng_v0d1;
        k = sys->N_edge;

        for (i = 0; i < sys->v0d2ColIdo.size(); i++){
            v0d2[sys->v0d2ColIdo[i] * sys->N_edge + sys->v0d2RowId[i]] = sys->v0d2val[i];
        }

        for (i = 1; i <= nrhs; i++){
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1aval[0], &sys->v0d1aRowId[0], &sys->v0d1aColId[0], &sys->v0d1aColId[1], &v0d2[(i - 1)*sys->N_edge], &beta, &b0d2[(i - 1)*leng_v0d1]);
        }

        for (i = 0; i < nrhs; i++){
            //clock_t t1 = clock();
            status = interativeSolver(leng_v0d1, 1, &b0d2[i*leng_v0d1], &sys->v0d1ColId[0], &sys->v0d1RowId[0], &sys->v0d1val[0], &sys->v0d1aColId[0], &sys->v0d1aRowId[0], &sys->v0d1aval[0], &solution_d2[i*leng_v0d1], sys);
            if (status != 0)
                return status;
            //cout << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
        }
        cout << "The first iteration is done!\n";

        temp1 = (double*)malloc(sys->N_edge*nrhs*sizeof(double));
        m = nrhs;
        transa = 'N';
        for (i = 0; i < nrhs; i++){
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
    }

    double *y0c2;
    sourcePort = 0;

    while (sourcePort < sys->numPorts){
        freq = 10.e+09;
        y0c2 = (double*)calloc(leng_v0c2, sizeof(double));
        sys->v0c2y0c2 = (double*)malloc(sys->N_edge*sizeof(double));
        sys->v0c2y0c2o = (double*)malloc(sys->N_edge*sizeof(double));
        y0c2[sourcePort] = 1;
        transa = 'T';
        m = leng_v0c2;
        k = sys->N_edge;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0c2val[0], &sys->v0c2RowId[0], &sys->v0c2ColId[0], &sys->v0c2ColId[1], y0c2, &beta, sys->v0c2y0c2);
        for (i = 0; i < sys->N_edge; i++){
            //if (isnan(sys->sig[i])){
                //cout << sys->sig[i] << endl;}
            sys->v0c2y0c2o[i] = sys->v0c2y0c2[i];
            sys->v0c2y0c2[i] = -sys->v0c2y0c2[i] * sqrt(sys->sig[i]);
            /*if (sys->v0c2y0c2[i] != 0)
            {
                cout << sys->v0c2y0c2[i] << " " << sys->sig[i] << "#" << endl;
            }*/
        }
        cout << endl;

        /* Compute C right hand side */
        xc = (double*)calloc(leng_v0c, sizeof(double));

        // Solve conductor equation
        //clock_t t = clock();
        /*vector<int> rowIda;
        vector<int> colIda;
        vector<double> vala;
        for (i = 0; i < sys->numCdt; i++){
            rowId.clear();
            colId.clear();
            val.clear();
            rowId.assign(sys->v0cRowId.begin() + cindex[i] + 1, sys->v0cRowId.begin() + cindex[i + 1]);
            colId.assign(sys->v0cColIdo.begin() + cindex[i] + 1, sys->v0cColIdo.begin() + cindex[i + 1]);
            val.assign(sys->v0cval.begin() + cindex[i] + 1, sys->v0cval.begin() + cindex[i + 1]);
            mark = colId[0];
            for (j = 0; j < rowId.size(); j++){
                colId[j] = colId[j] - mark;
            }

            rowIda.clear();
            colIda.clear();
            vala.clear();
            rowIda.assign(sys->v0caRowId.begin() + cindex[i] + 1, sys->v0caRowId.begin() + cindex[i + 1]);
            colIda.assign(sys->v0caColId.begin() + cindex[i] + 1, sys->v0caColId.begin() + cindex[i + 1]);
            vala.assign(sys->v0caval.begin() + cindex[i] + 1, sys->v0caval.begin() + cindex[i + 1]);
            mark = colIda[0];
            for (j = 0; j < rowIda.size(); j++){
                colIda[j] = colIda[j] - mark;
            }

            status = COO2CSR(colId, rowId, val);
            status = COO2CSR(colIda, rowIda, vala);
            transa = 'N';
            m = acu_cnno[i + 1] - acu_cnno[i];
            k = sys->N_edge;
            crhs = (double*)malloc(m * sizeof(double));
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &vala[0], &rowIda[0], &colIda[0], &colIda[1], sys->v0c2y0c2, &beta, crhs);
            status = interativeSolver(m, 1, crhs, &colId[0], &rowId[0], &val[0], &colIda[0], &rowIda[0], &vala[0], &xc[acu_cnno[i + 1] - m], sys);
            if (status != 0)
                return status;
            free(crhs);
            crhs = NULL;
        }
        cout << "The conductor part iteration is done!\n";*/

        //cout << (clock() - t) * 1.0 / CLOCKS_PER_SEC << endl;
        crhs = (double*)malloc(leng_v0c * sizeof(double));
        transa = 'N';
        m = leng_v0ca;
        k = sys->N_edge;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0caval[0], &sys->v0caRowId[0], &sys->v0caColId[0], &sys->v0caColId[1], sys->v0c2y0c2, &beta, crhs);
        /*for (i = 0; i < leng_v0c; i++){
            cout << crhs[i] << " ";
        }
        cout << endl;*/

        pardisoinit(pt, &mtype, iparm);
        nrhs = 1;
        iparm[38] = 1;
        iparm[34] = 1;    //0-based indexing
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &leng_v0c, a, ia, ja, &perm, &nrhs, iparm, &msglvl, crhs, xc, &error);

        /*for (i = 0; i < leng_v0c; i++){
            cout << xc[i] << " ";
        }
        cout << endl;*/

        /* Compute the dielectric right hand side */
        transa = 'T';
        m = leng_v0c;
        k = sys->N_edge;
        mdone = 1;
        ione = 1;
        n = leng_v0d1;
        nrhs = 1;
        solution_b1 = (double*)calloc(n, sizeof(double));
        bd1 = (double*)calloc(n, sizeof(double));
        bd2 = (double*)calloc(leng_v0d2, sizeof(double));
        sys->yc = (double*)calloc(sys->N_edge, sizeof(double));
        //if (temp == NULL)
        temp = (double*)calloc(sys->N_edge, sizeof(double));
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0cvalo[0], &sys->v0cRowId[0], &sys->v0cColId[0], &sys->v0cColId[1], xc, &beta, temp);

        sys->v0cy0c = sys->v0c2y0c2o;
        cblas_daxpy(sys->N_edge, mdone, temp, ione, sys->v0cy0c, ione);
        for (i = 0; i < sys->N_edge; i++){
            sys->yc[i] = sys->v0cy0c[i];
            sys->v0cy0c[i] = -sys->v0cy0c[i] * sqrt(sys->eps[i]);
        }

        transa = 'N';
        m = leng_v0d1;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1aval[0], &sys->v0d1aRowId[0], &sys->v0d1aColId[0], &sys->v0d1aColId[1], sys->v0cy0c, &beta, bd1);
        m = leng_v0d2;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2aval[0], &sys->v0d2aRowId[0], &sys->v0d2aColId[0], &sys->v0d2aColId[1], sys->v0cy0c, &beta, bd2);


        //ofstream outfile;
        /*outfile.open("yc.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < sys->N_edge; i++){
            outfile << sys->yc[i] << "\n";
        }
        outfile.close();*/


        /* Solve (V0d1'*D_eps*V0d1)\b1 */
        status = interativeSolver(leng_v0d1, 1, bd1, &sys->v0d1ColId[0], &sys->v0d1RowId[0], &sys->v0d1val[0], &sys->v0d1aColId[0], &sys->v0d1aRowId[0], &sys->v0d1aval[0], solution_b1, sys);
        if (status != 0)
            return status;
        cout << "The second iteration is done!\n";
        /* Solve (V0d1'*D_eps*V0d1)\(V0d1'*D_eps*V0d2) */


        /* Get the solution of xd2 */
        v0d2epsv0d1b1 = (double*)calloc(leng_v0d2, sizeof(double));
        xd2 = (double*)calloc(leng_v0d2, sizeof(double));
        if (leng_v0d2 != 0){

            transa = 'T';
            m = leng_v0d1;
            temp = (double*)malloc(sys->N_edge*sizeof(double));
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1val[0], &sys->v0d1RowId[0], &sys->v0d1ColId[0], &sys->v0d1ColId[1], solution_b1, &beta, temp);
            transa = 'N';
            m = nrhs;
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2aval[0], &sys->v0d2aRowId[0], &sys->v0d2aColId[0], &sys->v0d2aColId[1], temp, &beta, v0d2epsv0d1b1);
            mdone = -1;
            cblas_daxpy(nrhs, mdone, v0d2epsv0d1b1, ione, bd2, ione);
            free(temp);
            temp = NULL;

            for (i = 0; i < leng_v0d2; i++){
                for (j = 0; j < leng_v0d2; j++){
                    xd2[i] += v0d2epsv0d2[i + j*leng_v0d2] * bd2[j];
                }
            }
        }

        xd1 = (double*)calloc(leng_v0d1, sizeof(double));
        for (i = 0; i < leng_v0d2; i++){
            for (j = 0; j < leng_v0d1; j++){
                xd1[j] -= solution_d2[j + i*leng_v0d1] * xd2[i];
            }
        }

        /* Get the solution xd1 */
        mdone = 1;
        cblas_daxpy(leng_v0d1, mdone, solution_b1, ione, xd1, ione);

        /* Form the final yd */
        sys->yd = (double*)malloc(sys->N_edge*sizeof(double));
        transa = 'T';
        m = leng_v0d1;
        k = sys->N_edge;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d1valo[0], &sys->v0d1RowId[0], &sys->v0d1ColId[0], &sys->v0d1ColId[1], xd1, &beta, sys->yd);
        /*temp = (double*)malloc(sys->N_edge*sizeof(double));
        m = leng_v0d2;
        mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &sys->v0d2valo[0], &sys->v0d2RowId[0], &sys->v0d2ColId[0], &sys->v0d2ColId[1], xd2, &beta, temp);
        mdone = 1;
        cblas_daxpy(sys->N_edge, mdone, temp, ione, sys->yd, ione);

        free(temp);
        temp = NULL;*/

        /*outfile.open("yc.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < sys->N_edge; i++){
        outfile << sys->yc[i] << "\n";
        }
        outfile.close();*/


        /* Form the final y */
        sys->y = (double*) malloc(sys->N_edge * sizeof(double));
        for (i = 0; i < sys->N_edge; i++){
            sys->y[i] = sys->yd[i] + sys->yc[i];
            //cout << sys->y[i] << " ";
        }
        /*outfile.open("yd.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < sys->N_edge; i++){
            outfile << sys->yd_re[i] <<" "<<sys->yd_im[i]<< "\n";
        }
        outfile.close();*/

        xcol++;

        for (i = 0; i < sys->numPorts; i++){
            port = i + 1;
            v0csRowId.clear();
            v0csColId.clear();
            v0csval.clear();
            v0csin = 0;

            for (j = 0; j < sys->portCoor[i].nodenum; j++){
                inz = sys->portCoor[i].node[j] / sys->N_node_s;
                inx = (sys->portCoor[i].node[j] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->portCoor[i].node[j] % (sys->N_cell_y + 1);
                if (inz != 0 && inz != sys->N_cell_z){
                    v0csRowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                }
                else if (inz == sys->N_cell_z){
                    v0csRowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                }

                if (iny != 0 && iny != sys->N_cell_y){
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                }
                else if (iny == sys->N_cell_y){
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                }
                else{
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                }

                if (inx != 0 && inx != sys->N_cell_x){
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                }
                else if (inx == sys->N_cell_x){
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                }
                else{
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                }

                if (inz != sys->N_cell_z && inz != 0){
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                }
                else if (inz == 0){
                    v0csRowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    v0csColId.push_back(v0csin);
                    v0csval.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                }
                v0csin++;
            }

            leng_v0cs = v0csin;
            epsv0dy0d = (double*)malloc(sys->N_edge*sizeof(double));
            epsv0cy0c = (double*)malloc(sys->N_edge*sizeof(double));
            sigv0cy0c = (double*)malloc(sys->N_edge*sizeof(double));
            v0csy0d = (double*)malloc(leng_v0cs*sizeof(double));
            v0csy0c = (double*)malloc(leng_v0cs*sizeof(double));
            v0cssigy0c = (double*)malloc(leng_v0cs*sizeof(double));
            status = COO2CSR(v0csColId, v0csRowId, v0csval);
            if (status != 0)
                return status;

            for (j = 0; j < sys->N_edge; j++){
                epsv0dy0d[j] = 2 * PI*freq*sys->eps[j] * sys->yd[j];
            }
            transa = 'N';
            m = leng_v0cs;
            k = sys->N_edge;
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &v0csval[0], &v0csRowId[0], &v0csColId[0], &v0csColId[1], epsv0dy0d, &beta, v0csy0d);

            for (j = 0; j < sys->N_edge; j++){
                epsv0cy0c[j] = 2 * PI*freq*sys->eps[j] * sys->yc[j];
            }
            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &v0csval[0], &v0csRowId[0], &v0csColId[0], &v0csColId[1], epsv0cy0c, &beta, v0csy0c);

            mdone = 1;
            cblas_daxpy(leng_v0cs, mdone, v0csy0c, ione, v0csy0d, ione);    // the addition of v0csy0d and v0csy0c


            for (j = 0; j < sys->N_edge; j++){
                if (sys->sig[j] != 0)
                    sigv0cy0c[j] = sys->sig[j] * sys->yc[j];
            }

            mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, &v0csval[0], &v0csRowId[0], &v0csColId[0], &v0csColId[1], sigv0cy0c, &beta, v0cssigy0c);

            sys->v0csJ = (complex<double>*)malloc(leng_v0cs * sizeof(complex<double>));
            for (j = 0; j < leng_v0cs; j++){
                sys->v0csJ[j] = -v0cssigy0c[j] - 1i * v0csy0d[j]; // Complex<double> literal
            }
            node = sys->portCoor[i].node[0];
            int node1;
            double a = 0, area, areasum;
            current = complex<double> (0., 0.);

if (sys->portCoor[i].x1 == sys->portCoor[i].x2){

    inx = xi[sys->portCoor[i].x1];

    for (l = zi[sys->portCoor[i].z1]; l <= zi[sys->portCoor[i].z2]; l++){
        for (k = yi[sys->portCoor[i].y1]; k <= yi[sys->portCoor[i].y2]; k++){
            mark = (k - yi[sys->portCoor[i].y1])
                + (l - zi[sys->portCoor[i].z1])*(yi[sys->portCoor[i].y2] - yi[sys->portCoor[i].y2] + 1);
            a = 0;
            areasum = 0;
            if (k == 0){
                a += 1 / (sys->yn[k + 1] - sys->yn[k]);
                area = 1;
                if (inx == 0){
                    area *= (sys->xn[inx + 1] - sys->xn[inx]);
                }
                else if (inx == sys->N_cell_x){
                    area *= (sys->xn[inx] - sys->xn[inx - 1]);
                }
                else{
                    area *= (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area;
            }
            else if (k == sys->N_cell_y){
                a += 1 / (sys->yn[k] - sys->yn[k - 1]);
                area = 1;
                if (inx == 0){
                    area *= (sys->xn[inx + 1] - sys->xn[inx]);
                }
                else if (inx == sys->N_cell_x){
                    area *= (sys->xn[inx] - sys->xn[inx - 1]);
                }
                else{
                    area *= (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area;
            }
            else{
                a += 4 / (sys->yn[k + 1] - sys->yn[k - 1]);
                area = 1;
                if (inx == 0){
                    area *= (sys->xn[inx + 1] - sys->xn[inx]);
                }
                else if (inx == sys->N_cell_x){
                    area *= (sys->xn[inx] - sys->xn[inx - 1]);
                }
                else{
                    area *= (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area * 2;
            }
            if (inx == 0){
                a += 1 / (sys->xn[inx + 1] - sys->xn[inx]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *= (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *= (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area;
            }
            else if (inx == sys->N_cell_x){
                a += 1 / (sys->xn[inx] - sys->xn[inx - 1]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *= (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *= (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area;
            }
            else{
                a += 4 / (sys->xn[inx + 1] - sys->xn[inx - 1]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *= (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *= (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area * 2;
            }
            if (l == 0){
                a = a + 1 / (sys->zn[l + 1] - sys->zn[l]);
                area = 1;
                if (k == 0){
                    area *=  (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *= (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *= (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (inx == 0){
                    area *= (sys->xn[inx + 1] - sys->xn[inx]);
                }
                else if (inx == sys->N_cell_x){
                    area *= (sys->xn[inx] - sys->xn[inx - 1]);
                }
                else{
                    area *= (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2;
                }
                areasum += area;
            }
            else if (l == sys->N_cell_z){
                a += 1 / (sys->zn[l] - sys->zn[l - 1]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *= (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *= (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (inx == 0){
                    area *= (sys->xn[inx + 1] - sys->xn[inx]);
                }
                else if (inx == sys->N_cell_x){
                    area *= (sys->xn[inx] - sys->xn[inx - 1]);
                }
                else{
                    area *= (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2;
                }
                areasum +=  area;
            }
            else{
                a += 4 / (sys->zn[l + 1] - sys->zn[l - 1]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *= (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *= (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (inx == 0){
                    area *= (sys->xn[inx + 1] - sys->xn[inx]);
                }
                else if (inx == sys->N_cell_x){
                    area *= (sys->xn[inx] - sys->xn[inx - 1]);
                }
                else{
                    area *= (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2;
                }
                areasum += area * 2;
            }
            sys->v0csJ[mark] = sys->v0csJ[mark] / (-a);

            current +=  (sys->v0csJ[mark]) * areasum;

        }
    }
}
else if (sys->portCoor[i].y1 == sys->portCoor[i].y2){

    iny = yi[sys->portCoor[i].y1];

    for (l = zi[sys->portCoor[i].z1]; l <= zi[sys->portCoor[i].z2]; l++){
        for (j = xi[sys->portCoor[i].x1]; j <= xi[sys->portCoor[i].x2]; j++){
            mark = (j - xi[sys->portCoor[i].x1])
                + (l - zi[sys->portCoor[i].z1])*(xi[sys->portCoor[i].x2] - xi[sys->portCoor[i].x1] + 1);
            node1 = l*(sys->N_node_s) + (sys->N_cell_y + 1) * j + yi[sys->portCoor[i].y1];
            a = 0;
            areasum = 0;
            if (iny == 0){
                a +=1 / (sys->yn[iny + 1] - sys->yn[iny]);
                area = 1;
                if (j == 0){
                    area *= (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *= (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *= (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum +=  area;
            }
            else if (iny == sys->N_cell_y){
                a += 1 / (sys->yn[iny] - sys->yn[iny - 1]);
                area = 1;
                if (j == 0){
                    area *= (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *= (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *= (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area;
            }
            else{
                a += 4 / (sys->yn[iny + 1] - sys->yn[iny - 1]);
                area = 1;
                if (j == 0){
                    area *= (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *= (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *= (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum += area * 2;
            }
            if (j == 0){
                a += 1 / (sys->xn[j + 1] - sys->xn[j]);
                area = 1;
                if (iny == 0){
                    area *=  (sys->yn[iny + 1] - sys->yn[iny]);
                }
                else if (iny == sys->N_cell_y){
                    area *= (sys->yn[iny] - sys->yn[iny - 1]);
                }
                else{
                    area *= (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *= (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum +=  area;
            }
            else if (j == sys->N_cell_x){
                a += 1 / (sys->xn[j] - sys->xn[j - 1]);
                area = 1;
                if (iny == 0){
                    area *= (sys->yn[iny + 1] - sys->yn[iny]);
                }
                else if (iny == sys->N_cell_y){
                    area *=  (sys->yn[iny] - sys->yn[iny - 1]);
                }
                else{
                    area *= (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2;
                }
                if (l == 0){
                    area *= (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *= (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *=  (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum +=  area;
            }
            else{
                a += 4 / (sys->xn[j + 1] - sys->xn[j - 1]);
                area = 1;
                if (iny == 0){
                    area *=  (sys->yn[iny + 1] - sys->yn[iny]);
                }
                else if (iny == sys->N_cell_y){
                    area *=  (sys->yn[iny] - sys->yn[iny - 1]);
                }
                else{
                    area *=  (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2;
                }
                if (l == 0){
                    area *=  (sys->zn[l + 1] - sys->zn[l]);
                }
                else if (l == sys->N_cell_z){
                    area *=  (sys->zn[l] - sys->zn[l - 1]);
                }
                else{
                    area *=  (sys->zn[l + 1] - sys->zn[l - 1]) / 2;
                }
                areasum +=  area * 2;
            }
            if (l == 0){
                a +=1 / (sys->zn[l + 1] - sys->zn[l]);
                area = 1;
                if (iny == 0){
                    area *=  (sys->yn[iny + 1] - sys->yn[iny]);
                }
                else if (iny == sys->N_cell_y){
                    area *=  (sys->yn[iny] - sys->yn[iny - 1]);
                }
                else{
                    area *= (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2;
                }
                if (j == 0){
                    area *=  (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                areasum +=  area;
            }
            else if (l == sys->N_cell_z){
                a += 1 / (sys->zn[l] - sys->zn[l - 1]);
                area = 1;
                if (iny == 0){
                    area *=  (sys->yn[iny + 1] - sys->yn[iny]);
                }
                else if (iny == sys->N_cell_y){
                    area *=  (sys->yn[iny] - sys->yn[iny - 1]);
                }
                else{
                    area *=  (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2;
                }
                if (j == 0){
                    area *=  (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                areasum +=  area;
            }
            else{
                a += 4 / (sys->zn[l + 1] - sys->zn[l - 1]);
                area = 1;
                if (iny == 0){
                    area *=  (sys->yn[iny + 1] - sys->yn[iny]);
                }
                else if (iny == sys->N_cell_y){
                    area *=  (sys->yn[iny] - sys->yn[iny - 1]);
                }
                else{
                    area *=  (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2;
                }
                if (j == 0){
                    area *= (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                areasum +=  area * 2;
            }

            sys->v0csJ[mark] = sys->v0csJ[mark] / (-a);

            current += (sys->v0csJ[mark]) * areasum;
        }
    }
}
else{

    inz = zi[sys->portCoor[i].z1];

    for (j = xi[sys->portCoor[port - 1].x1]; j <= xi[sys->portCoor[port - 1].x2]; j++){
        for (k = yi[sys->portCoor[port - 1].y1]; k <= yi[sys->portCoor[port - 1].y2]; k++){
            mark = (k - yi[sys->portCoor[port - 1].y1])
                + (j - xi[sys->portCoor[port - 1].x1])*(yi[sys->portCoor[port - 1].y2] - yi[sys->portCoor[port - 1].y1] + 1);

            a = 0;
            areasum = 0;
            if (k == 0){
                a +=1 / (sys->yn[k + 1] - sys->yn[k]);
                area = 1;
                if (j == 0){
                    area *=  (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                if (l == 0){
                    area *=  (sys->zn[inz + 1] - sys->zn[inz]);
                }
                else if (inz == sys->N_cell_z){
                    area *=  (sys->zn[inz] - sys->zn[inz - 1]);
                }
                else{
                    area *=  (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2;
                }
                areasum += area;
            }
            else if (k == sys->N_cell_y){
                a += 1 / (sys->yn[k] - sys->yn[k - 1]);
                area = 1;
                if (j == 0){
                    area *= (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *= (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                if (inz == 0){
                    area = area * (sys->zn[inz + 1] - sys->zn[inz]);
                }
                else if (inz == sys->N_cell_z){
                    area *= (sys->zn[inz] - sys->zn[inz - 1]);
                }
                else{
                    area *= (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2;
                }
                areasum +=  area;
            }
            else{
                a += 4 / (sys->yn[k + 1] - sys->yn[k - 1]);
                area = 1;
                if (j == 0){
                    area *=  (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *= (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                if (inz == 0){
                    area *=  (sys->zn[inz + 1] - sys->zn[inz]);
                }
                else if (inz == sys->N_cell_z){
                    area *=  (sys->zn[inz] - sys->zn[inz - 1]);
                }
                else{
                    area *=  (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2;
                }
                areasum +=  area * 2;
            }
            if (j == 0){
                a += 1 / (sys->xn[j + 1] - sys->xn[j]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *=  (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *=  (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (l == 0){
                    area*= (sys->zn[inz + 1] - sys->zn[inz]);
                }
                else if (inz == sys->N_cell_z){
                    area *=  (sys->zn[inz] - sys->zn[inz - 1]);
                }
                else{
                    area *=  (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2;
                }
                areasum +=  area;
            }
            else if (j == sys->N_cell_x){
                a += 1 / (sys->xn[j] - sys->xn[j - 1]);
                area = 1;
                if (k == 0){
                    area *=  (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *=  (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *=  (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (inz == 0){
                    area *= (sys->zn[inz + 1] - sys->zn[inz]);
                }
                else if (inz == sys->N_cell_z){
                    area *=  (sys->zn[inz] - sys->zn[inz - 1]);
                }
                else{
                    area *=  (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2;
                }
                areasum +=  area;
            }
            else{
                a +=4 / (sys->xn[j + 1] - sys->xn[j - 1]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *=  (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *=  (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (inz == 0){
                    area *=  (sys->zn[inz + 1] - sys->zn[inz]);
                }
                else if (inz == sys->N_cell_z){
                    area *=  (sys->zn[inz] - sys->zn[inz - 1]);
                }
                else{
                    area *=  (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2;
                }
                areasum += area * 2;
            }
            if (inz == 0){
                a += 1 / (sys->zn[inz + 1] - sys->zn[inz]);
                area = 1;
                if (k == 0){
                    area *= (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *=  (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *=  (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (j == 0){
                    area *=  (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                areasum +=  area;
            }
            else if (inz == sys->N_cell_z){
                a +=  1 / (sys->zn[inz] - sys->zn[inz - 1]);
                area = 1;
                if (k == 0){
                    area *=  (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *=  (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *=  (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (j == 0){
                    area *=  (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                areasum +=  area;
            }
            else{
                a +=4 / (sys->zn[inz + 1] - sys->zn[inz - 1]);
                area = 1;
                if (k == 0){
                    area *=  (sys->yn[k + 1] - sys->yn[k]);
                }
                else if (k == sys->N_cell_y){
                    area *= (sys->yn[k] - sys->yn[k - 1]);
                }
                else{
                    area *=  (sys->yn[k + 1] - sys->yn[k - 1]) / 2;
                }
                if (j == 0){
                    area *=  (sys->xn[j + 1] - sys->xn[j]);
                }
                else if (j == sys->N_cell_x){
                    area *=  (sys->xn[j] - sys->xn[j - 1]);
                }
                else{
                    area *=  (sys->xn[j + 1] - sys->xn[j - 1]) / 2;
                }
                areasum +=  area * 2;
            }

            sys->v0csJ[mark] = sys->v0csJ[mark] / (-a);

            current += (sys->v0csJ[mark]) * areasum;
        }
    }
}
            sys->Y[sourcePort * sys->numPorts + i] = current;
        }
        free(y0c2); y0c2 = NULL;
        free(sys->v0c2y0c2); sys->v0c2y0c2 = NULL;
        free(sys->v0c2y0c2o); sys->v0c2y0c2o = NULL;
        free(xc); xc = NULL;
        free(solution_b1); solution_b1 = NULL;
        free(bd1); bd1 = NULL;
        free(bd2); bd2 = NULL;
        free(xd2); xd2 = NULL;
        free(xd1); xd1 = NULL;
        free(v0d2epsv0d1b1); v0d2epsv0d1b1 = NULL;
        free(sys->yd); sys->yd = NULL;
        free(sys->yc); sys->yc = NULL;
        free(sys->y); sys->y = NULL;
        free(epsv0dy0d); epsv0dy0d = NULL;
        free(epsv0cy0c); epsv0cy0c = NULL;
        free(sigv0cy0c); sigv0cy0c = NULL;
        free(v0csy0d); v0csy0d = NULL;
        free(v0csy0c); v0csy0c = NULL;
        free(v0cssigy0c); v0cssigy0c = NULL;

        sourcePort++;

    }
    free(v0d2epsv0d2); v0d2epsv0d2 = NULL;
    free(v0d2epsv0d1); v0d2epsv0d1 = NULL;
    free(v0d2); v0d2 = NULL;
    free(b0d2); b0d2 = NULL;
    free(solution_d2); solution_d2 = NULL;
    free(cindex); cindex = NULL;

    for (i = 0; i < sys->numPorts; i++){
        for (j = 0; j < sys->numPorts; j++){
            //cout << sys->Y[i*sys->numPorts + j].real() << " " << sys->Y[i*sys->numPorts + j].imag() << endl;
            cout << sys->Y[i*sys->numPorts + j] << endl;
        }
    }
    //cout << "leng_v0c: " << leng_v0c << endl;

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
    cout << itercount << "\n";
    free(tmp);
    free(temp);
    free(ipar);
    free(dpar);
    free(y);

    return 0;
}



int matrixMulti(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> bRowId, vector<int> bColId, vector<double> bval, vector<int> &cRowId, vector<int> &cColId, vector<double> &cval){
    //the first matrix is row by row, the second matrix is column by column

    int i=0, j=0;
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


int mvMulti(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> &bRowId, vector<int> &bColId, vector<double> &bval, double *index_val, int size){
    //the same sequence in aColId and index
    double *v;
    int i;

    i = 0;
    v = (double*)calloc(size, sizeof(double));
    while (i < aColId.size()){
        v[aRowId[i]] += index_val[aColId[i]]*aval[i];
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


int nodeAdd(vector<int> &rowId, vector<int> &colId, vector<double> &val, vector<int> index, int size, fdtdMesh *sys){
    /*int i1, i2, j;
    double *v;

    rowId.clear();
    colId.clear();
    val.clear();
    for (j = 0; j < sys->nodeEdge[index[0]].size(); j++){
        rowId.push_back(sys->nodeEdge[index[0]][j].first);
        colId.push_back(0);
        val.push_back(sys->nodeEdge[index[0]][j].second);
    }
    for (j = 1; j < index.size(); j++){
        i1 = 0;
        i2 = 0;
        while (i1 < rowId.size() && i2 < sys->nodeEdge[index[j]].size()){
            if (rowId[i1] < sys->nodeEdge[index[j]][i2].first){
                i1++;
            }
            else if (rowId[i1] > sys->nodeEdge[index[j]][i2].first){
                rowId.insert(rowId.begin() + i1, sys->nodeEdge[index[j]][i2].first);
                colId.push_back(0);
                val.insert(val.begin() + i1, sys->nodeEdge[index[j]][i2].second);
                i1++;
                i2++;
            }
            else{
                val[i1] += sys->nodeEdge[index[j]][i2].second;
                if (abs(val[i1]) < 1e-10){    // if the value is 0, remove this element in the sparse form
                    rowId.erase(rowId.begin() + i1);
                    colId.pop_back();
                    val.erase(val.begin() + i1);
                    i1--;
                }
                i1++;
                i2++;
            }
        }
        while (i2 < sys->nodeEdge[index[j]].size()){
            rowId.push_back(sys->nodeEdge[index[j]][i2].first);
            colId.push_back(0);
            val.push_back(sys->nodeEdge[index[j]][i2].second);
            i2++;
        }
    }

    return 0;*/
    int i, j;
    double *v;

    v = (double*)calloc(size, sizeof(double));
    for (i = 0; i < index.size(); i++){
        for (j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            v[sys->nodeEdge[index[i]][j].first] += sys->nodeEdge[index[i]][j].second;
        }
    }
    for (i = 0; i < size; i++){
        if (abs(v[i]) > 1e-5) {
            rowId.push_back(i);
            colId.push_back(0);
            val.push_back(v[i]);
        }
    }
    free(v);
    v = NULL;
    return 0;

}

int nodeAddAvg(vector<int> &rowId, vector<int> &colId, vector<double> &val, int index, int size, fdtdMesh *sys){    // Get the average V0d2 (around the conductor)
    int i, j;
    double *v;
    int inx, iny, inz;

    v = (double*)calloc(size, sizeof(double));
    inz = index / sys->N_node_s;
    inx = (index - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    iny = index % (sys->N_cell_y + 1);
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

    rowId.clear();
    colId.clear();
    val.clear();


    int *visited;
    vector<int> st;
    double ratio;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    int mark;
    int count;
    st.clear();
    st.push_back(index);
    visited[index] = 1;
    while (!st.empty()){
        mark = 0;
        for (j = 0; j < sys->nodeEdge[st.back()].size(); j++){
            if (sys->markEdge[sys->nodeEdge[st.back()][j].first] != 0){
                if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] == 0)){
                    visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] = 1;
                    
                    inz = sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] / sys->N_node_s;
                    inx = (sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iny = sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] % (sys->N_cell_y + 1);
                    if (iny == 0){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                        colId.push_back(1);
                        val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else if (iny == sys->N_cell_y){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                        colId.push_back(1);
                        val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else{
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                        colId.push_back(1);
                        val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                        colId.push_back(1);
                        val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }

                    if (inx == 0){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                        colId.push_back(1);
                        val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else if (inx == sys->N_cell_x){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                        colId.push_back(1);
                        val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else{
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                        colId.push_back(1);
                        val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                        colId.push_back(1);
                        val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }

                    if (inz == 0){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else if (inz == sys->N_cell_z){
                        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                        if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else{
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                        if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    
                    for (i = 0; i < rowId.size(); i++){
                        v[rowId[i]] = v[rowId[i]] +  ratio * val[i];
                    }
                    rowId.clear();
                    colId.clear();
                    val.clear();

                    st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]);
                    mark = 1;

                    break;
                }
                else if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] == 0)){
                    visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] = 1;
                    
                    inz = sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] / sys->N_node_s;
                    inx = (sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iny = sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] % (sys->N_cell_y + 1);
                    if (iny == 0){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                        colId.push_back(1);
                        val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else if (iny == sys->N_cell_y){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                        colId.push_back(1);
                        val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else{
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                        colId.push_back(1);
                        val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                        colId.push_back(1);
                        val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }

                    if (inx == 0){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                        colId.push_back(1);
                        val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else if (inx == sys->N_cell_x){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                        colId.push_back(1);
                        val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else{
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                        colId.push_back(1);
                        val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                        colId.push_back(1);
                        val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }

                    if (inz == 0){
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else if (inz == sys->N_cell_z){
                        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                        if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    else{
                        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                        if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                        colId.push_back(1);
                        val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                        if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.back()][j].first){
                            ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.back()][j].first];
                        }
                    }
                    for (i = 0; i < rowId.size(); i++){
                        v[rowId[i]] += ratio * val[i];
                    }
                    rowId.clear();
                    colId.clear();
                    val.clear();

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

    for (i = 0; i < size; i++){
        if (abs(v[i]) > 1e-5) {
            rowId.push_back(i);
            colId.push_back(0);
            val.push_back(v[i]);
        }
    }

    free(visited);
    visited = NULL;
    free(v);
    v = NULL;

    return 0;
}
