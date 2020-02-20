#ifndef GDS2PARA_MAP_INDEX_H_
#define GDS2PARA_MAP_INDEX_H_

#include "fdtd.hpp"

class mapIndex {
public:
    // num of bricks Nx*Ny*Nz
    myint Nx, Ny, Nz;
    myint N_totEdges;                   // number of all {e} or edges without PEC removal
    myint N_totEdges_rmPEC;             // number of all {e} or edges after PEC removal

    // num of edges after removing PEC
    myint N_edgesAtPEC;                 // number of {e} or edges at PEC surfaces
    myint N_surfExEz_rmPEC;             // number of {e}_surface at one surface e.g. 0s, mode (growY, removed PEC)
    myint N_volEy_rmPEC;                // number of {e}_volume at one layer, e.g. 0v, mode (growY, removed PEC)

    // Upper PEC edge indices (z=zmax) and lower PEC edge indices (z=zmin), index mode (growY, no PEC removal)
    unordered_set<myint> edges_upperPEC;
    unordered_set<myint> edges_lowerPEC;

    // maps
    vector<myint> eInd_map_z2y;         // map {e} index: (growZ, no PEC removal) -> (growY, no PEC removal)
    vector<myint> eInd_map_y2z;         // map {e} index: (growY, no PEC removal) -> (growZ, no PEC removal)

    vector<myint> eInd_map_y2rmPEC;     // map {e} index: (growY, no PEC removal) -> (growY, removed PEC)
    vector<myint> eInd_map_rmPEC2y;     // map {e} index: (growY, removed PEC) -> (growY, no PEC removal)

    // sys.mapEdge;                     // map {e} index: (growZ, no PEC removal) -> (growZ, removed PEC)
    // sys.mapEdgeR;                    // map {e} index: (growZ, removed PEC) -> (growZ, no PEC removal)

    vector<myint> eInd_map_rmPEC_z2y;   // map {e} index: (growZ, removed PEC) -> (growY, removed PEC)
    vector<myint> eInd_map_rmPEC_y2z;   // map {e} index: (growY, removed PEC) -> (growZ, removed PEC)

    // Constructor
    mapIndex(myint Nx, myint Ny, myint Nz) {
        this->Nx = Nx;
        this->Ny = Ny;
        this->Nz = Nz;

        this->N_totEdges = Nx*(Ny + 1)*(Nz + 1) + (Nx + 1)*Ny*(Nz + 1) + (Nx + 1)*(Ny + 1)*Nz;

        // For different PEC BCs
        this->N_edgesAtPEC = 0;
        this->N_surfExEz_rmPEC = Nx*(Nz + 1) + Nz*(Nx + 1);
        this->N_volEy_rmPEC = (Nx + 1)*(Nz + 1);
#ifdef LOWER_BOUNDARY_PEC
        this->N_edgesAtPEC += Nx*(Ny + 1) + (Nx + 1)*Ny;    // add num of Ex & Ey at PEC surface
        this->N_surfExEz_rmPEC -= Nx;                       // remove Ex at PEC
        this->N_volEy_rmPEC -= (Nx + 1);                    // remove Ey at PEC
#endif
#ifdef UPPER_BOUNDARY_PEC
        this->N_edgesAtPEC += Nx*(Ny + 1) + (Nx + 1)*Ny;    // add num of Ex & Ey at PEC surface
        this->N_surfExEz_rmPEC -= Nx;                       // remove Ex at PEC
        this->N_volEy_rmPEC -= (Nx + 1);                    // remove Ey at PEC
#endif
        this->N_totEdges_rmPEC = this->N_totEdges - this->N_edgesAtPEC;
    }

    // Set global {e} index map between mode (growZ, no PEC removal) and mode (growY, no PEC removal)
    void setEdgeMap_growZgrowY() {
    /*  The index follows y - x - z ordering, and {e} is stacked layer by layer as
        {e_surface, e_volume}.T at each layer. For the two cases here:
        - case 1 ~ grow along Z: {e} = {e_s, | e_v}.T = {ey,ex, | ez}.T, in which vector
                {ey} first runs through y for each x, then runs through x as by
                y - x - z ordering. Same for {ex} ~y->x and {ez} ~y->x.
        - case 2 ~ grow along Y: {e} = {ex,ez, | ey}.T at every layer.
                {ex}, {ez} and {ey} ~x->z, frist all x for each z then all z

        Results: 
            - 2 index maps are stored in attributes "eInd_map_z2y" and "eInd_map_y2z"
            - Each map is a vector of size N_e, eInd_map_?2?[oldInd] = newInd

        Yee's grid is used, with E at edge center and H at face center. Outmost
        boundaries are all E edges. Removal of {e} due to PEC BC has not been considered. */

        // Num of surface or volumetric unknown e at each layer
        myint Nx = this->Nx;
        myint Ny = this->Ny;
        myint Nz = this->Nz;

        myint n_surfEy_growZ = Ny*(Nx + 1);
        myint n_surfEyEx = n_surfEy_growZ + Nx*(Ny + 1);
        myint n_volEz = (Nx + 1)*(Ny + 1);
        myint n_layerE_growZ = n_surfEyEx + n_volEz;

        myint n_surfEx_growY = Nx*(Nz + 1);
        myint n_surfExEz = n_surfEx_growY + Nz*(Nx + 1);
        myint n_volEy = (Nx + 1)*(Nz + 1);
        myint n_layerE_growY = n_surfExEz + n_volEy;

        myint N_tot_E = n_surfEyEx*(Nz + 1) + n_volEz*Nz;

        // y - x - z ordering
        // Grow along Z : {e} = { ey,ex, | ez }.T, y->x frist all y for each x then all x
        // Grow along Y : {e} = { ex,ez, | ey }.T, x->z frist all x for each z then all z

        // Map Grow Z to Grow Y
        this->eInd_map_z2y.reserve(N_tot_E);
        this->eInd_map_z2y.assign(N_tot_E, -1);

        for (myint iz = 0; iz < Nz + 1; iz++) {     // Map the edges along y direction, Ey
            for (myint iy = 0; iy < Ny; iy++) {
                for (myint ix = 0; ix < Nx + 1; ix++) {
                    this->eInd_map_z2y[iz*n_layerE_growZ + ix*(Ny)+iy] =
                        iy*n_layerE_growY + n_surfExEz + iz*(Nx + 1) + ix;
                }
            }
        }
        for (myint iz = 0; iz < Nz + 1; iz++) {     // Map the edges along x direction, Ex
            for (myint iy = 0; iy < Ny + 1; iy++) {
                for (myint ix = 0; ix < Nx; ix++) {
                    this->eInd_map_z2y[iz*n_layerE_growZ + n_surfEy_growZ + ix*(Ny + 1) + iy] =
                        iy*n_layerE_growY + iz*(Nx)+ix;
                }
            }
        }
        for (myint iz = 0; iz < Nz; iz++) {         // Map the edges along z direction, Ez
            for (myint iy = 0; iy < Ny + 1; iy++) {
                for (myint ix = 0; ix < Nx + 1; ix++) {
                    this->eInd_map_z2y[iz*n_layerE_growZ + n_surfEyEx + ix*(Ny + 1) + iy] =
                        iy*n_layerE_growY + n_surfEx_growY + iz*(Nx + 1) + ix;
                }
            }
        }

        // Map Grow Y to Grow Z
        this->eInd_map_y2z.reserve(N_tot_E);
        this->eInd_map_y2z.assign(N_tot_E, -1);
        for (myint id = 0; id < N_tot_E; id++) {
            this->eInd_map_y2z[this->eInd_map_z2y[id]] = id;
        }
    }

    // Set global {e} index map between mode (growY, no PEC removal) and (growY, removed PEC)
    void setEdgeMap_growYremovePEC(const set<myint> &ubde, const set<myint> &lbde, myint bden) {
    /*  Inputs:
            - ubde = sys.ubde: upper PEC edge indices (z=zmax) with index mode (growZ, no PEC removal)
            - lbde = sys.lbde: lower PEC edge indices (z=zmin) with index mode (growZ, no PEC removal)
            - bden = sys.bden: number of PEC edges that was counted at meshing stage
        Results:
            - PEC edge indices of mode (growY, no PEC removal) stored in attributes "edges_upperPEC" & "edges_lowerPEC"
            - 2 index maps are stored in attributes "eInd_map_y2rmPEC" and "eInd_map_rmPEC2y"       */

        if (this->eInd_map_z2y.size() == 0) {
            cout << "Failure! Run function setEdgeMap_growZgrowY() first." << endl;
            exit(2);
        }
        if (this->N_edgesAtPEC != bden) {
            cout << "Failure at mapping PEC removal! Inconsistent PEC edges." << endl <<
                "Try to comment out code: sys->findBoundNodeEdge(inx, iny, inz); at mesh.cpp " << endl;
            exit(2);
        }

        // Map the PEC edge indices from (growZ, no PEC removal) to (growY, no PEC removal)
        for (auto eInd_growZ : ubde) {
            myint eInd_growY = this->eInd_map_z2y[eInd_growZ];
            this->edges_upperPEC.insert(eInd_growY);
        }
        for (auto eInd_growZ : lbde) {
            myint eInd_growY = this->eInd_map_z2y[eInd_growZ];
            this->edges_lowerPEC.insert(eInd_growY);
        }

        // Map mode (growY, no PEC removal) and (growY, removed PEC)
        this->eInd_map_y2rmPEC.reserve(this->N_totEdges);
        this->eInd_map_y2rmPEC.assign(this->N_totEdges, -1);        // removed PEC edges will be mapped to index -1
        this->eInd_map_rmPEC2y.reserve(this->N_totEdges - this->N_edgesAtPEC);
        this->eInd_map_rmPEC2y.assign(this->N_totEdges - this->N_edgesAtPEC, -1);
        
        myint count_edgesAtPEC = 0;                                 // number of identified PEC edges
        for (myint eInd = 0; eInd < this->N_totEdges; eInd++) {     // from smallest edge index in mode (growY, no PEC removal)
            bool edgeIsAtUpperPEC = this->edges_upperPEC.find(eInd) != this->edges_upperPEC.end();
            bool edgeIsAtLowerPEC = this->edges_lowerPEC.find(eInd) != this->edges_lowerPEC.end();
            if (edgeIsAtUpperPEC || edgeIsAtLowerPEC) {             // if this edge is at upper or lower PEC boundary
                count_edgesAtPEC++;
            }
            else {
                this->eInd_map_y2rmPEC[eInd] = eInd - count_edgesAtPEC;
                this->eInd_map_rmPEC2y[eInd - count_edgesAtPEC] = eInd;
            }
        }

    }

    // Set global {e} index map between mode (growZ, removed PEC) and (growY, removed PEC)
    void setEdgeMap_rmPEC_growZgrowY(myint *eInd_map_rmPEC2z) {
    /*  Inputs:
            - eInd_map_rmPEC2z = sys.mapEdgeR: map {e} index from (growZ, removed PEC) to (growZ, no PEC removal)
        Results:
            - 2 index maps are stored in attributes "eInd_map_rmPEC_z2y" and "eInd_map_rmPEC_y2z"       */

        if (this->eInd_map_y2rmPEC.size() == 0) {
            cout << "Failure! Run function setEdgeMap_growYremovePEC first." << endl;
            exit(2);
        }

        // Map mode (growZ, removed PEC) and (growY, removed PEC)
        myint N_edgesNotAtPEC = this->N_totEdges - this->N_edgesAtPEC;
        this->eInd_map_rmPEC_z2y.reserve(N_edgesNotAtPEC);
        this->eInd_map_rmPEC_z2y.assign(N_edgesNotAtPEC, -1);
        this->eInd_map_rmPEC_y2z.reserve(N_edgesNotAtPEC);
        this->eInd_map_rmPEC_y2z.assign(N_edgesNotAtPEC, -1);

        for (myint eInd_rmPECz = 0; eInd_rmPECz < N_edgesNotAtPEC; eInd_rmPECz++) { // (growZ, removed PEC)
            myint eInd_z = eInd_map_rmPEC2z[eInd_rmPECz];                           // -> (growZ, no PEC removal)
            myint eInd_y = this->eInd_map_z2y[eInd_z];                              // -> (growY, no PEC removal)
            myint eInd_rmPECy = this->eInd_map_y2rmPEC[eInd_y];                     // -> (growY, removed PEC)
            this->eInd_map_rmPEC_z2y[eInd_rmPECz] = eInd_rmPECy;                    // (growZ, removed PEC) -> (growY, removed PEC)
            this->eInd_map_rmPEC_y2z[eInd_rmPECy] = eInd_rmPECz;                    // (growY, removed PEC) -> (growZ, removed PEC)
        }
    }

    // Map a (Block_rowId, Block_colId) pair to one unique Block index
    myint mapBlockRowColToBlockInd(myint B_rowId, myint B_colId) const {
    /*  Explanation to notations here:

        Matrix S is partitioned into blocks as surface-vol-surface-... along row (col)
            - Bolck ordering:     0s,0v,1s,1v,..., Nlayer_s
            - B_rowId or B_colId: 0, 1, 2, 3, ..., 2*Nlayer
            - Retrun BlockId:     counted nnz block Id, row major.

        Example: (Nlayer = 3, Nblock = 8*Nlayer + 1 = 25)
                          0s  0v  1s  1v  2s  2v  3s
            B_colId->     0   1   2   3   4   5   6
                                                        B_rowId (below)
            BlockId:    | 0   1   2                  |    0   0s
                        | 3   4   5                  |    1   0v
                        | 6   7   8   9   10         |    2   1s
                        |         11  12  13         |    3   1v
                        |         14  15  16  17  18 |    4   2s
                        |                 19  20  21 |    5   2v
                        |                 22  23  24 |    6   3s        */

        myint BlockId = 0;
        myint layerId = B_rowId / 2;
        myint isVol = B_rowId % 2;

        if (layerId == 0) {     // B_rowId ~ 0s, 0v
            BlockId = B_rowId * 3 + B_colId;
        }
        else if (isVol == 1) {  // The B_row is n-v
            BlockId = 6 * layerId + B_colId + 3;
        }
        else {                  // The B_row is n-s
            BlockId = 6 * layerId + B_colId;
        }

        return BlockId;
    }
};



#endif