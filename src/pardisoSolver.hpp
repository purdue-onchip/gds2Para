#ifndef GDS2PARA_PARDISO_SOLVER_H_
#define GDS2PARA_PARDISO_SOLVER_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

using namespace std;

typedef vector< tuple<myint, myint, double> > BlockType;    // store rows-cols-vals of each block matrix

vector<myint> Map_eInd_GrowZ2Y(const myint Nx, const myint Ny, const myint Nz) {
    /* This function changes the global {e} index from layer growth along Z to that along Y.

    The index follows y - x - z ordering, and {e} is stacked layer by layer as
    {e_surface, e_volumn}.T at each layer. For the two cases here :
    case 1 ~ grow along Z: {e} = {e_s, | e_v}.T = {ey,ex, | ez}.T, in which vector
    {ey} first runs through y for each x, then runs through x as by
    y - x - z ordering. Same for {ex} ~y->x and {ez} ~y->x.
    case 2 ~ grow along Y: {e} = {ex,ez, | ey}.T at every layer.
    {ex}, {ez} and {ey} ~x->z, frist all x for each z then all z

    Input: num of bricks Nx*Ny*Nz
    Output: vector "eInd_map_z2y" of size N_e mapping global indices of {e}, eInd_map_z2y[oldInd] = newInd

    Yee's grid is used, with E at edge center and H at face center. Outmost
    boundaries are all E edges. Removal of {e} due to PEC BC has not been considered. */

    // Num of surface or volumetric unknown e at each layer
    myint n_surfEy_growZ        = Ny*(Nx + 1);
    myint n_surfEyEx            = n_surfEy_growZ + Nx*(Ny + 1);
    myint n_volEz               = (Nx + 1)*(Ny + 1);
    myint n_layerE_growZ        = n_surfEyEx + n_volEz;

    myint n_surfEx_growY        = Nx*(Nz + 1);
    myint n_surfExEz            = n_surfEx_growY + Nz*(Nx + 1);
    myint n_volEy               = (Nx + 1)*(Nz + 1);
    myint n_layerE_growY        = n_surfExEz + n_volEy;

    myint N_tot_E               = n_surfEyEx*(Nz + 1) + n_volEz*Nz;


    // y - x - z ordering
    // Grow along Z : {e} = { ey,ex, | ez }.T, y->x frist all y for each x then all x
    // Grow along Y : {e} = { ex,ez, | ey }.T, x->z frist all x for each z then all z

    vector<myint> eInd_map_z2y(N_tot_E, 0);

    // Map the edges along y direction, Ey
    for (myint iz = 0; iz < Nz + 1; iz++) {
        for (myint iy = 0; iy < Ny; iy++) {
            for (myint ix = 0; ix < Nx + 1; ix++) {
                eInd_map_z2y[iz*n_layerE_growZ + ix*(Ny)+iy] =
                    iy*n_layerE_growY + n_surfExEz + iz*(Nx + 1) + ix;
            }
        }
    }

    // Map the edges along x direction, Ex
    for (myint iz = 0; iz < Nz + 1; iz++) {
        for (myint iy = 0; iy < Ny + 1; iy++) {
            for (myint ix = 0; ix < Nx; ix++) {
                eInd_map_z2y[iz*n_layerE_growZ + n_surfEy_growZ + ix*(Ny + 1) + iy] =
                    iy*n_layerE_growY + iz*(Nx)+ix;
            }
        }
    }

    // Map the edges along z direction, Ez
    for (myint iz = 0; iz < Nz; iz++) {
        for (myint iy = 0; iy < Ny + 1; iy++) {
            for (myint ix = 0; ix < Nx + 1; ix++) {
                eInd_map_z2y[iz*n_layerE_growZ + n_surfEyEx + ix*(Ny + 1) + iy] =
                    iy*n_layerE_growY + n_surfEx_growY + iz*(Nx + 1) + ix;
            }
        }
    }

    return eInd_map_z2y;
}

// Reverse the map by swapping the index and value of the mapping array
vector<myint> Reverse_Map_eInd(const vector<myint> &eInd_map_z2y) {
    
    myint N_tot_E = eInd_map_z2y.size();
    vector<myint> eInd_map_y2z(N_tot_E, 0);

    for (myint id = 0; id < N_tot_E; id++) {
        eInd_map_y2z[eInd_map_z2y[id]] = id;
    }

    return eInd_map_y2z;
}

// Map specific Block_rowId and Block_colId to an unique Block index
inline myint Map_RowCol_2BlockId(const myint B_rowId, const myint B_colId) {
    /* Explanation to notations here:

    Matrix S is partitioned into blocks as surface-vlo-surface-... along row/col
    Bolck ordering:     0s,0v,1s,1v,..., Nlayer_s
    Block row/col Id:   0, 1, 2, 3, ..., 2*Nlayer
    BlockId:            counted nnz block Id, row major. 

    Example: (Nlayer = 3, Nblock = 8*Nlayer + 1 = 25)
                0s  0v  1s  1v  2s  2v  3s
    B_colId->   0   1   2   3   4   5   6
                                                B_rowId (below) 
    BlockId:  | 0   1   2                  |    0   0s
              | 3   4   5                  |    1   0v
              | 6   7   8   9   10         |    2   1s
              |         11  12  13         |    3   1v
              |         14  15  16  17  18 |    4   2s
              |                 19  20  21 |    5   2v
              |                 22  23  24 |    6   3s

    Input: Block row/col Id
    Output: BlockId of this block    */
    
    myint BlockId = 0;
    myint layerId = B_rowId / 2;
    myint isVol = B_rowId % 2;

    if (layerId == 0) {     // B_rowId ~ 0s, 0v
        BlockId = B_rowId * 3 + B_colId;
    }
    else if(isVol == 1) {   // The B_row is n-v
        BlockId = 6 * layerId + B_colId + 3;
    }
    else {                  // The B_row is n-s
        BlockId = 6 * layerId + B_colId;
    }
    
    return BlockId;
}

void EliminateVolumE(const vector<BlockType> &layerS) {
    /* This function eliminates e_vol from whole S matrix (all layers coupled) and obtain the
    2*2 block matrix regarding the coupled e_surf at an isolated layer.

    Input:
        layerS: a vector containing 9 blocks of the whole matrix S, ~ within one layer
    Output:
        reducedS: 

    Each isolated layer here contains 2 surfaces and 1 middle volumn e, namely 0s-0v-1s. 
    A symbolic form is (each number represents a block matrix):
                
                 0s  0v  1s                                    0s  1s
    layerS =   | 0   1   2 |    0s      ==>     reducedS =   | 0'  1'|    0s
               | 3   4   5 |    0v                           | 2'  3'|    1s
               | 6   7   8 |    1s
    Notes:
    1)  9 blocks of the input layerS is directly truncated from system matrix S, so block 0 and block 9 
        might be double counted compared to isolated single-layer S.
    2)  PEC BCs on top or bottom (z-direction) will be considered here
    */
    
}

void Cascade_matrixS(fdtdMesh *psys) {

    // Num of e at each surface or each layer. Layer growth along y.
    myint Nx                    = psys->N_cell_x;
    myint Nz                    = psys->N_cell_z;
    myint n_surfExEz            = Nx*(Nz + 1) + Nz*(Nx + 1);
    myint n_volEy               = (Nx + 1)*(Nz + 1);
    myint n_layerE_growY        = n_surfExEz + n_volEy;

    // Store all Bolck matrices in a vector (2-D vector)
    vector< BlockType > Blocks(8* psys->N_cell_y + 1);

    // Determine the block id of each nnz element of S and store in corresponding block matrix
    myint nnzS_rowId, nnzS_colId, B_rowId, B_colId, BlockId;
    for (myint i_nnz = 0; i_nnz < psys->leng_S; i_nnz++) {
        nnzS_rowId = psys->SRowId[i_nnz];
        nnzS_colId = psys->SColId[i_nnz];

        B_rowId = nnzS_rowId / n_layerE_growY * 2 + (nnzS_rowId % n_layerE_growY) / n_surfExEz;
        B_colId = nnzS_colId / n_layerE_growY * 2 + (nnzS_colId % n_layerE_growY) / n_surfExEz;

        BlockId = Map_RowCol_2BlockId(B_rowId, B_colId);

        Blocks[BlockId].push_back(make_tuple(nnzS_rowId, nnzS_colId, psys->Sval[i_nnz]));
    }

    // Free original matrix S to save memory
    free(psys->SRowId);
    free(psys->SColId);
    free(psys->Sval);

    // Eliminate e_volumn at each layer
    myint N_layers = psys->N_cell_y;
    for (myint i_layer = 0; i_layer < N_layers; i_layer++) {
        
        vector< BlockType > layerS(Blocks.begin() + 8* i_layer, Blocks.begin() + 8 * i_layer + 9);
        EliminateVolumE(layerS);

        layerS.clear();
    }
    Blocks.clear();     // free blocks in matrix S to save memory

}

// Cal all the computed freq points and store in a vector
vector<double> CalAllFreqPointsHz(const fdtdMesh &sys) {
    vector<double> vFreqHz;

    vFreqHz.push_back(sys.freqStart * sys.freqUnit);          // First frequency in sweep

    for (int id = 1; id < sys.nfreq; id++) {                  // When nfreq > 1
        if (sys.freqScale == 1) {
            vFreqHz.push_back((sys.freqStart + id * (sys.freqEnd - sys.freqStart) / (sys.nfreq - 1)) * sys.freqUnit);
        }   // Linear interpolation of frequency sweep
        else {
            vFreqHz.push_back(sys.freqStart * sys.freqUnit * pow(sys.freqEnd / sys.freqStart, (id * 1.0 / (sys.nfreq - 1))));
        }   // Logarithmic interpolation of frequency sweep
    }

    return vFreqHz;
}

// Assign source current density for a source port
void Assign_J_eachPort(fdtdMesh *psys, int sourcePort) {
    if (psys->J != nullptr) {
        free(psys->J);
    }   // delete previous J assignment

    psys->J = (double*)calloc(psys->N_edge, sizeof(double));
    for (int sourcePortSide = 0; sourcePortSide < psys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
        for (int indEdge = 0; indEdge < psys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
            psys->J[psys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]] = psys->portCoor[sourcePort].portDirection[sourcePortSide];
        }
    }
}

// Solve E field and Z-parameters in Pardiso, solve layer by layer. (under developing)
void Solve_E_Zpara_InPardiso_layered(fdtdMesh *psys) {

    Cascade_matrixS(psys);

    // All computed freq points
    vector<double> vFreqHz = CalAllFreqPointsHz(*psys);

    // S matrix in CSR format
    const myint *pcsrS_colId = psys->SColId;

    // Right-Hand-Side term -iwJ (A * m^-2 / s) and electric field eField (V/m), SI unit
    complex<double> *pRhsJ_SI = new complex<double>[psys->N_edge - 2 * psys->N_edge_s]();    // -iw{j}
    complex<double> *peField_SI = new complex<double>[psys->N_edge - 2 * psys->N_edge_s]();  // {e}
    
    // Use complex double constructor to assign initial output matrix for single-frequency solve
    psys->x.assign(psys->numPorts * psys->numPorts * psys->nfreq, complex<double>(0., 0.));

    for (int indFreq = 0; indFreq < vFreqHz.size(); indFreq++) {
        //for (int sourcePort = 0; sourcePort < psys->numPorts; sourcePort++) {
        //    AssignSourceCurrentForSourcePort(psys, sourcePort);
        //    Solve_E_InPardiso(psys, peField_SI, pRhsJ_SI, vFreqHz[indFreq], psys->SRowId, psys->SColId, psys->Sval);
        //    psys->Construct_Z_V0_Vh(peField_SI, indFreq, sourcePort, bdl, bdu);
        //}   // solve for each port
        reference(psys, indFreq, psys->SRowId, psys->SColId, psys->Sval);
        
    }   // solve for each frequency

    // Print Z-parameters
    psys->print_z_V0_Vh();

    delete[] pRhsJ_SI, peField_SI;
}

// Solve E field and Z-parameters in Pardiso, solve the whole structure as reference
void Solve_E_Zpara_InPardiso_reference(fdtdMesh *psys) {

    // All computed freq points
    vector<double> vFreqHz = CalAllFreqPointsHz(*psys);

    // Use complex double constructor to assign initial output matrix for single-frequency solve
    psys->x.assign(psys->numPorts * psys->numPorts * psys->nfreq, complex<double>(0., 0.));

    // At each frequency, solve all ports together
    for (int indFreq = 0; indFreq < vFreqHz.size(); indFreq++) {
        reference(psys, indFreq, psys->SRowId, psys->SColId, psys->Sval);
    }   

    // Print Z-parameters
    psys->print_z_V0_Vh();

}

#endif