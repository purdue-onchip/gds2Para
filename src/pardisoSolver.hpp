#ifndef GDS2PARA_PARDISO_SOLVER_H_
#define GDS2PARA_PARDISO_SOLVER_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

typedef vector< tuple<myint, myint, double> > BlockType;    // store rows-cols-vals of each block matrix

#define BUFF_ALIGN 64   // parameter of mkl_malloc, alignment of the buffer

// 0-based row-major 3-array CSR format of a matrix
class csrFormatOfMatrix {
public:
    // matrix information
    myint N_rows;
    myint N_cols;
    myint N_nnz;
    myint *rows         = nullptr;      // CSR row indices, size = N_rows + 1
    myint *cols         = nullptr;      // CSR col indices, size = N_nnz
    double *vals        = nullptr;      // CSR nnz values,  size = N_nnz

    // Constructor
    csrFormatOfMatrix(myint N_rows, myint N_cols, myint N_nnz) {
        /* Inputs: 
            N_rows*N_cols:  matrix size;
            N_nnz:          num of nonzeros.     */

        this->N_rows = N_rows;
        this->N_cols = N_cols;
        this->N_nnz = N_nnz;

        this->rows = (myint *)mkl_malloc(sizeof(myint) * (N_rows + 1), BUFF_ALIGN);
        this->cols = (myint *)mkl_malloc(sizeof(myint) * N_nnz, BUFF_ALIGN);
        this->vals = (double *)mkl_malloc(sizeof(double) * N_nnz, BUFF_ALIGN);
    }

    // Destructor
    ~csrFormatOfMatrix() {
        mkl_free(this->rows);
        mkl_free(this->cols);
        mkl_free(this->vals);
    }

    // Declaration of functions
    int convertBlockTypeToCsr(const BlockType &block);
};

int csrFormatOfMatrix::convertBlockTypeToCsr(const BlockType &block) {
    /* This function convert BlockType (COO) to CSR format. 
    The input BlockType must have been sorted by row index. */

    // The first index of rows is always 0 and the last is always N_nnz
    this->rows[0] = 0;
    this->rows[this->N_rows] = block.size();

    myint nnz_ind = 0;
    myint thisRow_ind = 0;
    for (const auto &nnz_tuple : block) {
        
        // cols and vals are direct copy of the tuple for every nnz
        this->cols[nnz_ind] = get<1>(nnz_tuple);
        this->vals[nnz_ind] = get<2>(nnz_tuple);

        // When saw nnz at a new row
        while (thisRow_ind < get<0>(nnz_tuple)) {
            thisRow_ind++;
            this->rows[thisRow_ind] = nnz_ind;
        }

        nnz_ind++;

        // Extrame case: when the last few rows are all blank
        while (nnz_ind == block.size() && thisRow_ind < this->N_rows-1) {
            thisRow_ind++;
            this->rows[thisRow_ind] = nnz_ind;
        }
    }

    return 0;
}

// 0-based column-major dense format of a matrix stored col by col in 1-D array
class denseFormatOfMatrix {
public:
    // matrix information
    myint N_rows;
    myint N_cols;
    myint matrixSize;
    double *vals = nullptr;

    // Constructor
    denseFormatOfMatrix(myint N_rows, myint N_cols) {
        this->N_rows = N_rows;
        this->N_cols = N_cols;
        this->matrixSize = N_rows * N_cols;

        this->vals = (double *)mkl_malloc(sizeof(double) * this->matrixSize, BUFF_ALIGN);
    }

    // Destructor
    ~denseFormatOfMatrix() {
        mkl_free(this->vals);
    }

    // Convert BlockType (COO) to dense format.
    void convertBlockTypeToDense(const BlockType &block) {
        for (const auto &nnz_tuple : block) {
            myint row_ind = get<0>(nnz_tuple);
            myint col_ind = get<1>(nnz_tuple);
            myint ind_inArray = col_ind * this->N_rows + row_ind;
            this->vals[ind_inArray] = get<2>(nnz_tuple);
        }
    }
};

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
inline myint mapRowColIdToBlockId(const myint B_rowId, const myint B_colId) {
    /* Explanation to notations here:

    Matrix S is partitioned into blocks as surface-vol-surface-... along row/col
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

// Solve D0sD1s = inv(B22)*B21B23 in Pardiso
int solveInverseInPardiso(const csrFormatOfMatrix &csrB22, 
    const denseFormatOfMatrix &denseB21B23, denseFormatOfMatrix *pdenseD0sD1s) {

    // Pardiso parameters, see https://software.intel.com/en-us/mkl-developer-reference-c-pardiso
    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 11;                 /* Real unsymmetric matrix */
    MKL_INT phase = 13;
    MKL_INT perm;
    myint nrhs = denseB21B23.N_cols;    /* Number of right hand sides */
    MKL_INT msglvl = 0;                 /* If msglvl=1, print statistical information */
    MKL_INT error = 0;

    void *pt[64];
    myint iparm[64];
    pardisoinit(pt, &mtype, iparm);
    iparm[38] = 1;
    iparm[34] = 1;         /* 0-based indexing */
    iparm[3] = 0;          /* No iterative-direct algorithm */
    //iparm[59] = 2;       /* out of core version to solve very large problem */
    //iparm[10] = 0;       /* Use nonsymmetric permutation and scaling MPS */

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &(csrB22.N_rows), csrB22.vals, csrB22.rows, csrB22.cols,
        &perm, &nrhs, iparm, &msglvl, denseB21B23.vals, pdenseD0sD1s->vals, &error);
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }

    // Release internal memory
    phase = -1;
    double ddum;           /* Double dummy */
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &(csrB22.N_rows), &ddum, csrB22.rows, csrB22.cols,
        &perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    return 0;
}

// Compute y = alpha * A * x + beta * y and store computed dense matrix in y
int csrMultiplyDense(const sparse_matrix_t &csrA_mklHandle,     // mkl handle of csr A
                        myint N_rows_x, double *x,              // x,y ~ dense matrix stored in array 
                        denseFormatOfMatrix *y) {
    /* see doc: https://software.intel.com/en-us/onemkl-developer-reference-c-mkl-sparse-mm
    Example: 
            y   = alpha * A   * x   + beta * y
            C11 = -1.0  * B12 * D0s + 1.0  * B11    */


    double alpha = -1.0, beta = 1.0;
    struct matrix_descr descrA;             // Descriptor of main sparse matrix properties
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Compute y = alpha * A * x + beta * y
    mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE,
        alpha,
        csrA_mklHandle,
        descrA,
        SPARSE_LAYOUT_COLUMN_MAJOR,         // Describes the storage scheme for the dense matrix
        x,
        y->N_cols,
        N_rows_x,
        beta,
        y->vals,
        y->N_rows);

    return 0;
}

void eliminateVolumE(const vector<BlockType> &layerS, myint N_surfE, myint N_volE) {
    /* This function eliminates e_vol from whole S matrix (all layers coupled) and obtain the
    2*2 block matrix related only to e_surf at single layer.

    Inputs:
        layerS: a vector containing 9 blocks of the whole matrix S, ~ within one layer
        N_surfE: number of {e}_surface at one surface, e.g. 0s
        N_volE: number of {e}_volumn at one layer, e.g. 0v
    Output:
        reducedS: 

    Each isolated layer here contains 2 surfaces and 1 middle volumn e, namely 0s-0v-1s. 
    A symbolic form is (each number represents the blockId in vector<BlockType>):
                
                 0s  0v  1s                                    0s  1s
    layerS =   | 0   1   2 |    0s      ==>     reducedS =   | 0'  1'|    0s
               | 3   4   5 |    0v                           | 2'  3'|    1s
               | 6   7   8 |    1s
    
    corresponding to block variable names in this function:
                 0s  0v  1s                                    0s  1s
    layerS ~   | B11 B12 B13 |  0s      ==>     reducedS ~   | C11 C12 |  0s
               | B21 B22 B23 |  0v                           | C21 C22 |  1s
               | B31 B32 B33 |  1s
    
    Due to the structure of matrix S, the central row of layerS is full, which could 
    be used to eliminate e0v as:  
        B21*e0s+B22*e0v+B23*e1s = 0  =>  e0v = -D0s*e0s - D1s*e1s, 
        where 
            D0s = inv(B22)*B21, D1s = inv(B22)*B23.
        Then the relation is: 
            C11 = B11 - B12*D0s,    C12 = B13 - B12*D1s
            C21 = B31 - B32*D0s,    C22 = B33 - B32*D1s
    */
    
    // Combine B21 and B23 in order to be sloved in one Pardiso run.
    BlockType vB21B23(layerS[3]);
    move(layerS[5].begin(), layerS[5].end(), back_inserter(vB21B23));
    denseFormatOfMatrix denseB21B23(N_volE, 2 * N_surfE);       // combined dense [B21, B23]
    denseB21B23.convertBlockTypeToDense(vB21B23);
    vB21B23.clear();

    // Solve D0s = inv(B22)*B21, D1s = inv(B22)*B23 in Pardiso
    csrFormatOfMatrix csrB22(N_volE, N_volE, layerS[4].size());
    csrB22.convertBlockTypeToCsr(layerS[4]);
    denseFormatOfMatrix denseD0sD1s(N_volE, 2 * N_surfE);       // combined dense [D0s, D1s]

    solveInverseInPardiso(csrB22, denseB21B23, &denseD0sD1s);
    csrB22.~csrFormatOfMatrix();                                // free CSR B22
    denseB21B23.~denseFormatOfMatrix();                         // free combined dense [B21, B23]

    // Convert B12, B32 to mkl internal CSR matrix handles
    sparse_matrix_t csrB12_mklHandle, csrB32_mklHandle;

    csrFormatOfMatrix csrB12(N_surfE, N_volE, layerS[1].size());
    csrB12.convertBlockTypeToCsr(layerS[1]);
    mkl_sparse_d_create_csr(&csrB12_mklHandle, SPARSE_INDEX_BASE_ZERO, csrB12.N_rows, csrB12.N_cols,
        csrB12.rows, csrB12.rows + 1, csrB12.cols, csrB12.vals);
    csrB12.~csrFormatOfMatrix();                                // free CSR B12

    csrFormatOfMatrix csrB32(N_surfE, N_volE, layerS[7].size());
    csrB32.convertBlockTypeToCsr(layerS[7]);
    mkl_sparse_d_create_csr(&csrB32_mklHandle, SPARSE_INDEX_BASE_ZERO, csrB32.N_rows, csrB32.N_cols,
        csrB32.rows, csrB32.rows + 1, csrB32.cols, csrB32.vals);
    csrB32.~csrFormatOfMatrix();                                // free CSR B32

    // Solve reducedS blocks C11 = B11 - B12*D0s, C12 = B13 - B12*D1s
    denseFormatOfMatrix denseB11(N_surfE, N_surfE);
    denseB11.convertBlockTypeToDense(layerS[0]);
    csrMultiplyDense(csrB12_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals, &denseB11);

    myint matrixSizeD0s = N_volE * N_surfE;
    denseFormatOfMatrix denseB13(N_surfE, N_surfE);
    denseB13.convertBlockTypeToDense(layerS[2]);
    csrMultiplyDense(csrB12_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals + matrixSizeD0s, &denseB13);

    mkl_sparse_destroy(csrB12_mklHandle);                       // free mkl CSR B12 handle

    // Solve reducedS blocks C21 = B31 - B32*D0s, C22 = B33 - B32*D1s
    denseFormatOfMatrix denseB31(N_surfE, N_surfE);
    denseB31.convertBlockTypeToDense(layerS[6]);
    csrMultiplyDense(csrB32_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals, &denseB31);

    denseFormatOfMatrix denseB33(N_surfE, N_surfE);
    denseB33.convertBlockTypeToDense(layerS[8]);
    csrMultiplyDense(csrB32_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals + matrixSizeD0s, &denseB33);

    mkl_sparse_destroy(csrB32_mklHandle);                       // free mkl CSR B32 handle

    denseD0sD1s.~denseFormatOfMatrix();                         // free combined dense [D0s, D1s]
}

void cascadeMatrixS(fdtdMesh *psys) {

    // Num of e at each surface or each layer. Layer growth along y.
    myint Nx                    = psys->N_cell_x;
    myint Nz                    = psys->N_cell_z;
    myint n_surfExEz            = Nx*(Nz + 1) + Nz*(Nx + 1);
    myint n_volEy               = (Nx + 1)*(Nz + 1);
    myint n_layerE_growY        = n_surfExEz + n_volEy;
    myint N_layers              = psys->N_cell_y;   // num of layers

    // Store all Bolck matrices in a vector (2-D vector)
    vector< BlockType > Blocks(8 * N_layers + 1);

    // Determine the block id of each nnz element of S and store in corresponding block matrix
    myint nnzS_rowId, nnzS_colId, B_rowId, B_colId, BlockId;
    for (myint i_nnz = 0; i_nnz < psys->leng_S; i_nnz++) {
        nnzS_rowId = psys->SRowId[i_nnz];
        nnzS_colId = psys->SColId[i_nnz];

        B_rowId = nnzS_rowId / n_layerE_growY * 2 + (nnzS_rowId % n_layerE_growY) / n_surfExEz;
        B_colId = nnzS_colId / n_layerE_growY * 2 + (nnzS_colId % n_layerE_growY) / n_surfExEz;

        BlockId = mapRowColIdToBlockId(B_rowId, B_colId);

        // Shift the start row index and col index to be 0 inside each block
        nnzS_rowId = nnzS_rowId - (B_rowId / 2) * n_layerE_growY - (B_rowId % 2) * n_surfExEz;
        nnzS_colId = nnzS_colId - (B_colId / 2) * n_layerE_growY - (B_colId % 2) * n_surfExEz;

        Blocks[BlockId].push_back(make_tuple(nnzS_rowId, nnzS_colId, psys->Sval[i_nnz]));
    }

    // Free original matrix S to save memory
    free(psys->SRowId);
    free(psys->SColId);
    free(psys->Sval);

    /******************** Start Cascading Matrix S from Blocks ***************************/
    
    // Consider PEC BCs on top or bottom (zmin & zmax) here
    // To be done!

    // Tell the size of sparse block matrices 
    myint N_surfE = n_surfExEz;
    myint N_volE = n_volEy;

    // Sort each block matrix by its row indices, 1st element of the tuple
    /* The purpose is to allow easy convert from COO to CSR format*/
    for (auto &block : Blocks) {
        sort(block.begin(), block.end());
    }

    // Half the value of overlapped ns-ns blocks between adjcent two layers
    /* The partitioned blocks no longer mean physical curl-curl opeartor, but mamatically, this is doable to
    cascade block matrix S. Any partition works like C = aC' + bC''. C = C' + C' makes codeing easier.*/
    for (myint nsns_BlockId = 8; nsns_BlockId < 8 * N_layers; nsns_BlockId += 8) {
        for (auto &nnz : Blocks[nsns_BlockId]) {
            get<2>(nnz) *= 0.5;
        }
    }

    // Eliminate e_volumn at each layer
    for (myint i_layer = 0; i_layer < N_layers; i_layer++) {
        
        vector<BlockType> layerS(Blocks.begin() + 8* i_layer, Blocks.begin() + 8 * i_layer + 9);
        eliminateVolumE(layerS, N_surfE, N_volE);

        layerS.clear();
    }

    // free blocks of original whole matrix S to save memory
    Blocks.clear();

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

    cascadeMatrixS(psys);

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