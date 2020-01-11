#ifndef GDS2PARA_PARDISO_SOLVER_H_
#define GDS2PARA_PARDISO_SOLVER_H_

#include <iostream>
#include <fstream>
#include <string>

#include "fdtd.hpp"
#include "matrixTypeDef.hpp"
#include "mapIndex.hpp"
//#define DEBUG_SOLVE_REORDERED_S   // debug mode: directly solve the entire reordered S (growZ, rmPEC)

// Compute y = alpha * A * x + beta * y and store computed dense matrix in y
int csrMultiplyDense(const sparse_matrix_t *csrA_mklHandle,     // mkl handle of csr A
    myint N_rows_x, complex<double> *xval,     // x,y ~ dense matrix stored in continuous memory array 
    denseFormatOfMatrix *y) {
    /* see doc: https://software.intel.com/en-us/onemkl-developer-reference-c-mkl-sparse-mm
    Example:
            y   = alpha * A   * x   + beta * y
            C11 = -1.0  * B12 * D0s + 1.0  * B11    */


    MKL_Complex16 alpha = { -1.0, 0.0 };
    MKL_Complex16 beta = { 1.0, 0.0 };
    struct matrix_descr descrA;             // Descriptor of main sparse matrix properties
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    vector<MKL_Complex16> yval_mklComplex(y->matrixSize);
    y->copyToMKL_Complex16(yval_mklComplex.data());
    myint x_matrixSize = N_rows_x * y->N_rows;
    vector<MKL_Complex16> xval_mklComplex(x_matrixSize);
    for (myint ind = 0; ind < x_matrixSize; ind++) {
        xval_mklComplex[ind].real = xval[ind].real();
        xval_mklComplex[ind].imag = xval[ind].imag();
    }

    // Compute y = alpha * A * x + beta * y
    sparse_status_t returnStatus = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,
        alpha,                              // alpha = -1.0
        *csrA_mklHandle,                    // A
        descrA,
        SPARSE_LAYOUT_COLUMN_MAJOR,         // Describes the storage scheme for the dense matrix
        xval_mklComplex.data(),
        y->N_cols,
        N_rows_x,
        beta,                               // beta = 1.0
        yval_mklComplex.data(),             // Pointer to the memory array used by the vector to store its owned elements
        y->N_rows);
    if (returnStatus != SPARSE_STATUS_SUCCESS) {
        cout << "ERROR! Return from mkl_sparse_z_mm is: " << returnStatus << endl;
        exit(2);
    }

    y->copyFromMKL_Complex16(yval_mklComplex.data());

    return 0;
}

int eliminateVolumE(const vector<BlockType> &layerS, myint N_surfE, myint N_volE, denseFormatOfMatrix *preducedS) {
    /* This function eliminates e_vol from whole S matrix (all layers coupled) and obtain the
    2*2 block matrix related only to e_surf at single layer.

    Inputs:
        layerS: a vector containing 9 blocks of the whole matrix S, ~ within one layer
        N_surfE: number of {e}_surface at one surface, e.g. 0s
        N_volE: number of {e}_volume at one layer, e.g. 0v
    Output:
        preducedS: pointing to memory of a vector storing 4 reduced dense blocks as {C11, C12, C21, C22}

    Each isolated layer here contains 2 surfaces and 1 middle volume e, namely 0s-0v-1s. 
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
    BlockType vB23_newColInd(layerS[5]);
    for (auto &nnz : vB23_newColInd) {
        nnz.col_ind += N_surfE;
    }   // shift the col index of COO B23 in order to combine [B21, B23] correctly
    move(vB23_newColInd.begin(), vB23_newColInd.end(), back_inserter(vB21B23));
    denseFormatOfMatrix denseB21B23(N_volE, 2 * N_surfE);       // combined dense [B21, B23]
    denseB21B23.convertBlockTypeToDense(vB21B23);
    vB21B23.clear();

    // Solve D0s = inv(B22)*B21, D1s = inv(B22)*B23 in Pardiso
    csrFormatOfMatrix csrB22(N_volE, N_volE, layerS[4].size()); // CSR B22
    csrB22.convertBlockTypeToCsr(layerS[4]);
    denseFormatOfMatrix denseD0sD1s = 
        csrB22.backslashDense(denseB21B23);                     // combined dense [D0s, D1s]
    
    /*denseB21B23.writeToFile("blockB21B23.txt");
    denseD0sD1s.writeToFile("blockD.txt");*/

    denseB21B23.~denseFormatOfMatrix();                         // free combined dense [B21, B23]

    // Convert B12, B32 to mkl internal CSR matrix handles
    sparse_matrix_t csrB12_mklHandle, csrB32_mklHandle;
    csrFormatOfMatrix csrB12(N_surfE, N_volE, layerS[1].size());
    csrB12.convertBlockTypeToCsr(layerS[1]);
    vector<MKL_Complex16> csrB12vals_mklCompl(csrB12.N_nnz);
    csrB12.copyToMKL_Complex16(csrB12vals_mklCompl.data());
    mkl_sparse_z_create_csr(&csrB12_mklHandle, SPARSE_INDEX_BASE_ZERO, csrB12.N_rows, csrB12.N_cols,
        csrB12.rows.data(), csrB12.rows.data() + 1, csrB12.cols.data(), csrB12vals_mklCompl.data());
    csrFormatOfMatrix csrB32(N_surfE, N_volE, layerS[7].size());
    csrB32.convertBlockTypeToCsr(layerS[7]);
    vector<MKL_Complex16> csrB32vals_mklCompl(csrB32.N_nnz);
    csrB32.copyToMKL_Complex16(csrB32vals_mklCompl.data());
    mkl_sparse_z_create_csr(&csrB32_mklHandle, SPARSE_INDEX_BASE_ZERO, csrB32.N_rows, csrB32.N_cols,
        csrB32.rows.data(), csrB32.rows.data() + 1, csrB32.cols.data(), csrB32vals_mklCompl.data());

    // Solve reducedS blocks C11 = B11 - B12*D0s, C12 = B13 - B12*D1s
    preducedS->convertBlockTypeToDense(layerS[0]);              // dense B11 & dense C11
    csrMultiplyDense(&csrB12_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals.data(), preducedS);
    myint matrixSizeD0s = N_volE * N_surfE;
    (preducedS + 1)->convertBlockTypeToDense(layerS[2]);        // dense B13 & dense C12
    csrMultiplyDense(&csrB12_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals.data() + matrixSizeD0s, preducedS + 1);
    mkl_sparse_destroy(csrB12_mklHandle);                       // free mkl CSR B12 handle

    // Solve reducedS blocks C21 = B31 - B32*D0s, C22 = B33 - B32*D1s
    (preducedS + 2)->convertBlockTypeToDense(layerS[6]);        // dense B31 & dense C21
    csrMultiplyDense(&csrB32_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals.data(), preducedS + 2);
    (preducedS + 3)->convertBlockTypeToDense(layerS[8]);        // dense B33 & dense C22
    csrMultiplyDense(&csrB32_mklHandle, denseD0sD1s.N_rows, denseD0sD1s.vals.data() + matrixSizeD0s, preducedS + 3);
    mkl_sparse_destroy(csrB32_mklHandle);                       // free mkl CSR B32 handle

    /*preducedS->writeToFile("block_C0.txt");
    (preducedS + 1)->writeToFile("block_C1.txt");
    (preducedS + 2)->writeToFile("block_C2.txt");
    (preducedS + 3)->writeToFile("block_C3.txt");*/

    denseD0sD1s.~denseFormatOfMatrix();                         // free combined dense [D0s, D1s]
    return 0;
}

// Reconstruct blocks stored in portportBlocks (or surfSurfBlocks) to S matrix in 1-D dense format
denseFormatOfMatrix reconstructBlocksToDense(const vector<vector<denseFormatOfMatrix>> &portportBlocks) {
    /* Example: N_layers = 2
            |D11     D12          |
            |D21   D22+D11'   D12'|     =   cascadedS
            |        D21'     D22'|
        Input: N_layers*4 2-D vector of blocks:  
            portportBlocks[0] = {D11, D12, D21, D22}
            portportBlocks[1] = {D11', D12', D21', D22'}
        Return:
            above block-tridiagonal square matrix cascadedS in 1-D dense format
    */

    // All blocks D and S are square matrix
    myint N_rows_perBlock = portportBlocks[0][0].N_rows;
    myint N_layers = portportBlocks.size();
    myint N_rows_S = N_rows_perBlock * (N_layers + 1);
    denseFormatOfMatrix cascadedS(N_rows_S, N_rows_S);
    denseFormatOfMatrix tempD(N_rows_perBlock, N_rows_perBlock);

    for (myint i_layer = 0; i_layer < N_layers; i_layer++) {
        for (myint j_col = 0; j_col < N_rows_perBlock; j_col++) {
            for (myint i_row = 0; i_row < N_rows_perBlock; i_row++) {
                myint ind_inBolck = j_col * N_rows_perBlock + i_row;

                // Upper triangular blocks Di i+1
                myint ind_inS_upperTri = (j_col + (i_layer+1) * N_rows_perBlock)* N_rows_S + i_row + i_layer * N_rows_perBlock;
                cascadedS.vals[ind_inS_upperTri] = portportBlocks[i_layer][1].vals[ind_inBolck];
                
                // Lower triangular blocks Di i-1
                myint ind_inS_lowerTri = (j_col + i_layer * N_rows_perBlock)* N_rows_S + i_row + (i_layer + 1) * N_rows_perBlock;
                cascadedS.vals[ind_inS_lowerTri] = portportBlocks[i_layer][2].vals[ind_inBolck];
            }
        }
    }

    // First diagonal block
    tempD = portportBlocks[0][0];
    for (myint j_col = 0; j_col < N_rows_perBlock; j_col++) {
        for (myint i_row = 0; i_row < N_rows_perBlock; i_row++) {
            myint ind_inBolck = j_col * N_rows_perBlock + i_row;
            myint ind_inS = j_col * N_rows_S + i_row;
            cascadedS.vals[ind_inS] = tempD.vals[ind_inBolck];
        }
    }

    // Last diagonal block
    tempD = portportBlocks[N_layers - 1][3];
    for (myint j_col = 0; j_col < N_rows_perBlock; j_col++) {
        for (myint i_row = 0; i_row < N_rows_perBlock; i_row++) {
            myint ind_inBolck = j_col * N_rows_perBlock + i_row;
            myint ind_inS = (j_col + N_layers * N_rows_perBlock)* N_rows_S + i_row + N_layers * N_rows_perBlock;
            cascadedS.vals[ind_inS] = tempD.vals[ind_inBolck];
        }
    }
    
    // The overlapped diagonal blocks (D22 and next D11')
    for (myint i_layer = 1; i_layer < N_layers; i_layer++) {
        tempD = portportBlocks[i_layer - 1][3].add(portportBlocks[i_layer][0]);
        for (myint j_col = 0; j_col < N_rows_perBlock; j_col++) {
            for (myint i_row = 0; i_row < N_rows_perBlock; i_row++) {
                myint ind_inBolck = j_col * N_rows_perBlock + i_row;
                myint ind_inS = (j_col + i_layer * N_rows_perBlock)* N_rows_S + i_row + i_layer * N_rows_perBlock;
                cascadedS.vals[ind_inS] = tempD.vals[ind_inBolck];
            }
        }
    }
    return cascadedS;
}

// Find the surface index of each port (layer growth along Y)
vector<myint> findSurfLocationOfPort(fdtdMesh *psys) {
/*  Method explained:
        - Each sourcePort has its position stored at psys->portCoor[sourcePort].y1, .y2, etc.
        - All node positions are stored at psys->yn, xn, zn.
        - Search psys->portCoor[sourcePort].y1 in psys->yn to find the surface index for that node.
    Notes: 
        - Many ports at 1 surf are counted once
        - For layer growth along y only and source current along x or z only

    Example: Nlayer = psys->N_cell_y = 3
        layer:      0s 0v 1s 1v 2s 2v 3s
        surfId:     0     1     2     3
        if port at: 0           2
        returned surfLocationOfPort = { 0, 2 }      */

    vector<myint> surfLocationOfPort;
    double dyErrorUpperBound = (psys->yn[1] - psys->yn[0]) * MINDISFRACY; // used to compare two equal y coordinates

    for (myint sourcePort = 0; sourcePort < psys->numPorts; sourcePort++) {
        double y_thisPort = psys->portCoor[sourcePort].y1[0];
        for (myint indy = 0; indy < psys->N_cell_y + 1; indy++) {   // search over all y nodes
            if (fabs(y_thisPort - psys->yn[indy]) < dyErrorUpperBound) {
                surfLocationOfPort.push_back(indy);
                break;
            }
        }
    }

    if (surfLocationOfPort.size() != psys->numPorts) {
        cout << "ERROR! Missed some ports." << endl;
        exit(2);
    }

    // Sort vector and erase duplicated surface index
    sort(surfLocationOfPort.begin(), surfLocationOfPort.end());
    surfLocationOfPort.erase( unique(surfLocationOfPort.begin(), surfLocationOfPort.end()), surfLocationOfPort.end() );

    if (surfLocationOfPort.size() != psys->numPorts) {
        cout << "ERROR! Unsupported port setting! Ports should not locate at the same surface" << endl;
        exit(2);
    }

    return surfLocationOfPort;
}

denseFormatOfMatrix cascadeMatrixS(fdtdMesh *psys, double omegaHz, const mapIndex &indexMap) {
    /* This function cascades original S matrix in COO format to a dense matrix with only port surfaces left.
    Inputs:
        - omegaHz (w): objective angular frequency in unit Hz
        - matrix ShSe/mu:
                (psys->SRowId, psys->SColId, psys->Sval) ~ COO format, index mode (growZ, removed PEC)
                psys->leng_S: number of nnz in ShSe/mu. Note that nnz at PEC has already been removed
        - indexMap: contains maps between different index modes
    Return:
        - cascadedS: matrix "(-w^2*D_eps+iw*D_sig+ShSe/mu)" with only {e}_portSurf left    */

    // Num of {e} at each surface or each layer, index mode (growY, removed PEC)
    myint n_surfExEz            = indexMap.N_surfExEz_rmPEC;
    myint n_volEy               = indexMap.N_volEy_rmPEC;
    myint n_layerE_growY        = n_surfExEz + n_volEy;
    myint N_layers              = psys->N_cell_y;       // num of layers

#ifdef DEBUG_SOLVE_REORDERED_S
    BlockType coo_reorderedS;
    coo_reorderedS.reserve(psys->leng_S);
    for (myint i_nnz = 0; i_nnz < psys->leng_S; i_nnz++) {
        // (rowId, colId, val) of this nnz element in ShSe/mu, index mode (growZ, removed PEC)
        myint nnzS_rowId_rmPECz = psys->SRowId[i_nnz];
        myint nnzS_colId_rmPECz = psys->SColId[i_nnz];
        complex<double> cascadedS_val = psys->Sval[i_nnz];

        // Map row and col index from (growZ, removed PEC) to (growY, removed PEC)
        myint nnzS_rowId = indexMap.eInd_map_rmPEC_z2y[nnzS_rowId_rmPECz];
        myint nnzS_colId = indexMap.eInd_map_rmPEC_z2y[nnzS_colId_rmPECz];

        // For diagonal nnz, add "-w^2*eps+iw*sig" to psys->Sval (ShSe/mu)
        if (nnzS_rowId == nnzS_colId) {     // if at diagonal
            myint nnzS_rowId_growZ = psys->mapEdgeR[nnzS_rowId_rmPECz];
            myint layerInd_alongZ = (nnzS_rowId_growZ + psys->N_edge_v) / (psys->N_edge_s + psys->N_edge_v);
            double epsi_thisnnz = psys->stackEpsn[layerInd_alongZ] * EPSILON0;
            double sigma_thisnnz = 0;
            if (psys->markEdge[nnzS_rowId_growZ] != 0) {
                sigma_thisnnz = SIGMA;
            }   // if inside a conductor 

            complex<double> epsi_sigma = { -omegaHz*omegaHz*epsi_thisnnz, omegaHz*sigma_thisnnz };
            cascadedS_val += epsi_sigma;    // -w^2*D_eps+iw*D_sig+ShSe/mu
        }
        coo_reorderedS.push_back({ nnzS_rowId, nnzS_colId, cascadedS_val });
    }

    ofstream file_obj;
    file_obj.open("sys_cooS_growYrmPEC.txt", ios::out);
    file_obj << "rowInd  colInd  val.real val.imag \n";
    sort(coo_reorderedS.begin(), coo_reorderedS.end(), ascendByRowIndThenByColInd);
    for (const auto &nnz : coo_reorderedS) {
        file_obj << nnz.row_ind << ' ' << nnz.col_ind << ' ' << nnz.val.real() << ' ' << nnz.val.imag() << endl;
    }
    file_obj.close();

    //csrFormatOfMatrix csrS_reordered(indexMap.N_totEdges_rmPEC, indexMap.N_totEdges_rmPEC, psys->leng_S);
    //sort(coo_reorderedS.begin(), coo_reorderedS.end(), ascendByRowIndThenByColInd);
    //csrS_reordered.convertBlockTypeToCsr(coo_reorderedS);
    //return csrS_reordered;
    denseFormatOfMatrix denseS_reordered(indexMap.N_totEdges_rmPEC, indexMap.N_totEdges_rmPEC);
    denseS_reordered.convertBlockTypeToDense(coo_reorderedS);
    return denseS_reordered;

#endif

    // Determine the block id of each nnz element of S and store in corresponding block matrix
    vector< BlockType > Blocks(8 * N_layers + 1);       // store all bolck matrices of S in a 2D vector
    for (myint i_nnz = 0; i_nnz < psys->leng_S; i_nnz++) {
        // (rowId, colId, val) of this nnz element in ShSe/mu, index mode (growZ, removed PEC)
        myint nnzS_rowId_rmPECz = psys->SRowId[i_nnz];
        myint nnzS_colId_rmPECz = psys->SColId[i_nnz];
        complex<double> cascadedS_val = psys->Sval[i_nnz];

        // Map row and col index from (growZ, removed PEC) to (growY, removed PEC)
        myint nnzS_rowId = indexMap.eInd_map_rmPEC_z2y[nnzS_rowId_rmPECz];
        myint nnzS_colId = indexMap.eInd_map_rmPEC_z2y[nnzS_colId_rmPECz];

        // For diagonal nnz, add "-w^2*eps+iw*sig" to psys->Sval (ShSe/mu)
        if (nnzS_rowId == nnzS_colId) {     // if at diagonal
            myint nnzS_rowId_growZ = psys->mapEdgeR[nnzS_rowId_rmPECz];
            myint layerInd_alongZ = (nnzS_rowId_growZ + psys->N_edge_v) / (psys->N_edge_s + psys->N_edge_v);
            double epsi_thisnnz = psys->stackEpsn[layerInd_alongZ] * EPSILON0;
            double sigma_thisnnz = 0;
            if (psys->markEdge[nnzS_rowId_growZ] != 0) {
                sigma_thisnnz = SIGMA;
            }   // if inside a conductor 

            complex<double> epsi_sigma = { -omegaHz*omegaHz*epsi_thisnnz, omegaHz*sigma_thisnnz };
            cascadedS_val += epsi_sigma;    // -w^2*D_eps+iw*D_sig+ShSe/mu
        }

        // Determine which block this nnz is at
        myint B_rowId = nnzS_rowId / n_layerE_growY * 2 + (nnzS_rowId % n_layerE_growY) / n_surfExEz;
        myint B_colId = nnzS_colId / n_layerE_growY * 2 + (nnzS_colId % n_layerE_growY) / n_surfExEz;
        myint BlockId = indexMap.mapBlockRowColToBlockInd(B_rowId, B_colId);

        // Shift the start row index and col index to be 0 inside each block
        myint nnzS_rowId_inBlock = nnzS_rowId - (B_rowId / 2) * n_layerE_growY - (B_rowId % 2) * n_surfExEz;
        myint nnzS_colId_inBlock = nnzS_colId - (B_colId / 2) * n_layerE_growY - (B_colId % 2) * n_surfExEz;

        // Store this nnz at corresponding block matrix
        Blocks[BlockId].push_back({ nnzS_rowId_inBlock, nnzS_colId_inBlock, cascadedS_val });
    }

    //// Free original matrix S to save memory
    //free(psys->SRowId);
    //free(psys->SColId);
    //free(psys->Sval);

    /******************** Start Cascading Matrix S from Blocks ***************************/
    // Tell the size of sparse block matrices after removing edges at PEC
    myint N_surfE = indexMap.N_surfExEz_rmPEC;
    myint N_volE = indexMap.N_volEy_rmPEC;

    // Sort nnzs in each block matrix by row index first then by col index
    /* The purpose is to allow easy convert from COO to CSR format*/
    for (auto &block : Blocks) {
        sort(block.begin(), block.end(), ascendByRowIndThenByColInd);
    }

    /*ofstream file_obj;
    for (int i_block = 0; i_block < Blocks.size(); i_block++) {
        string filename = "block_B" + to_string(i_block) + ".txt";
        file_obj.open(filename, ios::out);
        file_obj << "rowInd  colInd  val.real val.imag \n";
        for (const auto &nnz : Blocks[i_block]) {
            file_obj << nnz.row_ind << ' ' << nnz.col_ind << ' ' << nnz.val.real() << ' ' << nnz.val.imag() << endl;
        }
        file_obj.close();
    }*/
    

    // Half the value of overlapped ns-ns blocks between adjcent two layers
    /* The partitioned blocks no longer mean physical curl-curl opeartor, but mamatically, this is doable to
    cascade block matrix S. Any partition works like C = aC' + bC''. C = C' + C' makes codeing easier.*/
    for (myint nsns_BlockId = 8; nsns_BlockId < 8 * N_layers; nsns_BlockId += 8) {
        for (auto &nnz : Blocks[nsns_BlockId]) {
            nnz.val *= 0.5;
        }
    }

    // Eliminate e_volume at each layer and store 4 reduced dense blocks of each layer
    vector<vector<denseFormatOfMatrix>> surfSurfBlocks(N_layers);   // N_layers*4 dense blocks
    for (myint i_layer = 0; i_layer < N_layers; i_layer++) {
        // Init and allocate memory for 4 dense blocks of each layer
        surfSurfBlocks[i_layer].reserve(4);
        for (myint i_block = 0; i_block < 4; i_block++) {           
            surfSurfBlocks[i_layer].push_back(denseFormatOfMatrix(N_surfE, N_surfE));
        }

        // From 9 blocks at this layer to 4, surfSurfBlocks[i_layer] = {C11, C12, C21, C22}
        vector<BlockType> layerS(Blocks.begin() + 8 * i_layer, Blocks.begin() + 8 * i_layer + 9);
        eliminateVolumE(layerS, N_surfE, N_volE, surfSurfBlocks[i_layer].data());
        layerS.clear();
    }
    Blocks.clear();     // free blocks of original whole matrix S to save memory

    // Cascade surf-surf blocks and only keep layers where ports are
    vector<myint> surfLocationOfPort = findSurfLocationOfPort(psys);
    
    /* Each layer: surfSurfBlocks[i_layer] = {C11, C12, C21, C22} */

    // Left -> first port: cascaded C22' = C22 - C21*inv(C11)*C12
    for (myint i_layer = 0; i_layer < surfLocationOfPort.front(); i_layer++) {
        surfSurfBlocks[i_layer][3] = surfSurfBlocks[i_layer][3].minus(
            surfSurfBlocks[i_layer][2].dot(surfSurfBlocks[i_layer][0].backslash(surfSurfBlocks[i_layer][1])));

        if (i_layer != N_layers) {  // if not reaching the right most layer
            // Add cascaded C22' at this layer to C11 at next layer
            surfSurfBlocks[i_layer + 1][0] = surfSurfBlocks[i_layer + 1][0].add(surfSurfBlocks[i_layer][3]);
        }
    }

    // Right -> last port: cascaded C11' = C11 - C12*inv(C22)*C21
    for (myint i_layer = N_layers-1; i_layer >= surfLocationOfPort.back(); i_layer--) {
        surfSurfBlocks[i_layer][0] = surfSurfBlocks[i_layer][0].minus(
            surfSurfBlocks[i_layer][1].dot(surfSurfBlocks[i_layer][3].backslash(surfSurfBlocks[i_layer][2])));

        if (i_layer != 0) {         // if not reaching the left most layer
            // Add cascaded C11' at this layer to C22 at previous layer
            surfSurfBlocks[i_layer - 1][3] = surfSurfBlocks[i_layer - 1][3].add(surfSurfBlocks[i_layer][0]);
        }
    }

    // Middle: 
    /* For the middle, cascaded the middle surface ({e}_midLayer) between 2 layers
         For example:
            |C11     C12          | {e}_thisPortLayer
            |C21   C22+C11'   C12'| {e}_midLayer       , where  4 Cij at thisPortLayer
            |        C21'     C22'| {e}_midLayer+1              4 Cij' at the midLayer
            To cascade the middle surface {e}_midLayer, the precedures are:
                1) tempC22 = C22+C11'
                2) tempD21 = tempC22\C21, tempD12 = tempC22\C12'
                3) cascaded 4 blocks still store in Cij at thisPortLayer:
                    cascaded C11 = C11 - C12*tempD21,  cascaded C12 = 0 - C12*tempD12
                    cascaded C21 =  0 - C21'*tempD21,  cascaded C22 = C22' - C21'*tempD12
      
      Then, the cascaded 4 blocks for each port (index 0 <= i_port < N_ports-1, except last port) will be stored in corresponding layer.
            let: thisPortLayer = surfLocationOfPort[i_port];
                 nextPortLayer = surfLocationOfPort[i_port + 1];
            i_port's 4 blocks are stored in:
                 surfSurfBlocks[thisPortLayer] = {D11, D12, D21, D22}
            next port's 4 blocks are stored in:
                 surfSurfBlocks[nextPortLayer] = {D11', D12', D21', D22'}
            To reconstruct cascaded S matrix, we need to organize above cascaded blocks as
                |D11     D12          | {e}_thisPort
                |D21   D22+D11'   D12'| {e}_nextPort
                |        D21'     D22'| {e}_nextnextPort
    */
    denseFormatOfMatrix tempC22(N_surfE, N_surfE);              // tempC22 = C22+C11'
    denseFormatOfMatrix tempD21(N_surfE, N_surfE);              // tempD21 = tempC22\C21
    denseFormatOfMatrix tempD12(N_surfE, N_surfE);              // tempD12 = tempC22\C12'
    for (myint i_port = 0; i_port < surfLocationOfPort.size() - 1; i_port++) {                  // all the intervals between ports
        myint thisPortLayer = surfLocationOfPort[i_port];       // surface index and layer index of this port (i_port)
        myint nextPortLayer = surfLocationOfPort[i_port + 1];   // surface index and layer index of next port
        for (myint i_midLayer = thisPortLayer + 1; i_midLayer < nextPortLayer; i_midLayer++) {    // all middle layers between 2 ports
            tempC22 = surfSurfBlocks[thisPortLayer][3].add(surfSurfBlocks[i_midLayer][0]);
            tempD21 = tempC22.backslash(surfSurfBlocks[thisPortLayer][2]);
            tempD12 = tempC22.backslash(surfSurfBlocks[i_midLayer][1]);
            surfSurfBlocks[thisPortLayer][0] = surfSurfBlocks[thisPortLayer][0].minus(
                surfSurfBlocks[thisPortLayer][1].dot(tempD21));
            surfSurfBlocks[thisPortLayer][1] = surfSurfBlocks[thisPortLayer][1].dot(tempD12).multiplyScalar(-1.0);
            surfSurfBlocks[thisPortLayer][2] = surfSurfBlocks[i_midLayer][2].dot(tempD21).multiplyScalar(-1.0);
            surfSurfBlocks[thisPortLayer][3] = surfSurfBlocks[i_midLayer][3].minus(
                surfSurfBlocks[i_midLayer][2].dot(tempD12));
        }
    }
    
    // Collect all cascaded blocks only related to port surfaces
    myint N_ports = surfLocationOfPort.size();
    if (N_ports == 1) {     // when only one port, return the block C11
        myint thisPortLayer = surfLocationOfPort[0];
        return surfSurfBlocks[thisPortLayer][0];
    }
    vector<vector<denseFormatOfMatrix>> portportBlocks(N_ports-1);   // (N_ports-1)*4 dense blocks
    for (myint i_port = 0; i_port < N_ports - 1; i_port++) {
        myint thisPortLayer = surfLocationOfPort[i_port];
        portportBlocks[i_port] = surfSurfBlocks[thisPortLayer];
    }
    surfSurfBlocks.clear(); // free surf-surf blocks 

    return reconstructBlocksToDense(portportBlocks);
}

// Calculate all the simulated freq points and store in a vector
vector<double> calAllFreqPointsHz(const fdtdMesh &sys) {
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

// Assign source current density (-iw{j}) for every port excitation, index mode (growY, removed PEC, cascaded)
denseFormatOfMatrix assignRhsJForAllPorts(fdtdMesh *psys, double omegaHz, const mapIndex &indexMap) {
    /* -iw{j} for each port is stored continusly in denseFormatOfMatrix (col by col)
    
    - step 1: get psys->portCoor[sourcePort].portEdge[sourcePortSide][], which stored the global edge index under
              mode (growZ, no PEC removal) of all excitation current lines when excited at port index "sourcePort"
    - step 2: map to (growY, removed PEC)
    - step 3: cascade to port surfaces    */

    // Map -iw{j} at at edges from (growZ, no PEC removal) to (growY, removed PEC)
    myint N_edgesNotAtPEC = indexMap.N_totEdges - indexMap.N_edgesAtPEC;                        // num of edges not at PEC
    denseFormatOfMatrix RhsJ_SI(N_edgesNotAtPEC, psys->numPorts);                               // store -iw{j} of (growY, removed PEC)
    for (myint sourcePort = 0; sourcePort < psys->numPorts; sourcePort++) {
        for (myint sourcePortSide = 0; sourcePortSide < psys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
            for (myint inde = 0; inde < psys->portCoor[sourcePort].portEdge[sourcePortSide].size(); inde++) {
                myint ind_j_growZ = psys->portCoor[sourcePort].portEdge[sourcePortSide][inde];  // j edge index (growZ, no PEC removal)
                myint ind_j_growY = indexMap.eInd_map_z2y[ind_j_growZ];                         // -> (growY, no PEC removal)
                myint ind_j_growYrmPEC = indexMap.eInd_map_y2rmPEC[ind_j_growY];                // -> (growY, removed PEC)
                ind_j_growYrmPEC += sourcePort * N_edgesNotAtPEC;                               // index shifted by port number
                RhsJ_SI.vals[ind_j_growYrmPEC] = { 0.0, -1.0 * omegaHz * (psys->portCoor[sourcePort].portDirection[sourcePortSide]) };
            }
        }
    }

#ifdef DEBUG_SOLVE_REORDERED_S
    return RhsJ_SI;
#endif

    // Only keep -iw{j} at port surfaces (cascading)
    myint N_surfExEz = indexMap.N_surfExEz_rmPEC;                               // number of edges at one surface
    myint n_volEy = indexMap.N_volEy_rmPEC;                                     // number of volume edges in one layer
    myint N_allCascadedEdges = psys->numPorts * N_surfExEz;                     // number of all edges in kept surfaces
    denseFormatOfMatrix cascadedRhsJ_SI(N_allCascadedEdges, psys->numPorts);    // -iw{j} only at port surfaces
    vector<myint> surfLocationOfPort = findSurfLocationOfPort(psys);
    for (myint excitedPort = 0; excitedPort < psys->numPorts; excitedPort++) {
        for (myint ind_thisPort = 0; ind_thisPort < psys->numPorts; ind_thisPort++) {           // for each kept port surface
            for (myint j_ind = 0; j_ind < N_surfExEz; j_ind++) {    // for all e at this surface, index starting from 0
                // Cascaded index = indJ at this surface + shift by previous port surfaces + shift by previous port excitation
                myint indJ_cascaded = j_ind + ind_thisPort * N_surfExEz + excitedPort * N_allCascadedEdges;

                // For original noncascaded index
                myint surfInd_thisPort = surfLocationOfPort[ind_thisPort];
                myint shift_surfvol = surfInd_thisPort * (N_surfExEz + n_volEy);    // index shift from all surf-vol {e} before this surface
                myint indJ_noncasc = j_ind + shift_surfvol + excitedPort * N_edgesNotAtPEC;

                cascadedRhsJ_SI.vals[indJ_cascaded] = RhsJ_SI.vals[indJ_noncasc];
            }   
        }   
    }

    return cascadedRhsJ_SI;
}

// Reconstruct {e} (growZ, removed PEC) from {e} at index mode (growY, removed PEC, cascaded)
denseFormatOfMatrix reconstruct_e(fdtdMesh *psys, const mapIndex &indexMap, const denseFormatOfMatrix &cascadedeField_SI, double excitedPort) {
    /*  Inputs:
            - indexMap: contains maps between different index modes
            - cascadedeField_SI: of size N_allCascadedEdges * psys->numPorts, mode (growY, removed PEC, cascaded)
            - excitedPort: excited port index
        Return:
            - eField_oneExcit: of size (psys->N_edge - psys->bden) * 1, mode (growZ, removed PEC)    */

    myint N_surfExEz = indexMap.N_surfExEz_rmPEC;                               // number of edges at one surface
    myint n_volEy = indexMap.N_volEy_rmPEC;                                     // number of volume edges in one layer
    myint N_allCascadedEdges = psys->numPorts * N_surfExEz;                     // number of all edges in kept surfaces
    myint N_edgesNotAtPEC = psys->N_edge - psys->bden;                          // num of all edges after removing PEC
    vector<myint> surfLocationOfPort = findSurfLocationOfPort(psys);            // surf index of the port's location

    // Reconstruct {e} solution when excited at port index "excitedPort"
    denseFormatOfMatrix eField_oneExcit(N_edgesNotAtPEC, 1);
    for (myint ind_thisPort = 0; ind_thisPort < psys->numPorts; ind_thisPort++) {           // for each kept port surface
        for (myint e_ind = 0; e_ind < N_surfExEz; e_ind++) {    // for all e at this surface, index starting from 0
            // Cascaded index = e_ind at this surface + shift by previous port surfaces + shift by previous port excitation
            myint eInd_cascaded = e_ind + ind_thisPort * N_surfExEz + excitedPort * N_allCascadedEdges;

            // For original noncascaded index, mode (growZ, removed PEC)
            myint surfInd_thisPort = surfLocationOfPort[ind_thisPort];
            myint shift_surfvol = surfInd_thisPort * (N_surfExEz + n_volEy);    // index shift from all surf-vol {e} before this surface
            myint eInd_growYrmPEC = e_ind + shift_surfvol;                      // Step 0: -> index at (growY, removed PEC)
            myint eInd_growY = indexMap.eInd_map_rmPEC2y[eInd_growYrmPEC];      // Step 1: -> index at (growY, no PEC removal)
            myint eInd_growZ = indexMap.eInd_map_y2z[eInd_growY];               // Step 2: map (growY, no PEC removal) to (growZ, no PEC removal)
            myint eInd_growZrmPEC = psys->mapEdge[eInd_growZ];                  // Step 3: map (growZ, no PEC removal) to (growZ, removed PEC)
            
            eField_oneExcit.vals[eInd_growZrmPEC] = cascadedeField_SI.vals[eInd_cascaded];
        }
    }

#ifdef DEBUG_SOLVE_REORDERED_S
    for (myint e_ind = 0; e_ind < N_allCascadedEdges; e_ind++) {
        myint eInd_cascaded = e_ind + excitedPort * N_allCascadedEdges;

        myint eInd_growY = indexMap.eInd_map_rmPEC2y[e_ind];                // Step 1: -> index at (growY, no PEC removal)
        myint eInd_growZ = indexMap.eInd_map_y2z[eInd_growY];               // Step 2: map (growY, no PEC removal) to (growZ, no PEC removal)
        myint eInd_growZrmPEC = psys->mapEdge[eInd_growZ];                  // Step 3: map (growZ, no PEC removal) to (growZ, removed PEC)

        eField_oneExcit.vals[eInd_growZrmPEC] = cascadedeField_SI.vals[eInd_cascaded];
    }
    return eField_oneExcit;
#endif

    return eField_oneExcit;
}

// Solve E field and Z-parameters in Pardiso, solve layer by layer. (under developing)
void solveE_Zpara_layered(fdtdMesh *psys) {

    // Initilize Z-parameters for all frequencies
    psys->x.assign(psys->numPorts * psys->numPorts * psys->nfreq, complex<double>(0., 0.));

    // Set up necessary index map
    mapIndex indexMap(psys->N_cell_x, psys->N_cell_y, psys->N_cell_z);
    indexMap.setEdgeMap_growZgrowY();
    indexMap.setEdgeMap_growYremovePEC(psys->ubde, psys->lbde, psys->bden);
    indexMap.setEdgeMap_rmPEC_growZgrowY(psys->mapEdgeR);

    vector<double> vFreqHz = calAllFreqPointsHz(*psys);
    for (int indFreq = 0; indFreq < vFreqHz.size(); indFreq++) {    // for each computed freq point
        double omegaHz = 2.0 * M_PI * vFreqHz[indFreq];

        // Cascaded system matrix (-w^2*D_eps+iw*D_sig+ShSe/mu)
        denseFormatOfMatrix cascadedS = cascadeMatrixS(psys, omegaHz, indexMap);

        // Cascaded -iwJ and {e}. All excitations at each port are solved together
        denseFormatOfMatrix cascadedRhsJ_SI = assignRhsJForAllPorts(psys, omegaHz, indexMap);   // -iw{j} in unit (A * m^-2 / s)
        denseFormatOfMatrix cascadedeField_SI = cascadedS.backslash(cascadedRhsJ_SI);           // {e} in unit (V/m)

        // For each port excitation, reconstruct {e} and solve Z-parameters
        for (myint excitedPort = 0; excitedPort < psys->numPorts; excitedPort++) {
            denseFormatOfMatrix eField_oneExcit = reconstruct_e(psys, indexMap, cascadedeField_SI, excitedPort);
            psys->Construct_Z_V0_Vh(eField_oneExcit.vals.data(), indFreq, excitedPort);
        }
    }

    // Print Z-parameters
    psys->print_z_V0_Vh();

}

// Solve E field and Z-parameters in Pardiso, solve the whole structure as reference
void solveE_Zpara_reference(fdtdMesh *psys) {

    // All computed freq points
    vector<double> vFreqHz = calAllFreqPointsHz(*psys);

    // Initilize Z-parameters for all frequencies
    psys->x.assign(psys->numPorts * psys->numPorts * psys->nfreq, complex<double>(0., 0.));

    // At each frequency, solve all ports together
    for (int indFreq = 0; indFreq < vFreqHz.size(); indFreq++) {
        reference(psys, indFreq, psys->SRowId, psys->SColId, psys->Sval);
    }   

    // Print Z-parameters
    psys->print_z_V0_Vh();

}

#endif