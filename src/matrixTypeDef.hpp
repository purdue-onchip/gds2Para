#ifndef GDS2PARA_MATRIX_TYPEDEF_H_
#define GDS2PARA_MATRIX_TYPEDEF_H_

#include "fdtd.hpp"

// COO format
typedef struct {
    myint row_ind;
    myint col_ind;
    complex<double> val;
}onennzType;                            // store the row-col-val of one nnz element
typedef vector<onennzType> BlockType;   // store rows-cols-vals of all nnzs inside each block matrix
bool ascendByRowIndThenByColInd(const onennzType& nnz1, const onennzType& nnz2) {
    return tie(nnz1.row_ind, nnz1.col_ind) < tie(nnz2.row_ind, nnz2.col_ind);
}                                       // operator for sorting nnz by row_ind first then by col_ind

// 0-based column-major dense format of a matrix stored col by col in 1-D vector
class denseFormatOfMatrix {
public:
    // matrix information
    myint N_rows;
    myint N_cols;
    myint matrixSize;
    vector<complex<double>> vals;

    // Constructor
    denseFormatOfMatrix(myint N_rows, myint N_cols) {
        this->N_rows = N_rows;
        this->N_cols = N_cols;
        this->matrixSize = N_rows * N_cols;

        this->vals.reserve(this->matrixSize);
        this->vals.assign(this->matrixSize, { 0.0, 0.0 });
    }

    void writeToFile(string filename) {
        ofstream file_obj;
        file_obj.open(filename, ios::out);
        file_obj << "N_rows  N_cols \n";
        file_obj << this->N_rows << ' ' << this->N_cols << endl;
        file_obj << "val.real  val.imag \n";
        for (const auto val : this->vals) {
            file_obj << val.real() << ' ' << val.imag() << endl;
        }
        file_obj.close();
    }

    // Convert BlockType (COO) to dense format.
    void convertBlockTypeToDense(const BlockType &block) {
        for (const auto &onennz : block) {
            myint row_ind = onennz.row_ind;
            myint col_ind = onennz.col_ind;
            myint ind_inArray = col_ind * this->N_rows + row_ind;
            this->vals[ind_inArray] = onennz.val;
        }
    }

    // Copy this->vals to a MKL_Complex16 type memory
    void copyToMKL_Complex16(MKL_Complex16 *ptr) const {
        /* Input: a pointer to assigned MKL_Complex16 memory of size this->matrixSize */
        for (myint ind = 0; ind < this->matrixSize; ind++) {
            ptr[ind].real = this->vals[ind].real();
            ptr[ind].imag = this->vals[ind].imag();
        }
    }

    // Copy this->vals from a MKL_Complex16 type memory
    void copyFromMKL_Complex16(MKL_Complex16 *ptr) {
        /* Input: a pointer to assigned MKL_Complex16 memory of size this->matrixSize */
        for (myint ind = 0; ind < this->matrixSize; ind++) {
            this->vals[ind] = { ptr[ind].real , ptr[ind].imag };
        }
    }

    // C = A (this) + B
    denseFormatOfMatrix add(const denseFormatOfMatrix &B) const {
        if (this->N_rows != B.N_rows && this->N_cols != B.N_cols) {
            cout << "Failure to add matrices, dimensions not match!" << endl;
            exit(2);
        }

        denseFormatOfMatrix C(this->N_rows, this->N_cols);
        for (myint ind = 0; ind < this->matrixSize; ind++) {
            C.vals[ind] = this->vals[ind] + B.vals[ind];
        }
        return C;
    }

    // C = A (this) - B
    denseFormatOfMatrix minus(const denseFormatOfMatrix &B) {
        if (this->N_rows != B.N_rows && this->N_cols != B.N_cols) {
            cout << "Failure to minus matrice, dimensions not match!" << endl;
            exit(2);
        }

        denseFormatOfMatrix C(this->N_rows, this->N_cols);
        for (myint ind = 0; ind < this->matrixSize; ind++) {
            C.vals[ind] = this->vals[ind] - B.vals[ind];
        }
        return C;
    }

    // C = A (this) * scalar
    denseFormatOfMatrix multiplyScalar(double scalar) {
        denseFormatOfMatrix C(this->N_rows, this->N_cols);
        for (myint ind = 0; ind < this->matrixSize; ind++) {
            C.vals[ind] = this->vals[ind] * scalar;
        }
        return C;
    }

    // C = A (this) dot B (matrix-matrix dot multiply, not element-wise multiply)
    denseFormatOfMatrix dot(const denseFormatOfMatrix &B) {
        if (this->N_cols != B.N_rows) {
            cout << "Failure to dot multiply matrices, dimensions not match!" << endl;
            exit(2);
        }

        // Use mkl cblas_?gemm (C := alpha*op(A)*op(B) + beta*C) to get C=A*B
        denseFormatOfMatrix C(this->N_rows, B.N_cols);
        complex<double> alpha(1.0, 0.0);
        complex<double> beta(0.0, 0.0);
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            this->N_rows, B.N_cols, this->N_cols,       // matrix dimensions ~ m, n, k
            &alpha,                                     // alpha = 1
            this->vals.data(), this->N_rows,            // A
            B.vals.data(), this->N_cols,                // B
            &beta,                                      // beta = 0
            C.vals.data(), this->N_rows);               // C
        return C;
    }

    // C = A (this) \ B = inv(A) dot B
    denseFormatOfMatrix backslash(const denseFormatOfMatrix &B) {
        if (this->N_rows != this->N_cols) {
            cout << "Failure, matrix is not invertable!" << endl;
            exit(2);
        }
        if (this->N_cols != B.N_rows) {
            cout << "Failure to backslash matrices, dimensions not match!" << endl;
            exit(2);
        }

        // Compute the LU factorization of a general m-by-n matrix A (this) and store in A_LU
        vector<MKL_Complex16> A_LU(this->matrixSize);
        this->copyToMKL_Complex16(A_LU.data());         // copy matrix A to a MKL_Complex16 type A_LU
        vector<myint> ipiv(this->N_rows);
        lapack_int info = LAPACKE_zgetrf(LAPACK_COL_MAJOR,
            this->N_rows, this->N_cols,                 // matrix dimensions ~ m, n
            A_LU.data(),                                // matrix A
            this->N_rows,
            ipiv.data()                                 // pivot indices
        );
        if (info != 0) {
            cout << "Issue on LU factorization, LAPACKE_?getrf returns: " << info << endl;
            exit(2);
        }

        // Solve C = A_LU \ B with LAPACKE_?getrs
        vector<MKL_Complex16> C_mklComplex(B.matrixSize);
        B.copyToMKL_Complex16(C_mklComplex.data());     // copy matrix B to a MKL_Complex16 type C_mklComplex
        info = LAPACKE_zgetrs(LAPACK_COL_MAJOR, 'N',
            B.N_rows, B.N_cols,                         // n & nrhs
            A_LU.data(), B.N_rows,                      // A_LU
            ipiv.data(),                                // pivot indices of A_LU
            C_mklComplex.data(), B.N_rows               // C (copied B) will be overwritten by the solution matrix
        );
        if (info != 0) {
            cout << "Issue on mkl backslash, LAPACKE_?getrs returns: " << info << endl;
            exit(2);
        }

        //// Refines the solution C = A_LU \ B and estimate its error (Error < 1e-12, negligible at this step)
        //vector<MKL_Complex16> A_mklComplex(this->matrixSize);
        //this->copyToMKL_Complex16(A_mklComplex.data()); // copy matrix A to a MKL_Complex16 type A_mklComplex
        //vector<MKL_Complex16> B_mklComplex(B.matrixSize);
        //B.copyToMKL_Complex16(B_mklComplex.data());     // copy matrix B to a MKL_Complex16 type B_mklComplex
        //vector<double> ferr(B.N_cols, 0.0), berr(B.N_cols, 0.0);
        //info = LAPACKE_zgerfs(LAPACK_COL_MAJOR, 'N',
        //    B.N_rows, B.N_cols,                         // n & nrhs
        //    A_mklComplex.data(), B.N_rows,              // original A in MKL_Complex16 type
        //    A_LU.data(), B.N_rows,                      // A_LU
        //    ipiv.data(),                                // pivot indices of A_LU
        //    B_mklComplex.data(), B.N_rows,              // B
        //    C_mklComplex.data(), B.N_rows,              // C
        //    ferr.data(), berr.data()                    // forward and backward errors for each solution vector
        //);
        //if (info != 0) {
        //    cout << "Issue on backlash refinement, LAPACKE_?gerfs returns: " << info << endl;
        //    exit(2);
        //}

        // Store the solution at denseFormatOfMatrix
        denseFormatOfMatrix C(B.N_rows, B.N_cols);
        C.copyFromMKL_Complex16(C_mklComplex.data());   // copy C.vals (complex<double>) from MKL_Complex16

        return C;
    }
};

// 0-based row-major 3-array CSR format of a matrix
class csrFormatOfMatrix {
public:
    // matrix information
    myint N_rows;
    myint N_cols;
    myint N_nnz;
    vector<myint> rows;                  // CSR row indices, size = N_rows + 1
    vector<myint> cols;                  // CSR col indices, size = N_nnz
    vector<complex<double>> vals;        // CSR nnz values,  size = N_nnz

    // Constructor
    csrFormatOfMatrix(myint N_rows, myint N_cols, myint N_nnz) {
        /* Inputs:
        N_rows*N_cols:  matrix size;
        N_nnz:          num of nonzeros.     */

        this->N_rows = N_rows;
        this->N_cols = N_cols;
        this->N_nnz = N_nnz;

        this->rows.assign(N_rows + 1, 0);
        this->cols.assign(N_nnz, 0);
        this->vals.assign(N_nnz, { 0.0, 0.0 });
    }

    // Copy this->vals to a MKL_Complex16 type memory
    void copyToMKL_Complex16(MKL_Complex16 *ptr) const {
        /* Input: a pointer to assigned MKL_Complex16 memory of size this->N_nnz */
        for (myint ind = 0; ind < this->N_nnz; ind++) {
            ptr[ind].real = this->vals[ind].real();
            ptr[ind].imag = this->vals[ind].imag();
        }
    }

    // Convert BlockType (COO) to CSR
    void convertBlockTypeToCsr(const BlockType &block) {
        /* This function convert BlockType (COO) to CSR format.
        The input BlockType must have been sorted by row index first then col index
        within each row. All in increasing order. */

        // The first index of rows is always 0 and the last is always N_nnz
        this->rows[0] = 0;
        this->rows[this->N_rows] = block.size();

        myint nnz_ind = 0;
        myint thisRow_ind = 0;
        for (const auto &onennz : block) {

            // cols and vals are direct copy of the tuple for every nnz
            this->cols[nnz_ind] = onennz.col_ind;
            this->vals[nnz_ind] = onennz.val;

            // When saw nnz at a new row
            while (thisRow_ind < onennz.row_ind) {
                thisRow_ind++;
                this->rows[thisRow_ind] = nnz_ind;
            }

            nnz_ind++;

            // Extrame case: when the last few rows are all blank
            while (nnz_ind == block.size() && thisRow_ind < this->N_rows - 1) {
                thisRow_ind++;
                this->rows[thisRow_ind] = nnz_ind;
            }
        }
    }

    // Return denseD0sD1s = csrB22 (this) \ denseB21B23 = inv(csrB22)*denseB21B23
    denseFormatOfMatrix backslashDense(const denseFormatOfMatrix &denseB21B23) {
        // combined dense [D0s, D1s]
        denseFormatOfMatrix denseD0sD1s(denseB21B23.N_rows, denseB21B23.N_cols);

        // Pardiso parameters, see https://software.intel.com/en-us/mkl-developer-reference-c-pardiso
        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT mtype = 13;                 /* Complex and nonsymmetric matrix */
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
        //iparm[26] = 1;       /* Check whether column indices are sorted in increasing order within each row */
        iparm[3] = 0;          /* No iterative-direct algorithm */
        //iparm[59] = 2;       /* out of core version to solve very large problem */
        //iparm[10] = 0;       /* Use nonsymmetric permutation and scaling MPS */

        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &(this->N_rows), this->vals.data(), this->rows.data(), this->cols.data(),
            &perm, &nrhs, iparm, &msglvl, (complex<double>*)denseB21B23.vals.data(), denseD0sD1s.vals.data(), &error);
        if (error != 0) {
            printf("\nERROR during PARDISO backslash: %d", error);
            exit(2);
        }

        // Release internal memory
        phase = -1;
        complex<double> ddum;               /* Complex<double> dummy */
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &(this->N_rows), &ddum, this->rows.data(), this->cols.data(),
            &perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

        return denseD0sD1s;
    }
};

#endif