#ifndef GDS2PARA_MATRIX_TYPEDEF_H_
#define GDS2PARA_MATRIX_TYPEDEF_H_

#include "fdtd.hpp"

typedef vector< tuple<myint, myint, double> > BlockType;    // store rows-cols-vals of each block matrix

// 0-based row-major 3-array CSR format of a matrix
class csrFormatOfMatrix {
public:
    // matrix information
    myint N_rows;
    myint N_cols;
    myint N_nnz;
    myint *rows = nullptr;      // CSR row indices, size = N_rows + 1
    myint *cols = nullptr;      // CSR col indices, size = N_nnz
    double *vals = nullptr;      // CSR nnz values,  size = N_nnz

                                 // Constructor
    csrFormatOfMatrix(myint N_rows, myint N_cols, myint N_nnz) {
        /* Inputs:
        N_rows*N_cols:  matrix size;
        N_nnz:          num of nonzeros.     */

        this->N_rows = N_rows;
        this->N_cols = N_cols;
        this->N_nnz = N_nnz;

        this->rows = new myint[N_rows + 1]();
        this->cols = new myint[N_nnz]();
        this->vals = new double[N_nnz]();
    }

    // Destructor
    ~csrFormatOfMatrix() {
        delete[] this->rows;
        delete[] this->cols;
        delete[] this->vals;
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
        while (nnz_ind == block.size() && thisRow_ind < this->N_rows - 1) {
            thisRow_ind++;
            this->rows[thisRow_ind] = nnz_ind;
        }
    }

    return 0;
}

// 0-based column-major dense format of a matrix stored col by col in 1-D vector
class denseFormatOfMatrix {
public:
    // matrix information
    myint N_rows;
    myint N_cols;
    myint matrixSize;
    vector<double> vals;

    // Constructor
    denseFormatOfMatrix(myint N_rows, myint N_cols) {
        this->N_rows = N_rows;
        this->N_cols = N_cols;
        this->matrixSize = N_rows * N_cols;

        this->vals.reserve(this->matrixSize);
        this->vals.assign(this->matrixSize, 0.0);
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
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            this->N_rows, B.N_cols, this->N_cols,       // matrix dimensions ~ m, n, k
            1.0,                                        // alpha
            this->vals.data(), this->N_rows,            // A
            B.vals.data(), this->N_cols,                // B
            0.0,                                        // beta
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
        vector<double> A_LU(this->vals);
        vector<myint> ipiv(this->N_rows);
        lapack_int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR,
            this->N_rows, this->N_cols,     // matrix dimensions ~ m, n
            A_LU.data(),                    // matrix A
            this->N_rows,
            ipiv.data()                     // pivot indices
        );
        if (info != 0) {
            cout << "Issue on LU factorization, LAPACKE_dgetrf returns: " << info << endl;
            exit(2);
        }

        // Solve C = A_LU \ B with LAPACKE_dgetrs
        denseFormatOfMatrix C(B.N_rows, B.N_cols);
        C.vals = B.vals;                // copy matrix B to C
        info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N',
            B.N_rows, B.N_cols,         // n & nrhs
            A_LU.data(), B.N_rows,      // A_LU
            ipiv.data(),                // pivot indices of A_LU
            C.vals.data(), B.N_rows     // C (copied B) will be overwritten by the solution matrix
        );
        if (info != 0) {
            cout << "Issue on mkl backslash, LAPACKE_dgetrs returns: " << info << endl;
            exit(2);
        }

        return C;
    }
};

#endif