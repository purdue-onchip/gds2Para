#include <math.h>
#include <ctime>
#include "fdtd.h"
#include "hypreSolverh.h"

#include "vis.c"

int hypre_FlexGMRESModifyPCAMGExample(void *precond_data, HYPRE_Int iterations,
    double rel_residual_norm);


int hypreSolve(fdtdMesh *sys, HYPRE_IJMatrix A, HYPRE_ParCSRMatrix parcsr_A, myint leng_A, double *bin, myint leng_v0, double *solution){
    HYPRE_Int i;
    int myid, num_procs;
    HYPRE_Int N, n;

    HYPRE_Int ilower, iupper;
    HYPRE_Int local_size, extra;

    HYPRE_Int solver_id;
    HYPRE_Int vis, print_system;

    double h, h2;

    /*HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;*/
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;

    HYPRE_Solver solver, precond;
    
    /* Initialize MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    N = leng_v0;
    local_size = N / num_procs;
    extra = N - local_size * num_procs;

    ilower = local_size * myid;
    ilower += hypre_min(myid, extra);
    
    iupper = local_size * (myid + 1);
    iupper += hypre_min(myid + 1, extra);
    iupper = iupper - 1;
    
    /* How many rows do I have? */
    local_size = iupper - ilower + 1;
    
    ///* Create the matrix.
    //Note that this is a square matrix, so we indicate the row partition
    //size twice (since number of rows = number of cols) */
    //HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

    ///* Choose a parallel csr format storage (see the User's Manual) */
    //HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

    ///* Initialize before setting coefficients */
    //HYPRE_IJMatrixInitialize(A);
    //
    //{
    //    int nnz;
    //    vector<double> values;
    //    vector<int> cols;
    //    int index = 0;

    //    for (i = ilower; i <= iupper; i++)
    //    {
    //        nnz = 0;   // Number of non-zeros on row i

    //        while (ARowId[index] == i){
    //            cols.push_back(AColId[index]);
    //            values.push_back(Aval[index]);
    //            nnz++;
    //            index++;
    //        }
    //        
    //        /* Set the values for row i */
    //        HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, &cols[0], &values[0]);
    //        
    //        cols.clear();
    //        values.clear();
    //    }
    //}
    //
    ///* Assemble after setting the coefficients */
    //HYPRE_IJMatrixAssemble(A);
    //
    ///* Get the parcsr matrix object to use */
    //HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);
    
    /* Create the rhs and solution */
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    /* Set the rhs values to h^2 and the solution to zero */
    {
        double *rhs_values, *x_values;
        HYPRE_Int    *rows;

        rhs_values = (double*)calloc(local_size, sizeof(double));
        x_values = (double*)calloc(local_size, sizeof(double));
        rows = (HYPRE_Int*)calloc(local_size, sizeof(HYPRE_Int));

        for (i = 0; i<local_size; i++)
        {
            rhs_values[i] = bin[i];
            x_values[i] = 0.0;
            rows[i] = ilower + i;
        }

        HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);
        HYPRE_IJVectorSetValues(x, local_size, rows, x_values);

        free(x_values);
        free(rhs_values);
        free(rows);
    }
    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **)&par_b);

    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **)&par_x);
    /* AMG */
    //{
    //    int num_iterations;
    //    double final_res_norm;

    //    /* Create solver */
    //    HYPRE_BoomerAMGCreate(&solver);

    //    /* Set some parameters (See Reference Manual for more parameters) */
    //    HYPRE_BoomerAMGSetPrintLevel(solver, 0);  /* no printout */
    //    HYPRE_BoomerAMGSetOldDefault(solver); /* Falgout coarsening with modified classical interpolaiton */
    //    HYPRE_BoomerAMGSetRelaxType(solver, 4);   /* G-S/Jacobi hybrid relaxation */
    //    HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   /* uses C/F relaxation */
    //    HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeps on each level */
    //    HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
    //    HYPRE_BoomerAMGSetTol(solver, 1e-4);      /* conv. tolerance */
    //    HYPRE_BoomerAMGSetMaxIter(solver, 50);    /* maximum number of iterations */

    //    /* Now setup and solve! */
    //    HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
    //    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

    //    /* Run info - needed logging turned on */
    //    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
    //    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    //    if (myid == 0)
    //    {
    //        printf("\n");
    //        printf("Iterations = %d\n", num_iterations);
    //        printf("Final Relative Residual Norm = %e\n", final_res_norm);
    //        printf("\n");
    //    }

    //    /* Output the final solution */
    //    HYPRE_Complex v;
    //    for (i = ilower; i <= iupper; i++){
    //        HYPRE_IJVectorGetValues(x, 1, &i, &v);
    //        solution[i] = v;
    //    }

    //    /* Destroy solver */
    //    HYPRE_BoomerAMGDestroy(solver);

    //}

    /* PCG with AMG preconditioner */
    //{
    //    int num_iterations;
    //    double final_res_norm;

    //    /* Create solver */
    //    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

    //    /* Set some parameters (See Reference Manual for more parameters) */
    //    HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
    //    HYPRE_PCGSetTol(solver, 1e-4); /* conv. tolerance */
    //    HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
    //    HYPRE_PCGSetPrintLevel(solver, 0); /* print solve info */
    //    HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

    //    /* Now set up the AMG preconditioner and specify any parameters */
    //    HYPRE_BoomerAMGCreate(&precond);
    //    HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
    //    HYPRE_BoomerAMGSetCoarsenType(precond, 6);
    //    HYPRE_BoomerAMGSetOldDefault(precond);
    //    HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
    //    HYPRE_BoomerAMGSetNumSweeps(precond, 1);
    //    HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
    //    HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

    //    /* Set the PCG preconditioner */
    //    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
    //        (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

    //    /* Now setup and solve! */
    //    HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
    //    HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

    //    /* Run info - needed logging turned on */
    //    HYPRE_PCGGetNumIterations(solver, &num_iterations);
    //    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    //    if (myid == 0)
    //    {
    //        printf("\n");
    //        printf("Iterations = %d\n", num_iterations);
    //        printf("Final Relative Residual Norm = %e\n", final_res_norm);
    //        printf("\n");
    //    }

    //    /* Output the final solution */
    //    HYPRE_Complex v;
    //    for (i = ilower; i <= iupper; i++){
    //        HYPRE_IJVectorGetValues(x, 1, &i, &v);
    //        solution[i] = v;
    //    }

    //    /* Destroy solver and preconditioner */
    //    HYPRE_ParCSRPCGDestroy(solver);
    //    HYPRE_BoomerAMGDestroy(precond);
    //}
    

    

    /* Flexible GMRES with AMG Preconditioner */
    {
        HYPRE_Int    num_iterations;
        double final_res_norm;
        HYPRE_Int    restart = 10;
        HYPRE_Int    modify = 1;

        
        /* Create solver */
        HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_FlexGMRESSetKDim(solver, restart);
        HYPRE_FlexGMRESSetMaxIter(solver, 1000); /* max iterations */
        HYPRE_FlexGMRESSetTol(solver, 1e-4); /* conv. tolerance */
        HYPRE_FlexGMRESSetPrintLevel(solver, 0); /* print solve info */
        HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */


        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 6);
        HYPRE_BoomerAMGSetOldDefault(precond);
        HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
        HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

        /* Set the FlexGMRES preconditioner */
        HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
            (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

        
        if (modify)
            /* this is an optional call  - if you don't call it, hypre_FlexGMRESModifyPCDefault
            is used - which does nothing.  Otherwise, you can define your own, similar to
            the one used here */
            HYPRE_FlexGMRESSetModifyPC(solver,
            (HYPRE_PtrToModifyPCFcn)hypre_FlexGMRESModifyPCAMGExample);

        
        /* Now setup and solve! */
        HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
        HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);

        /* Run info - needed logging turned on */
        HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
        
        if (myid == 0)
        {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }

        /* Output the final solution */
        HYPRE_Complex v;
        for (i = ilower; i <= iupper; i++){
            HYPRE_IJVectorGetValues(x, 1, &i, &v);
            solution[i] = v;
        }

        /* Destory solver and preconditioner */
        HYPRE_ParCSRFlexGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);

    }
    
    /* Clean up */
    //HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
    

    return(0);


}

int hypre_FlexGMRESModifyPCAMGExample(void *precond_data, HYPRE_Int iterations,
    double rel_residual_norm)
{


    if (rel_residual_norm > .1)
    {
        HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver)precond_data, 10);
    }
    else
    {
        HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver)precond_data, 1);
    }


    return 0;
}