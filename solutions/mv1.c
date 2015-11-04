// -*- Mode: C; -*-
//-----------------------------------------------------------------
// Data distribution of matrix: rowwise block striped
// based on example by
// Daniel R. Reynolds
// SMU, Mathematics
// Math 6495
// 7 January 2009
//=================================================================

//======= Inclusions ===========
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double** dmatrix(int n_rows, int n_cols)
{
  int i, total_pts;
  double **m;

  if (n_rows<=0 || n_cols<=0) return(NULL);
  total_pts = n_rows*n_cols;
  
  if ( (m = malloc( n_rows * sizeof(*m))) == NULL) {
    printf("dmatrix: memory allocation failure! %zd bytes\n",n_rows * sizeof(*m));
  } else {
    if ( (m[0] = malloc( total_pts * sizeof(*m[0]))) == NULL) {
      free(m);
      m = NULL;
      printf("dmatrix: memory allocation failure! %zd bytes\n",
	     total_pts * sizeof(*m[0]));
    }
    else
      for (i=1; i<n_rows; i++) m[i] = m[i-1] + n_cols;
  }
  return(m);
}

int main(int argc, char **argv)
{
  //-----------------------------------------------------------------
  // Description: 
  //    Computes the product of an m*n matrix and an n-vector.
  //-----------------------------------------------------------------

  //======= Declarations =========
  int m, n, i, is, ie, j, numprocs, myid, repeat;
  int *counts, *displs;
  double **A, *x, *myv, *v;
  double norm2, stime, ftime;

  //======= Internals ============
  
  // initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // input the size of the system
  if (myid == 0) {
     printf("%s\n", "We will multiply an m*n matrix by an n-vector");
     scanf("%d %d", &m, &n);
     printf("m = %d, n = %d",m, n);
  }

  // root node sends m and n out to other processors
  MPI_Bcast(&m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

  if ((m < 1) || (n < 1)) {
     if (myid == 0) {
        printf("illegal input, m = %d and n = %d must both be positive\n", m, n);
     }
     MPI_Abort(MPI_COMM_WORLD, 1);
     exit(1);
  }
  
  // root node outputs parallelism information to screen
  if (myid == 0) {
     printf("%s%d%s\n", " starting MPI with ", numprocs, " processes");
  }

  // determine this processor's interval
  is = myid*m/numprocs;
  ie = (myid+1)*m/numprocs;

  // initialize the matrix and vectors 
  // (only store matrix & vector parts that you need on this proc)
  A = dmatrix(ie-is, n);
  x = malloc(n*sizeof(*x));
  myv = malloc((ie-is)*sizeof(*myv));
  v = malloc(m*sizeof(*v));
  for(j=0; j<n; j++)
     x[j] = 1.0;
  for(i=0; i<ie-is; i++) {
     for(j=0; j<n; j++)
        A[i][j] = 1.0/(1+(i+is-j)*(i+is-j));
  }

  // all nodes collect result
  counts = malloc(numprocs*sizeof(*counts));
  displs = malloc(numprocs*sizeof(*displs));
  for (i=0; i<numprocs; i++) {
     displs[i] = i*m/numprocs;
     counts[i] = (i+1)*m/numprocs - displs[i];
  }

  // start timer
  MPI_Barrier(MPI_COMM_WORLD);
  stime = MPI_Wtime();

  for(repeat=0; repeat<10000; repeat++) {

     // compute matrix-vector product
     for(i=0; i<ie-is; i++) {
        double ip = 0.0;
        for(j=0; j<n; j++) {
	   ip += A[i][j]*x[j];
        }
        myv[i] = ip;
     }

     MPI_Allgatherv(myv, counts[myid], MPI_DOUBLE,
       v, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
  }

  // stop timer
  MPI_Barrier(MPI_COMM_WORLD);
  ftime = MPI_Wtime();

  // output 2-norm of product and runtime to screen
  norm2 = 0;
  for (i=0; i<m; i++)
    norm2 += v[i] * v[i];
  
  if (myid == 0) {
     printf("2-norm of product = %g\n",sqrt(norm2));
     printf("runtime = %g\n",ftime-stime);
  }

  // finalize MPI
  MPI_Finalize();

  // free vectors
  free(*A); free(A); free(x); free(v); free(counts); free(displs);
  return 0;
}
//=================================================================
