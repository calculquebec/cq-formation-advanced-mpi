// -*- Mode: C; -*-
//-----------------------------------------------------------------
// Data distribution of matrix: checkerboard block
// based on example by
// Daniel R. Reynolds
// SMU, Mathematics
// Math 6495
// 7 January 2009
// and Michael J. Quinn, Parallel programming in C with MPI and OpenMP
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
  int m, n, i, is, ie, j, js, je, numprocs, myid, repeat;
  MPI_Comm col_comm, row_comm, grid_comm;
  int *counts, *displs, grid_id, grid_size[2], grid_coords[2], periodic[2];
  double **A, *x, *myv, *v, *v_part;
  double norm2, stime, ftime;

  //======= Internals ============
  
  // initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  grid_size[0] = grid_size[1] = sqrt(numprocs);
  if (grid_size[0]*grid_size[1] != numprocs) {
     printf("Number of processors %d not square\n", numprocs);
     MPI_Abort(MPI_COMM_WORLD, 1);
     exit(1);
  }
    
  // input the size of the system
  if (myid == 0) {
     printf("%s\n", "We will multiply an m*n matrix by an n-vector");
     scanf("%d %d", &m, &n);
     printf("m = %d, n = %d\n",m, n);
  }

  periodic[0] = periodic[1] = 0;
  MPI_Cart_create (MPI_COMM_WORLD, 2, grid_size, periodic, 1, &grid_comm);
  MPI_Comm_rank (grid_comm, &grid_id);
  MPI_Cart_coords (grid_comm, grid_id, 2, grid_coords);
  MPI_Comm_split (grid_comm, grid_coords[0], grid_coords[1], &row_comm);
  MPI_Comm_split (grid_comm, grid_coords[1], grid_coords[0], &col_comm);

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
  is = grid_coords[0]*m/grid_size[0];
  ie = (grid_coords[0]+1)*m/grid_size[0];
  js = grid_coords[1]*n/grid_size[1];
  je = (grid_coords[1]+1)*n/grid_size[1];

  // initialize the matrix and vectors 
  // (only store matrix & vector parts that you need on this proc)
  A = dmatrix(ie-is, je-js);
  for(i=0; i<ie-is; i++) {
     for(j=0; j<je-js; j++)
        A[i][j] = 1.0/(1+(i+is-(j+js))*(i+is-(j+js)));
  }
  x = malloc((je-js)*sizeof(*x));
  myv = malloc((ie-is)*sizeof(*myv));

  /* initialize x and v for the first column only */
  v = NULL;
  if (grid_coords[1] == 0) {
     for(j=0; j<je-js; j++)
       x[j] = 1.0;
     v = malloc((ie-is)*sizeof(*v));
  }

  // start timer
  MPI_Barrier(MPI_COMM_WORLD);
  stime = MPI_Wtime();

  for(repeat=0; repeat<10000; repeat++) {

     // Step 1: send elements from first column processes to first row processes
     if (grid_coords[1] == 0 && grid_coords[0] > 0) {
        int dest, coords[2];
        coords[0] = 0;
        coords[1] = grid_coords[0];
        MPI_Cart_rank(grid_comm, coords, &dest);
        MPI_Send(x, je-js, MPI_DOUBLE, dest, 0, grid_comm);
     }

     if (grid_coords[0] == 0 && grid_coords[1] > 0) {
        int src, coords[2];
        MPI_Status status;
        coords[0] = grid_coords[1];
        coords[1] = 0;
        MPI_Cart_rank(grid_comm, coords, &src);
        MPI_Recv(x, je-js, MPI_DOUBLE, src, MPI_ANY_TAG, grid_comm, &status);
     }

     // Step 2: Broadcast along columns
     MPI_Bcast(x, je-js, MPI_DOUBLE, 0, col_comm);

     // compute matrix-vector product
     for(i=0; i<ie-is; i++) {
        double ip = 0.0;
        for(j=0; j<je-js; j++) {
           ip += A[i][j] * x[j];
        }
        myv[i] = ip;
     }
     MPI_Reduce(myv, v, ie-is, MPI_DOUBLE, MPI_SUM, 0, row_comm);
  }

  // stop timer
  MPI_Barrier(MPI_COMM_WORLD);
  ftime = MPI_Wtime();

  // output 2-norm of product and runtime to screen
  if (grid_coords[1] == 0) {
     double mynorm2 = 0;
     for (i=0; i<ie-is; i++)
        mynorm2 += v[i] * v[i];
     MPI_Reduce(&mynorm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, 0, col_comm);
     if (grid_coords[0] == 0) {
        printf("rank=%d = (%d,%d), 2-norm of product = %g\n",myid,grid_coords[0],grid_coords[1],sqrt(norm2));
     }
  }

  if (myid == 0) {
     printf("runtime = %g\n",ftime-stime);
  }

  // finalize MPI
  MPI_Finalize();

  // free vectors
  free(*A); free(A); free(x); free(v);
  return 0;
  }
//=================================================================
