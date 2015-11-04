/*
 * (C) 2012 by University of Chicago.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

typedef struct {
  double *base;
  int     n;
  int     rank, nproc;
  MPI_Win win;
} dsa_t;

const int N    = 10000;
const int NUPD = 100000;


/** Create a distributed array.
  */
dsa_t *DSA_Create(MPI_Comm comm, int n) {
  int   i;
  dsa_t *dsa = malloc(sizeof(dsa_t));

  dsa->n = n;

  MPI_Comm_rank(comm, &dsa->rank);
  MPI_Comm_size(comm, &dsa->nproc);

  MPI_Alloc_mem(n*sizeof(double), MPI_INFO_NULL, &dsa->base);

  for (i = 0; i < n; i++)
    dsa->base[i] = 0.0;

  MPI_Win_create(dsa->base, n*sizeof(double), sizeof(double), 
      MPI_INFO_NULL, comm, &dsa->win);
  
  return dsa;
}


/** Destroy a distributed array.
  */
void DSA_Destroy(dsa_t *dsa) {
  MPI_Win_free(&dsa->win);
  MPI_Free_mem(dsa->base);
  
  free(dsa);
}


/** Read a location in the distributed array.
  */
double DSA_Get(dsa_t *dsa, int idx) {
  int rank, offset;
  double value;

  rank = idx / dsa->n;
  offset = idx % dsa->n;

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, dsa->win);
  if (rank == dsa->rank) {
    value = dsa->base[offset];
  } else {
    MPI_Get(&value, 1, MPI_DOUBLE, rank, offset, 1, MPI_DOUBLE, dsa->win);
  }
  MPI_Win_unlock(rank, dsa->win);

  return value;
}


/** Read a location in the distributed array.
  */
void DSA_Put(dsa_t *dsa, int idx, double value) {
  int rank, offset;

  rank = idx / dsa->n;
  offset = idx % dsa->n;

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, dsa->win);
  if (rank == dsa->rank) {
    dsa->base[offset] = value;
  } else {
    MPI_Put(&value, 1, MPI_DOUBLE, rank, offset, 1, MPI_DOUBLE, dsa->win);
  }
  MPI_Win_unlock(rank, dsa->win);
}


int main(int argc, char **argv) {
  int    i, rank, nproc;
  double t;
  dsa_t *dsa;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  dsa = DSA_Create(MPI_COMM_WORLD, N);

  MPI_Barrier(MPI_COMM_WORLD);
  t = MPI_Wtime();

  for (i = 0; i < NUPD; i++) {
    int    idx;
    double v;

    idx = rand() % (N*nproc);

    v   = DSA_Get(dsa, idx);
    v  += rank;

    DSA_Put(dsa, idx, v);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  t = MPI_Wtime() - t;

  if (rank == 0) 
    printf("%d updates in %f sec (%f GUPS)\n", NUPD*nproc, t, (NUPD*nproc)/1.0e9/t);

  DSA_Destroy(dsa);

  MPI_Finalize();

  return 0;
}
