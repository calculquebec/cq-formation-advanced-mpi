/*
 * Copyright (c) 2012 Torsten Hoefler. All rights reserved.
 *
 * Author(s): Torsten Hoefler <htor@illinois.edu>
 *
 */

#include "stencil_par.h"

int main(int argc, char **argv) {

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided); 
  int r,p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &r);
  MPI_Comm_size(comm, &p);
  int n, energy, niters, px, py;

  if (r==0) {
      // argument checking
      if(argc < 6) {
          if(!r) printf("usage: stencil_mpi <n> <energy> <niters> <px> <py>\n");
          MPI_Finalize();
          exit(1);       
      }
      
      n = atoi(argv[1]); // nxn grid
      energy = atoi(argv[2]); // energy to be injected per iteration
      niters = atoi(argv[3]); // number of iterations
      px=atoi(argv[4]); // 1st dim processes
      py=atoi(argv[5]); // 2nd dim processes
      
      if(px * py != p) MPI_Abort(comm, 1);// abort if px or py are wrong
      if(n % py != 0) MPI_Abort(comm, 2); // abort px needs to divide n
      if(n % px != 0) MPI_Abort(comm, 3); // abort py needs to divide n
      
      // distribute arguments
      int args[5] = {n, energy, niters, px,  py};
      MPI_Bcast(args, 5, MPI_INT, 0, comm);
  }
  else {
      int args[5];
      MPI_Bcast(args, 5, MPI_INT, 0, comm);
      n=args[0]; energy=args[1]; niters=args[2]; px=args[3]; py=args[4];
  }

  // determine my coordinates (x,y) -- r=x*a+y in the 2d processor array
  int rx = r % px;
  int ry = r / px;
  // determine my four neighbors
  int north = (ry-1)*px+rx; if(ry-1 < 0)   north = MPI_PROC_NULL;
  int south = (ry+1)*px+rx; if(ry+1 >= py) south = MPI_PROC_NULL;
  int west= ry*px+rx-1;     if(rx-1 < 0)   west = MPI_PROC_NULL;
  int east = ry*px+rx+1;    if(rx+1 >= px) east = MPI_PROC_NULL;
  // decompose the domain
  int bx = n/px; // block size in x
  int by = n/py; // block size in y
  int offx = rx*bx; // offset in x
  int offy = ry*by; // offset in y

  //printf("%i (%i,%i) - w: %i, e: %i, n: %i, s: %i\n", r, ry,rx,west,east,north,south);

  // allocate two work arrays
  double *aold = (double*)calloc(1,(bx+2)*(by+2)*sizeof(double)); // 1-wide halo zones!
  double *anew = (double*)calloc(1,(bx+2)*(by+2)*sizeof(double)); // 1-wide halo zones!
  double *tmp;

  // initialize three heat sources
  const int nsources=3;
  int sources[nsources][2] = {{n/2,n/2}, {n/3,n/3}, {n*4/5,n*8/9}};
  int locnsources=0; // number of sources in my area
  int locsources[nsources][2]; // sources local to my rank
  for (int i=0; i<nsources; ++i) { // determine which sources are in my patch
    int locx = sources[i][0] - offx;
    int locy = sources[i][1] - offy;
    if(locx >= 0 && locx < bx && locy >= 0 && locy < by) {
      locsources[locnsources][0] = locx+1; // offset by halo zone
      locsources[locnsources][1] = locy+1; // offset by halo zone
      locnsources++;
    }
  }

  double t=-MPI_Wtime(); // take time
  // allocate communication buffers
  double *sbufnorth = (double*)calloc(1,bx*sizeof(double)); // send buffers
  double *sbufsouth = (double*)calloc(1,bx*sizeof(double));
  double *sbufeast = (double*)calloc(1,by*sizeof(double));
  double *sbufwest = (double*)calloc(1,by*sizeof(double));
  double *rbufnorth = (double*)calloc(1,bx*sizeof(double)); // receive buffers
  double *rbufsouth = (double*)calloc(1,bx*sizeof(double));
  double *rbufeast = (double*)calloc(1,by*sizeof(double));
  double *rbufwest = (double*)calloc(1,by*sizeof(double));

  double heat; // total heat in system
  for(int iter=0; iter<niters; ++iter) {
    // refresh heat sources
    for(int i=0; i<locnsources; ++i) {
      aold[ind(locsources[i][0],locsources[i][1])] += energy; // heat source
    }

    // exchange data with neighbors
    for(int i=0; i<bx; ++i) sbufnorth[i] = aold[ind(i+1,1)]; // pack loop - last valid region
    for(int i=0; i<bx; ++i) sbufsouth[i] = aold[ind(i+1,by)]; // pack loop
    for(int i=0; i<by; ++i) sbufeast[i] = aold[ind(bx,i+1)]; // pack loop
    for(int i=0; i<by; ++i) sbufwest[i] = aold[ind(1,i+1)]; // pack loop
/*       int MPI_Sendrecv(void *sbuf, int scount, MPI_Datatype sdtype,
                       int dest, int stag, void *rbuf, int rcount,
                       MPI_Datatype rdtype, int src, int rtag,
                       MPI_Comm comm, MPI_Status *status)*/

    MPI_Sendrecv(sbufnorth, bx, MPI_DOUBLE, north, 9,
                 rbufsouth, bx, MPI_DOUBLE, south, 9, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sbufsouth, bx, MPI_DOUBLE, south, 9,
                 rbufnorth, bx, MPI_DOUBLE, north, 9, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sbufeast, by, MPI_DOUBLE, east, 9,
                 rbufwest, by, MPI_DOUBLE, west, 9, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sbufwest, by, MPI_DOUBLE, west, 9,
                 rbufeast, by, MPI_DOUBLE, east, 9, comm, MPI_STATUS_IGNORE);
    for(int i=0; i<bx; ++i) aold[ind(i+1,0)] = rbufnorth[i]; // unpack loop - into ghost cells
    for(int i=0; i<bx; ++i) aold[ind(i+1,by+1)] = rbufsouth[i]; // unpack loop
    for(int i=0; i<by; ++i) aold[ind(bx+1,i+1)] = rbufeast[i]; // unpack loop
    for(int i=0; i<by; ++i) aold[ind(0,i+1)] = rbufwest[i]; // unpack loop

    // update grid points
    heat = 0.0;
#pragma omp parallel for reduction(+:heat)
    for(int j=1; j<by+1; ++j) {
      for(int i=1; i<bx+1; ++i) {
        anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
        heat += anew[ind(i,j)];
      }
    }

    // swap arrays
    tmp=anew; anew=aold; aold=tmp;

    // optional - print image
    //if(iter == niters-1) printarr_par(iter, anew, n, px, py, rx, ry, bx, by, offx, offy, comm);
  }
  t+=MPI_Wtime();

  // get final heat in the system
  double rheat;
  MPI_Allreduce(&heat, &rheat, 1, MPI_DOUBLE, MPI_SUM, comm);
  if(!r) printf("[%i] last heat: %f time: %f\n", r, rheat, t);

  MPI_Finalize();
}
