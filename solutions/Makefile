CFLAGS=-xHost -O3
MPICXX=mpicxx
MPICC=mpicc
OBJS=printarr_par.o
BINS=stencil stencil_blocked stencil_mpi stencil_mpi_ddt stencil_mpi_ddt_overlap stencil_mpi_ddt_overlap_carttopo \
stencil_mpi_ddt_overlap_carttopo_neighcolls \
stencil_mpi_carttopo_neighcolls stencil_mpi_overlap_carttopo_neighcolls \
stencil_mpi_shmem shmtest stencil_mpi_ddt_rma dsa-gups

all: $(BINS)

%.o: %.cpp Makefile
	$(MPICXX) $(CFLAGS) $< -c -o $@

stencil: stencil.cpp Makefile
	$(MPICXX) $(CFLAGS) -o $@ $<

stencil_blocked: stencil_blocked.cpp Makefile
	$(MPICXX) $(CFLAGS) -o $@ $<

stencil_mpi: stencil_mpi.cpp Makefile $(OBJS)
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

stencil_mpi_carttopo_neighcolls: stencil_mpi_carttopo_neighcolls.cpp Makefile $(OBJS) 
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

stencil_mpi_overlap_carttopo_neighcolls: stencil_mpi_overlap_carttopo_neighcolls.cpp Makefile $(OBJS) 
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

stencil_mpi_ddt_overlap_carttopo_neighcolls: stencil_mpi_ddt_overlap_carttopo_neighcolls.cpp Makefile $(OBJS) 
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

stencil_mpi_ddt: stencil_mpi_ddt.cpp Makefile $(OBJS)
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

stencil_mpi_ddt_overlap: stencil_mpi_ddt_overlap.cpp Makefile $(OBJS)
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

stencil_mpi_ddt_overlap_carttopo: stencil_mpi_ddt_overlap_carttopo.cpp Makefile $(OBJS)
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

stencil_mpi_shmem: stencil_mpi_shmem.cpp Makefile $(OBJS)
	$(MPICXX) $(CFLAGS) -o $@ $<  $(OBJS)

shmtest: shmtest.cpp Makefile $(OBJS)
	$(MPICXX) $(CFLAGS) -o $@ $<  $(OBJS)

stencil_mpi_ddt_rma: stencil_mpi_ddt_rma.cpp Makefile $(OBJS)
	$(MPICXX) $(CFLAGS) -o $@ $< $(OBJS)

dsa-gups: dsa-gups.c Makefile $(OBJS)
	$(MPICC) $(CFLAGS) -o $@ $< $(OBJS)

clean:
	rm -f $(BINS) $(OBJS)
	rm -f output*.bmp
