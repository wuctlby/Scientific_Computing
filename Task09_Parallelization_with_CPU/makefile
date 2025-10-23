CXX = g++
MPICXX = mpic++
CXXGSLFLAGS += `gsl-config --cflags` -std=c++17 -O2
LDGSLFLAGS  += `gsl-config --libs` -fopenmp
MPIFLAGS = -fopenmp

TARGETS = daxpy daxpy_MPI

All: $(TARGETS)

daxpy: daxpy.cxx
	${CXX} ${CXXGSLFLAGS} -o $@ $< ${LDGSLFLAGS}

daxpy_MPI: daxpy_MPI.cxx
	${MPICXX} ${CXXGSLFLAGS} ${MPIFLAGS} -o $@ $< ${LDGSLFLAGS}
