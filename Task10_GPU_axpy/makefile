NVCC = nvcc

CXXFLAGS  = -std=c++17 -O2 $(shell gsl-config --cflags)
LDFLAGS   = $(shell gsl-config --libs)
OMPFLAGS  = -Xcompiler -fopenmp 

TARGET = daxpy
SRC    = daxpy.cu

all: $(TARGET)

$(TARGET): $(SRC)
	$(NVCC) $(CXXFLAGS) $(OMPFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET)
