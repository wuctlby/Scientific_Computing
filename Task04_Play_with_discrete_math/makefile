CXX = g++
# CXXFLAGS = -I/usr/include/hdf5/serial -std=c++17 -O2
# LDFLAGS  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5
CXXGSLFLAGS += `gsl-config --cflags` -std=c++17 -O2
LDGSLFLAGS  += `gsl-config --libs`

# DEPS = utils.h
# DEPS_GSL = utils_gsl.h

# OBJ_CALC = calculator.o
# OBJ_VEC = vector_generator.o

OBJ = function.o

# TARGETS = calculator vector_generator
# TARGETS_GSL = calculator_gsl vector_generator_gsl

TARGETS = function

function: function.cxx
	${CXX} ${CXXGSLFLAGS} -o $@ $< ${LDGSLFLAGS}

.PHONY: hdf5 clean gsl

hdf5: $(TARGETS)

calculator: calculator.cxx $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

vector_generator: vector_generator.cxx $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

gsl: $(TARGETS_GSL)

calculator_gsl: calculator_gsl.cxx $(DEPS_GSL)
	$(CXX) $(CXXGSLFLAGS) -o $@ $< $(LDGSLFLAGS)

vector_generator_gsl: vector_generator_gsl.cxx $(DEPS_GSL)
	$(CXX) $(CXXGSLFLAGS) -o $@ $< $(LDGSLFLAGS)

clean:
	rm -f $(TARGETS) $(TARGETS_GSL)