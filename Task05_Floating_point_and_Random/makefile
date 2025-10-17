CXX = g++
CXXGSLFLAGS += `gsl-config --cflags` -std=c++17 -O2
LDGSLFLAGS  += `gsl-config --libs`

TARGETS = $(patsubst %.cxx, %, $(wildcard *.cxx))

All: $(TARGETS)

%: %.cxx
	${CXX} ${CXXGSLFLAGS} -o $@ $< ${LDGSLFLAGS}
