CXX = g++
CXXFFTWFLAGS += -I/usr/include -std=c++17 -O2
LDFFTWFLAGS  += -lfftw3

TARGETS = $(patsubst %.cxx, %, $(wildcard *.cxx))

All: $(TARGETS)

%: %.cxx
	${CXX} ${CXXFFTWFLAGS} -o $@ $< ${LDFFTWFLAGS}
