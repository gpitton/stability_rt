CXX=g++
CXXFLAGS=-O3 -std=gnu++17 -Wall -Wextra -Wpedantic
EIGEN3_INCLUDE_DIR=/usr/include/eigen3
INCLUDES=-I${EIGEN3_INCLUDE_DIR}
LDFLAGS=-fopenmp


all: stability two_parameters


stability: stability.cc


two_parameters: two_parameters.cc


clean:
	rm stability two_parameters


.cc:
	$(CXX) $(CXXFLAGS) $@.cc $(INCLUDES) $(LDFLAGS) -o $@
