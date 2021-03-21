
CXX = g++
CXXFLAGS = -g --std=c++11 
OBJS = 
CTESTS = test_distance test_linear_scan ask_iris test_median test_partition \
	test_create_node test_defeatist test_backtracking compare_search_algos

all: $(CTESTS)

test_distance: test_distance.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -o $@ $@.cpp 

test_linear_scan: test_linear_scan.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -o $@ $@.cpp 


ask_iris: ask_iris.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -o $@ $@.cpp 


test_median: test_median.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -o $@ $@.cpp 

test_partition: test_partition.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -o $@ $@.cpp 

test_create_node: test_create_node.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -o $@ $@.cpp 

test_defeatist: test_defeatist.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -o $@ $@.cpp 

test_backtracking: test_backtracking.cpp retrieval.hpp
	${CXX} ${CXXFLAGS}  -o $@ $@.cpp 

compare_search_algos: compare_search_algos.cpp retrieval.hpp
	${CXX} ${CXXFLAGS} -O3  -o $@ $@.cpp 

clean:
	rm -rf *.o $(CTESTS) *.dSYM *~
