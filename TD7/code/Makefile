EIGEN=../../eigen-3.3.9
ANN_INCLUDE=../../ann_1.1.2/include
ANN_LIB=../../ann_1.1.2/lib

# EIGEN=/opt/local/include/eigen3
# ANN_HOME=/Users/elias/ET/et-data/teaching/X2019/ann_1.1.2/
# ANN_INCLUDE=$(ANN_HOME)/include
# ANN_LIB=$(ANN_HOME)/lib

CC=g++

all: test_dataset test_ex1q1 test_ex1q2 test_ex1q3 test_ex2q1 test_ex2q2 test_linear test_knn test_k_knn

INCLUDES=-I$(EIGEN) -I$(ANN_INCLUDE)
CXXFLAGS=$(INCLUDES) -std=c++11 -g -O2 -Wall -Wextra

LDFLAGS=
LDFLAGS_ANN=-L$(ANN_LIB) -lANN

%: %.o
	$(CC)  $(CXXFLAGS) $^ $(LDFLAGS) -o $@

Dataset.o: Dataset.cpp Dataset.hpp
KnnRegression.o: KnnRegression.cpp KnnRegression.hpp
LinearRegression.o: LinearRegression.cpp LinearRegression.hpp
Regression.o: Regression.cpp Regression.hpp Dataset.hpp
test_dataset.o: test_dataset.cpp LinearRegression.hpp
test_k_knn.o: test_k_knn.cpp Dataset.hpp
test_knn.o: test_knn.cpp Dataset.hpp
test_linear.o: test_linear.cpp LinearRegression.hpp
test_ex1q1.o: test_ex1q1.cpp LinearRegression.hpp
test_ex1q2.o: test_ex1q2.cpp LinearRegression.hpp
test_ex1q3.o: test_ex1q3.cpp LinearRegression.hpp
test_ex2q1.o: test_ex2q1.cpp
test_ex2q2.o: test_ex2q2.cpp Dataset.hpp

test_dataset: test_dataset.o Dataset.o
test_linear: test_linear.o LinearRegression.o Regression.o Dataset.o
test_ex1q1: test_ex1q1.o LinearRegression.o Regression.o Dataset.o
test_ex1q2: test_ex1q2.o LinearRegression.o Regression.o Dataset.o
test_ex1q3: test_ex1q3.o LinearRegression.o Regression.o Dataset.o

test_ex2q1: test_ex2q1.o KnnRegression.o Regression.o Dataset.o
	$(CC) $^ $(LDFLAGS_ANN) -o $@
test_ex2q2: test_ex2q2.o KnnRegression.o Regression.o Dataset.o
	$(CC) $^ $(LDFLAGS_ANN) -o $@

test_knn_constructor: test_knn_constructor.o KnnRegression.o Regression.o Dataset.o
	$(CC) $^ $(LDFLAGS_ANN) -o $@
test_knn: test_knn.o KnnRegression.o Regression.o Dataset.o
	$(CC) $^ $(LDFLAGS_ANN) -o $@
test_k_knn: test_k_knn.o KnnRegression.o Regression.o Dataset.o
	$(CC) $^ $(LDFLAGS_ANN) -o $@

clean:
	rm -f *.o test_ex1q2 test_ex1q3 test_ex1q1 test_ex2q1 test_ex2q2 
	rm -f test_linear test_knn test_k_knn test_dataset output.txt

.PHONY: all clean
