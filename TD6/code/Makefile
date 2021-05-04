EIGEN=../../eigen-3.3.9
ANN_INCLUDE=../../ann_1.1.2/include
ANN_LIB=../../ann_1.1.2/lib

all: test_dataset test_1q1 test_1q2 test_2q1 test_3q1 test_3q2 test_knn test_random_projection test_roc_curve

Dataset.o: Dataset.cpp
	g++ -c -std=c++11 -O2 -Wall -Wextra Dataset.cpp

Classification.o: Classification.cpp
	g++ -c -std=c++11 -O2 -I$(ANN_INCLUDE) -Wall -Wextra Classification.cpp
	
KnnClassification.o: KnnClassification.cpp
	g++ -c -std=c++11 -O2 -I$(ANN_INCLUDE) -Wall -Wextra KnnClassification.cpp	

ConfusionMatrix.o: ConfusionMatrix.cpp
	g++ -c -std=c++11 -O2 -I$(ANN_INCLUDE) -Wall -Wextra ConfusionMatrix.cpp	

RandomProjection.o: RandomProjection.cpp
	g++ -c -std=c++11 -O2 -I$(EIGEN) -std=c++11 -Wall -Wextra RandomProjection.cpp	

test_dataset: test_dataset.cpp Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -Wall -Wextra test_dataset.cpp Dataset.o -o test_dataset

test_1q1: test_1q1.cpp KnnClassification.o Classification.o Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -I$(ANN_INCLUDE) -Wall -Wextra -Wno-unused-value test_1q1.cpp KnnClassification.o Classification.o Dataset.o -o test_1q1 -L$(ANN_LIB) -lANN 

test_1q2: test_1q2.cpp KnnClassification.o Classification.o Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -I$(ANN_INCLUDE) -Wall -Wextra -Wno-unused-value test_1q2.cpp KnnClassification.o Classification.o Dataset.o -o test_1q2 -L$(ANN_LIB) -lANN 

test_2q1: test_2q1.cpp ConfusionMatrix.o
	g++ -std=c++11 -O2 -I$(EIGEN) -Wall -Wextra -Wno-unused-value test_2q1.cpp ConfusionMatrix.o -o test_2q1

test_3q1: test_3q1.cpp RandomProjection.o Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -Wall -Wextra -Wno-unused-value test_3q1.cpp RandomProjection.o Dataset.o -o test_3q1

test_3q2: test_3q2.cpp RandomProjection.o Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -Wall -Wextra -Wno-unused-value test_3q2.cpp RandomProjection.o Dataset.o -o test_3q2

test_knn: test_knn.cpp ConfusionMatrix.o KnnClassification.o Classification.o Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -I$(ANN_INCLUDE) -Wall -Wextra test_knn.cpp ConfusionMatrix.o KnnClassification.o Classification.o Dataset.o -o test_knn -L$(ANN_LIB) -lANN 

test_random_projection: test_random_projection.cpp RandomProjection.o ConfusionMatrix.o KnnClassification.o Classification.o Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -I$(ANN_INCLUDE) -Wall -Wextra test_random_projection.cpp RandomProjection.o ConfusionMatrix.o KnnClassification.o Classification.o Dataset.o -o test_random_projection -L$(ANN_LIB) -lANN

test_roc_curve: test_roc_curve.cpp ConfusionMatrix.o KnnClassification.o Classification.o Dataset.o
	g++ -std=c++11 -O2 -I$(EIGEN) -I$(ANN_INCLUDE) -Wall -Wextra test_roc_curve.cpp ConfusionMatrix.o KnnClassification.o Classification.o Dataset.o -o test_roc_curve -L$(ANN_LIB) -lANN

.PHONY: all clean

clean:
	rm -f *.o
	rm -f test_dataset
	rm -f test_1q1
	rm -f test_1q2
	rm -f test_2q1
	rm -f test_3q1
	rm -f test_3q2
	rm -f test_knn
	rm -f test_random_projection
	rm -f test_roc_curve
