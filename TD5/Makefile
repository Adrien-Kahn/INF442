CPPFLAGS=-O3 -std=c++11 `pkg-config gtkmm-3.0 --cflags --libs`

all: display test_radial test_flat estimate_flat test_gaussian estimate_gaussian test_knn estimate_knn meanshift

display: point.o cloud.o DisplayArea.o display.o
	g++ -o $@ $^ $(CPPFLAGS)

test_radial: point.o cloud.o kernel.o radial.o test_radial.o
	g++ -o $@ $^ $(CPPFLAGS)

test_flat: point.o cloud.o kernel.o radial.o flat.o test_flat.o
	g++ -o $@ $^ $(CPPFLAGS)

estimate_flat: point.o cloud.o kernel.o radial.o flat.o KernelArea.o estimate_flat.o
	g++ -o $@ $^ $(CPPFLAGS)

test_gaussian: point.o cloud.o kernel.o radial.o gaussian.o test_gaussian.o
	g++ -o $@ $^ $(CPPFLAGS)

estimate_gaussian: point.o cloud.o kernel.o radial.o gaussian.o KernelArea.o estimate_gaussian.o
	g++ -o $@ $^ $(CPPFLAGS)

test_knn: point.o cloud.o kernel.o knn.o test_knn.o
	g++ -o $@ $^ $(CPPFLAGS)

estimate_knn: point.o cloud.o kernel.o knn.o KernelArea.o estimate_knn.o
	g++ -o $@ $^ $(CPPFLAGS)

meanshift: meanshift.cpp
	g++ -o $@ $^ $(CPPFLAGS)

%.o: %.cpp
	g++ -c -o $@ $^ $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f display test_radial test_flat estimate_flat test_gaussian estimate_gaussian test_knn estimate_knn meanshift

.PHONY: all clean
