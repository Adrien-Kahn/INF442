all: cloud-tests cloud-run

.PHONY: all clean

point.o: point.cpp point.hpp
	g++ -c point.cpp

cloud.o: cloud.cpp point.hpp cloud.hpp
	g++ -c cloud.cpp

MyArea.o: MyArea.cpp MyArea.hpp
	g++ -std=c++11 -c MyArea.cpp `pkg-config gtkmm-3.0 --cflags`

cloud-tests: cloud-tests.cpp point.o cloud.o 
	g++ -o cloud-tests cloud-tests.cpp point.o cloud.o

cloud-run: cloud-run.cpp MyArea.o cloud.o point.o
	g++ -std=c++11 -o cloud-run cloud-run.cpp *.o `pkg-config gtkmm-3.0 --cflags --libs`

clean:
	rm -vf cloud-tests cloud-run *.o *~
