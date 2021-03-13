CXX = g++
CXXFLAGS = --std=c++11

TESTS = Test3 Test4 Test5 Test6
OBJS = ArrayToMatrix.o TestInfra.o

.PHONY: all clean

all: $(TESTS)

libraries: $(OBJS)

Test3: Stats.cpp $(OBJS)
Test4: Stats.cpp $(OBJS)
Test5: Stats.cpp $(OBJS)
Test6: Stats.cpp $(OBJS)

clean:
	rm -f $(TESTS) $(OBJS)
	rm -f *~
