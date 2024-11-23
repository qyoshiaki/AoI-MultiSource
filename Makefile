CXX = g++
CXXFLAGS = -O2 -Wall
SRC = main.cc matrix/matrix.cc matrix/myVector.cc distribution/distribution.cc
OBJ = main.o matrix/matrix.o matrix/myVector.o distribution/distribution.o

a.out: $(OBJ)
	$(CXX) -o $@ $(OBJ)

clean:
	rm a.out $(OBJ)

.cc.o:
	$(CXX) -o $@ -c $< $(CXXFLAGS)
