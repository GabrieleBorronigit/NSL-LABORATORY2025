CXX = g++
CXXFLAGS = -std=c++17 -O2
AFLAGS = -larmadillo
LDFLAGS = -Wl,-rpath,/usr/local/lib
OBJ = TSP.o random.o ACO.o
DEPS = random.h ACO.h

all: tsp

tsp: $(OBJ)
	$(CXX) $(CXXFLAGS) -o tsp $(OBJ) $(LIBS) $(LDFLAGS) $(AFLAGS)

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o tsp
