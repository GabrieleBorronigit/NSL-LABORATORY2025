CXX = mpic++
CXXFLAGS = -std=c++17 -O2
AFLAGS = -larmadillo
LDFLAGS = -Wl,-rpath,/usr/local/lib
OBJ = TSP.o random.o Genetic_algorithm.o
DEPS = random.h Genetic_algorithm.h

all: tsp

tsp: $(OBJ)
	$(CXX) $(CXXFLAGS) -o tsp $(OBJ) $(LDFLAGS) $(AFLAGS)

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o tsp
