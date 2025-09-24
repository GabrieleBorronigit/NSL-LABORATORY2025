all: main

main: main.o RandomWalk.o random.o
	g++ main.o RandomWalk.o random.o -o main

main.o: main.cpp RandomWalk.h random.h
	g++ -std=c++11 -c main.cpp

RandomWalk.o: RandomWalk.cpp RandomWalk.h random.h
	g++ -std=c++11 -c RandomWalk.cpp

random.o: random.cpp random.h
	g++ -std=c++11 -c random.cpp

run: main
	./main

clean:
	rm -f *.o main
