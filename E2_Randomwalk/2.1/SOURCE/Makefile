all: main

main: main.o random.o functions.o
	g++ main.o random.o functions.o -o main

main.o: main.cpp random.h functions.h
	g++ -std=c++11 -c main.cpp

RNG/random.o: random.cpp random.h
	g++ -std=c++11 -c random.cpp -o random.o

functions.o: functions.cpp functions.h
	g++ -std=c++11 -c functions.cpp

run: main
	./main

clean:
	rm -f *.o *.o main
