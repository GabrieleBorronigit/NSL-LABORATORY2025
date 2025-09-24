CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: main.exe


main.exe : main.o random.o functions.o
	$(CC) functions.o random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o :	functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)

run: main.exe 
	./main.exe
clean :
	rm *.o main.exe seed.out
