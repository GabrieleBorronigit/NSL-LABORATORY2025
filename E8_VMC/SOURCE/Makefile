CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11
OBJS = main.o VariationalMC.o random.o
TARGET = vmc

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

main.o: main.cpp VariationalMC.h random.h
	$(CXX) $(CXXFLAGS) -c main.cpp

VariationalMC.o: VariationalMC.cpp VariationalMC.h random.h
	$(CXX) $(CXXFLAGS) -c VariationalMC.cpp

random.o: random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp

clean:
	rm -f $(OBJS) $(TARGET)

run: all
	./$(TARGET) 0