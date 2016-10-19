CFLAGS=-ggdb3 -std=c++11

all: casim.o main.o params.h
	g++ $(CFLAGS) -o evol main.o casim.o

main.o: main.cpp
	g++ $(CFLAGS) -c -o main.o main.cpp 


casim.o: casim.cpp params.h
	g++ $(CFLAGS) -c -o casim.o casim.cpp

casim: casim_main.cpp casim.o
	g++ $(CFLAGS) -o casim  casim.o casim_main.cpp

test: casim.o catest.cpp
	g++  $(CFLAGS) -o test catest.cpp casim.o

run:
	./evol

clean:
	rm -f evol
	rm -f casim
	rm -f *.o
