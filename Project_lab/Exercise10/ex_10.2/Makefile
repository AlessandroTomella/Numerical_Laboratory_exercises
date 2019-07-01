CFLAGS = -Wall -O3 --std=c++11

main: main.o funz.o
	mpicxx main.o funz.o -o main

main.o: main.cpp
	mpicxx -c main.cpp -o main.o $(CFLAGS)

funz.o: funz.cpp funz.h
	g++ -c funz.cpp -o funz.o $(CFLAGS)
