#Basic makefile for C++ code
CC=g++ #compiler
CFLAGS=-g -c -Wall -I/usr/local/include/ -std=c++0x -O3 #compiler flags
LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lm #linker flags
SOURCES=ddrd_fluctuations.cpp #source files
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bin_ddrd_fluctuations #program name

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJECTS): $(SOURCES)
	$(CC) $(CFLAGS) $<

clean:
	rm $(OBJECTS)
