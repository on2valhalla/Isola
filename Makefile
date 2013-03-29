
CC  = gcc
CXX = g++


INCLUDES =


CFLAGS   = -g -Wall $(INCLUDES)
CXXFLAGS = -g -Wall $(INCLUDES)

LDFLAGS = -g $(LDLIBS)
LDLIBS = 

isola: isola.cpp
	g++ -std=c++0x -o isola isola.cpp

.PHONY: clean
clean:
	rm -f *.o *~ a.out core isola

.PHONY: all
all: clean isola
