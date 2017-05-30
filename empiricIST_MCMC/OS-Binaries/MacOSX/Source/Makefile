CC=g++
VERSION = $(shell cat ./version)
CFLAGS=-c -O3 -m64 -std=c++11 -DVERSION=\"${VERSION}\" -I/usr/local/include/ -L/usr/local/lib/
LDFLAGS=-lgsl -lgslcblas 
SOURCES=empiricIST_MCMC.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=empiricIST_MCMC.out

all: clean_all $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *.o
	
.PHONY: clean_all
clean_all:
	rm -f *.o $(EXECUTABLE)