CC=g++
VERSION=$(shell cat ./version)
CFLAGS=-c -O3 -DVERSION=\"${VERSION}\" -I/local/include/ -I/local/bin/ -I/local/lib/
LDFLAGS=-L/local/include/ -L/local/bin/ -L/local/lib/ -lgsl -lgslcblas
SOURCES=empiricIST_MCMC.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=empiricIST_MCMC.exe

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
