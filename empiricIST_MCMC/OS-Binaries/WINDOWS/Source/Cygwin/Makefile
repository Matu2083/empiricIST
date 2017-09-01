CC=g++
VERSION=$(shell cat ./version)
CFLAGS=-c -O3 -DVERSION=\"${VERSION}\"
LDFLAGS=-lgsl -lgslcblas
SOURCES=empiricIST_MCMC.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=empiricIST_MCMC.exe

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
