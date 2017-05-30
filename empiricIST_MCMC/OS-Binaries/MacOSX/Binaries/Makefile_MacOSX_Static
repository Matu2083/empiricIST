CC=g++-mp-4.9
VERSION = $(shell cat ./version)
CFLAGS=-O3 -m64 -DVERSION=\"${VERSION}\" -Wall libgsl.a libgslcblas.a libstdc++.a libgomp.a libgcc.a -I . 
SOURCES=empiricIST_MCMC.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=empiricIST_MCMC.out

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(SOURCES)
	$(CC) $(SOURCES) $(CFLAGS) -o $(EXECUTABLE)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *.o
	
.PHONY: clean_all
clean_all:
	rm -f *.o $(EXECUTABLE)