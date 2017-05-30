CC=g++
CFLAGS=-c -O3
LDFLAGS=-lgsl -lgslcblas
SOURCES=empiricIST_MCMC-Cluster.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=empiricIST_MCMC-Cluster.out

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


