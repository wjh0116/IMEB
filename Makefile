CXX = g++
CC = gcc
CXXFLAGS = -O3 -std=c++17 -DSFMT_MEXP=19937
CFLAGS = -O3 -DSFMT_MEXP=19937
SFMT_OBJ = sfmt/SFMT.o
TARGET = IMIN_EB

all: $(TARGET)

SABRE_HDRS = ForwardSample.h CriticalEdgeSet.h AncestralEdge.h Greedy.h SABRE.h

$(SFMT_OBJ): sfmt/SFMT.c sfmt/SFMT.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(TARGET): main.cpp head.h graph.h infgraph.h iheap.h DESCE.h EdgeBlock.h CriticalEdge.h $(SABRE_HDRS) $(SFMT_OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.cpp $(SFMT_OBJ)

clean:
	rm -f $(TARGET) $(SFMT_OBJ)

.PHONY: all clean
