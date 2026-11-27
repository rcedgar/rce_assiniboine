CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -pthread

.PHONY: all clean

all: vleness

vleness: vleness.cpp vleness.h ../rce_naview/asserta.hpp
	$(CXX) $(CXXFLAGS) -o vleness vleness.cpp

clean:
	rm -f vleness
