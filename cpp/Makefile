CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

all: rk4_nbody euler_nbody

rk4_nbody: rk4_nbody.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

euler_nbody: euler_nbody.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f rk4_nbody euler_nbody

.PHONY: all clean
