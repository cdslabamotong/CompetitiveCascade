_MAIN_=start.cpp
PARALLELFLAGS=-lboost_thread -lboost_system -lpthread
CXXFLAGS=-std=c++11 -O3

start.out: $(_MAIN_)
	$(CXX) -o $@ $^ $(PARALLELFLAGS) $(CXXFLAGS)

all: start.out

clean:
	rm -f start.out