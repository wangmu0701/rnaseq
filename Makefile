ADOLCPATH = $(HOME)/packages/adolc_edge
ADOLCINCL = $(ADOLCPATH)/include
ADOLCLIBS = $(ADOLCPATH)/lib64
CXX = /usr/local/bin/g++
all: adrfunc rfunc

adrfunc:	adrfunc.cpp
	$(CXX) -O3 -I$(ADOLCINCL) -I/usr/local/include/ -o adrfunc adrfunc.cpp -L$(ADOLCLIBS) -ladolc

rfunc:	rfunc.cpp
	$(CXX) -o rfunc rfunc.cpp

clean:
	rm adrfunc rfunc
