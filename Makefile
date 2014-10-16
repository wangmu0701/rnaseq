ADOLC_PATH = /home/zsdfe/adolc_edge
ADOLC_LIB = lib
all: adrfunc rfunc

adrfunc:	adrfunc.cpp
	g++ -O3 -I$(ADOLC_PATH)/include/ -I/usr/local/include/ -o adrfunc adrfunc.cpp -L$(ADOLC_PATH)/$(ADOLC_LIB) -ladolc

rfunc:	rfunc.cpp
	g++ -o rfunc rfunc.cpp

clean:
	rm adrfunc rfunc
