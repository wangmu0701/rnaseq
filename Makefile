all: adrfunc rfunc

adrfunc:	adrfunc.cpp
	/usr/local/bin/g++ -O3 -I/Users/muwang/packages/adolc_edge/include/ -I/usr/local/include/ -o adrfunc adrfunc.cpp -L/Users/muwang/packages/adolc_edge/lib64 -ladolc

rfunc:	rfunc.cpp
	/usr/local/bin/g++ -o rfunc rfunc.cpp

clean:
	rm adrfunc rfunc