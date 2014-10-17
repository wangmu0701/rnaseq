# CHANGEME: This should be the name of your executable
EXE = cpp_example

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  cpp_example.o \
	ADOL-C_NLP.o

# CHANGE THIS TO YOUR ADOL-C BASE DIRECTORY
ADPATH = $(HOME)/packages/adolc_edge
ADINCL = $(ADPATH)/include

# CHANGEME: Additional libraries
ADDLIBS =  -L$(ADPATH)/lib64 -ladolc -L$(HOME)/packages/ColPack/lib/ -lColPack

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I$(ADINCL)  

include ./MakeEnv.mk


