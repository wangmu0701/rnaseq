# CHANGEME: This should be the name of your executable
EXE = cpp_example

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  cpp_example.o \
	rapsodia_NLP.o \
        RALib1/libRapsodia.a \
        RALib2/libRapsodia.a

RAPSODIAHOME = $(HOME)/packages/Rapsodia/
RAPSODIAINCL = $(RAPSODIAHOME)/hotCpp/
RAPSODIALIBS = $(RAPSODIAHOME)/hotCpp/MultiIndex.o $(RAPSODIAHOME)/hotCpp/HigherOrderTensor.o ./RALib1/libRapsodia.a ./RALib2/libRapsodia.a

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I$(RAPSODIAINCL)  

ADDLIBS=$(RAPSODIALIBS)

include ./../MakeEnv.mk

RALib1/libRapsodia.a:
	make -C ./RALib1
RALib2/libRapsodia.a:
	make -C ./RALib2

clean:
	rm -rf $(EXE) $(OBJS); make clean -C ./RALib1; make clean -C ./RALib2
