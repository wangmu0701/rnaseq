# ----------------------------------------------------------------------------
# ADOL-C -- Automatic Differentiation by Overloading in C++
# File:     Makefile
# Revision: $$
# Contents: Makefile for class myADOLC_NPL to interface ADOL-C with Ipopt
# 
# Copyright (c) Andrea Walther
#   
# This file is part of ADOL-C. This software is provided as open source.
# Any use, reproduction, or distribution of the software constitutes 
# recipient's acceptance of the terms of the accompanying license file.
# 
# This code is based on the corresponding makefile contained in the
# Ipopt package with the author Andreas Waechter   
#---------------------------------------------------------------------------

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGE entries below.     #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = cpp_example

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  cpp_example.o \
	ADOL-C_NLP.o

# Definitions for ADOL-C
# CHANGE THIS TO YOUR ADOL-C BASE DIRECTORY
ADPATH = $(HOME)/packages/adolc_edge/include
ADLIBDIR = $(HOME)/packages/adolc_edge/lib64

# CHANGEME: Additional libraries
ADDLIBS =  -L$(ADLIBDIR) -ladolc -L$(HOME)/packages/ColPack/lib/ -lColPack

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I$(ADPATH)  

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C++ Compiler command
CXX = /usr/local/bin/g++

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD

# additional C++ Compiler options for linking
CXXLINKFLAGS = 

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
#INCL = `PKG_CONFIG_PATH=/Users/muwang/packages/ipopt/lib64/pkgconfig:/Users/muwang/packages/ipopt/lib/pkgconfig:/Users/muwang/packages/ipopt/share/pkgconfig:  --cflags ipopt` $(ADDINCFLAGS)
INCL = -I/Users/muwang/packages/ipopt/include/coin  -I/Users/muwang/packages/ipopt/include/coin/ThirdParty   $(ADDINCFLAGS)

# Linker flags
LIBS = -L/Users/muwang/packages/ipopt/lib -lipopt  -L/Users/muwang/packages/ipopt/lib -lcoinhsl -L/usr/local/lib/gcc/x86_64-apple-darwin13.0.2/4.8.2 -L/usr/local/lib/gcc/x86_64-apple-darwin13.0.2/4.8.2/../../.. -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/lib -lgfortran -lSystem -lquadmath -lm  -llapack -lblas -lm  -ldl

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `$(CYGPATH_W) '$<'`



