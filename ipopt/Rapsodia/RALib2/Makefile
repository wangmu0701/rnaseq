######################################################################
# This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia) #
######################################################################


# generated by: ./configure.py

ifndef CXX
  $(error "environment variable CXX is set required but not configured.")
endif

# generic rules
%.o: %.f90
	$(F90C) $(FFLAGS) $(IPATH) $(MPATH)  -c $< -o $@
        
%.o: %.F90
	$(F90C) $(FFLAGS) $(IPATH) $(MPATH)  -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(IPATH) -c $< -o $@ 

RA_GEN_NAMES= \
RA2types \
RA2add \
RA2sub \
RA2uminus \
RA2mult \
RA2div \
RA2sqrt \
RA2log \
RA2exp \
RA2pow \
RA2asin \
RA2acos \
RA2sin \
RA2cos \
RA2sinh \
RA2cosh \
RA2min \
RA2max \
RA2abs \
RA2tan \
RA2atan \
RA2gt \
RA2ge \
RA2ne \
RA2lt \
RA2le \
RA2eq \
RA2eqadd \
RA2eqsub \
RA2eqmult \
RA2eqdiv


OBJS=$(addsuffix .o, $(RA_GEN_NAMES))
default: libRapsodia.a
libRapsodia.a: $(OBJS)
	ar rcs libRapsodia.a $?


clean:
	rm -f *.o *.a


.PHONY: default clean

# dependencies

RA2div.o : RA2div.cpp RA2prec.hpp RA2divAA.ipp RA2div.hpp RA2types.hpp RA2divAP.ipp RA2divPA.ipp
RA2log.o : RA2prec.hpp RA2types.hpp RA2log.cpp RA2log.ipp RA2log.hpp
RA2abs.o : RA2types.hpp RA2prec.hpp RA2abs.ipp RA2abs.hpp RA2abs.cpp
RA2sub.o : RA2subPA.ipp RA2subAP.ipp RA2sub.cpp RA2sub.hpp RA2prec.hpp RA2types.hpp RA2subAA.ipp
RA2exp.o : RA2exp.hpp RA2exp.ipp RA2prec.hpp RA2exp.cpp RA2types.hpp
RA2acos.o : RA2acos.ipp RA2sqrt.hpp RA2prec.hpp RA2sub.hpp RA2mult.hpp RA2div.hpp RA2types.hpp RA2uminus.hpp RA2acos.hpp RA2acos.cpp
RA2uminus.o : RA2types.hpp RA2prec.hpp RA2uminus.cpp RA2uminus.ipp RA2uminus.hpp
RA2eqsub.o : RA2eqsub.cpp RA2eqsubAP.ipp RA2prec.hpp RA2eqsubAA.ipp RA2eqsub.hpp RA2types.hpp
RA2gt.o : RA2gtAP.ipp RA2gtPA.ipp RA2prec.hpp RA2gt.hpp RA2gt.cpp RA2gtAA.ipp RA2types.hpp
RA2eqadd.o : RA2eqaddAP.ipp RA2prec.hpp RA2eqaddAA.ipp RA2eqadd.hpp RA2eqadd.cpp RA2types.hpp
RA2asin.o : RA2asin.ipp RA2sub.hpp RA2asin.cpp RA2types.hpp RA2asin.hpp RA2div.hpp RA2mult.hpp RA2prec.hpp RA2sqrt.hpp
RA2sqrt.o : RA2sqrt.hpp RA2types.hpp RA2prec.hpp RA2sqrt.cpp RA2sqrt.ipp
RA2eq.o : RA2eqAA.ipp RA2eqPA.ipp RA2eqAP.ipp RA2prec.hpp RA2types.hpp RA2eq.cpp RA2eq.hpp
RA2min.o : RA2min.cpp RA2minAP.ipp RA2prec.hpp RA2minAA.ipp RA2min.hpp RA2types.hpp RA2minPA.ipp
RA2cosh.o : RA2types.hpp RA2prec.hpp RA2sinhcosh.ipp RA2cosh.hpp RA2cosh.cpp
RA2pow.o : RA2exp.hpp RA2powPA.ipp RA2powAP.ipp RA2pow.cpp RA2types.hpp RA2pow.hpp RA2prec.hpp RA2log.hpp RA2mult.hpp RA2powAA.ipp
RA2eqdiv.o : RA2eqdivAP.ipp RA2prec.hpp RA2eqdivAA.ipp RA2eqdiv.hpp RA2eqdiv.cpp RA2div.hpp RA2types.hpp
RA2sin.o : RA2sin.cpp RA2types.hpp RA2sincos.ipp RA2sin.hpp RA2prec.hpp
RA2tan.o : RA2tan.hpp RA2sin.hpp RA2prec.hpp RA2div.hpp RA2tan.cpp RA2cos.hpp RA2types.hpp RA2tan.ipp
RA2ne.o : RA2nePA.ipp RA2neAP.ipp RA2ne.cpp RA2prec.hpp RA2neAA.ipp RA2ne.hpp RA2types.hpp
RA2types.o : RA2asgnP.ipp RA2prec.hpp RA2toArray.ipp RA2get.ipp RA2asgnA.ipp RA2set.ipp RA2types.hpp RA2fromArray.ipp RA2types.cpp
RA2max.o : RA2prec.hpp RA2maxAA.ipp RA2types.hpp RA2maxAP.ipp RA2max.cpp RA2maxPA.ipp RA2max.hpp
RA2le.o : RA2leAA.ipp RA2prec.hpp RA2leAP.ipp RA2types.hpp RA2lePA.ipp RA2le.hpp RA2le.cpp
RA2add.o : RA2prec.hpp RA2addAP.ipp RA2addPA.ipp RA2addAA.ipp RA2types.hpp RA2add.hpp RA2add.cpp
RA2atan.o : RA2atan.ipp RA2prec.hpp RA2mult.hpp RA2div.hpp RA2types.hpp RA2atan.hpp RA2atan.cpp RA2add.hpp
RA2cos.o : RA2prec.hpp RA2types.hpp RA2sincos.ipp RA2cos.hpp RA2cos.cpp
RA2sinh.o : RA2sinh.hpp RA2types.hpp RA2sinhcosh.ipp RA2sinh.cpp RA2prec.hpp
RA2ge.o : RA2geAP.ipp RA2gePA.ipp RA2ge.hpp RA2ge.cpp RA2prec.hpp RA2types.hpp RA2geAA.ipp
RA2lt.o : RA2lt.hpp RA2prec.hpp RA2ltAA.ipp RA2lt.cpp RA2types.hpp RA2ltAP.ipp RA2ltPA.ipp
RA2mult.o : RA2multPA.ipp RA2prec.hpp RA2multAP.ipp RA2types.hpp RA2mult.hpp RA2multAA.ipp RA2mult.cpp
RA2eqmult.o : RA2prec.hpp RA2eqmultAA.ipp RA2eqmultAP.ipp RA2types.hpp RA2eqmult.hpp RA2eqmult.cpp
