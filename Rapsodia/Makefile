##########################################################
# This file is part of Rapsodia released under the LGPL. #
# The full COPYRIGHT notice can be found in the top      #
# level directory of the Rapsodia distribution           #
##########################################################
include ../../../Downloads/Rapsodia/MakeDefs.mk
include ../../../Downloads/Rapsodia/MakeRules.mk

default: driver
	./$^

GEN_DIR=RALib

RA_EXTRAS=../../../Downloads/Rapsodia/hotCpp
IPATH+=-I$(GEN_DIR) -I$(RA_EXTRAS)

OBJS= \
$(addprefix $(RA_EXTRAS)/, $(addsuffix .o, $(HOTCPPNAMES))) \
adrfunc.o \
driver.o

driver:  $(OBJS) $(GEN_DIR)/libRapsodia.a
	$(CXX) $(CXXFLAGS) $(IPATH) -o $@ $^

$(GEN_DIR)/libRapsodia.a: FORCE
	../../../Downloads/Rapsodia/Generator/generate.py -d 28 -o 2 -c $(GEN_DIR)
	cd $(GEN_DIR) && $(MAKE)

FORCE:

clean: 
	rm -rf $(GEN_DIR) *.o driver

.PHONY: default clean

adrfunc.cpp:$(GEN_DIR)/libRapsodia.a
driver.cpp:$(GEN_DIR)/libRapsodia.a
