// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA2sin.hpp"

RA2floatS sin(const RA2floatS& a) {
  RA2floatS s;
  RA2floatS c;
  RA2floatS t;
  #include "RA2sincos.ipp"
  return s;
}

RA2floatD sin(const RA2floatD& a) {
  RA2floatD s;
  RA2floatD c;
  RA2floatD t;
  #include "RA2sincos.ipp"
  return s;
}
