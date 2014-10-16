// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA2atan.hpp"

RA2floatS atan(const RA2floatS& a) {
  RA2floatS r;
  RA2floatS y;
  RA2floatS s;
  float one;
  #include "RA2atan.ipp"
  return r;
}

RA2floatD atan(const RA2floatD& a) {
  RA2floatD r;
  RA2floatD y;
  RA2floatD s;
  double one;
  #include "RA2atan.ipp"
  return r;
}
