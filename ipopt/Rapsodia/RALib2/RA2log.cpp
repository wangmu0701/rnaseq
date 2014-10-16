// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA2log.hpp"

RA2floatS log(const RA2floatS& a) {
  RA2floatS r;
  RA2floatS s;
  RA2floatS t;
  float recip;
  #include "RA2log.ipp"
  return r;
}

RA2floatD log(const RA2floatD& a) {
  RA2floatD r;
  RA2floatD s;
  RA2floatD t;
  double recip;
  #include "RA2log.ipp"
  return r;
}
