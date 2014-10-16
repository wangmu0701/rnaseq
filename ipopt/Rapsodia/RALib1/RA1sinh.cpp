// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA1sinh.hpp"

RA1floatS sinh(const RA1floatS& a) {
  RA1floatS s;
  RA1floatS c;
  RA1floatS t;
  #include "RA1sinhcosh.ipp"
  return s;
}

RA1floatD sinh(const RA1floatD& a) {
  RA1floatD s;
  RA1floatD c;
  RA1floatD t;
  #include "RA1sinhcosh.ipp"
  return s;
}

