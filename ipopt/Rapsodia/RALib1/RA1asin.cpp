// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA1asin.hpp"

RA1floatS asin(const RA1floatS& a) {
  RA1floatS r;
  RA1floatS h;
  RA1floatS t;
  #include "RA1asin.ipp"
  return r;
}

RA1floatD asin(const RA1floatD& a) {
  RA1floatD r;
  RA1floatD h;
  RA1floatD t;
  #include "RA1asin.ipp"
  return r;
}
