// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA1mult.hpp"

RA1floatS operator *(const RA1floatS& a, const RA1floatS& b) {
  RA1floatS r;
  #include "RA1multAA.ipp"
  return r;
}

RA1floatD operator *(const RA1floatS& a, const RA1floatD& b) {
  RA1floatD r;
  #include "RA1multAA.ipp"
  return r;
}

RA1floatD operator *(const RA1floatD& a, const RA1floatS& b) {
  RA1floatD r;
  #include "RA1multAA.ipp"
  return r;
}

RA1floatD operator *(const RA1floatD& a, const RA1floatD& b) {
  RA1floatD r;
  #include "RA1multAA.ipp"
  return r;
}

RA1floatS operator *(const RA1floatS& a, const int& b) {
  RA1floatS r;
  #include "RA1multAP.ipp"
  return r;
}

RA1floatS operator *(const RA1floatS& a, const float& b) {
  RA1floatS r;
  #include "RA1multAP.ipp"
  return r;
}

RA1floatD operator *(const RA1floatS& a, const double& b) {
  RA1floatD r;
  #include "RA1multAP.ipp"
  return r;
}

RA1floatD operator *(const RA1floatD& a, const int& b) {
  RA1floatD r;
  #include "RA1multAP.ipp"
  return r;
}

RA1floatD operator *(const RA1floatD& a, const float& b) {
  RA1floatD r;
  #include "RA1multAP.ipp"
  return r;
}

RA1floatD operator *(const RA1floatD& a, const double& b) {
  RA1floatD r;
  #include "RA1multAP.ipp"
  return r;
}

RA1floatS operator *(const int& a, const RA1floatS& b) {
  RA1floatS r;
  #include "RA1multPA.ipp"
  return r;
}

RA1floatD operator *(const int& a, const RA1floatD& b) {
  RA1floatD r;
  #include "RA1multPA.ipp"
  return r;
}

RA1floatS operator *(const float& a, const RA1floatS& b) {
  RA1floatS r;
  #include "RA1multPA.ipp"
  return r;
}

RA1floatD operator *(const float& a, const RA1floatD& b) {
  RA1floatD r;
  #include "RA1multPA.ipp"
  return r;
}

RA1floatD operator *(const double& a, const RA1floatS& b) {
  RA1floatD r;
  #include "RA1multPA.ipp"
  return r;
}

RA1floatD operator *(const double& a, const RA1floatD& b) {
  RA1floatD r;
  #include "RA1multPA.ipp"
  return r;
}

