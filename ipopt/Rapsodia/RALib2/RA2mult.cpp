// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA2mult.hpp"

RA2floatS operator *(const RA2floatS& a, const RA2floatS& b) {
  RA2floatS r;
  #include "RA2multAA.ipp"
  return r;
}

RA2floatD operator *(const RA2floatS& a, const RA2floatD& b) {
  RA2floatD r;
  #include "RA2multAA.ipp"
  return r;
}

RA2floatD operator *(const RA2floatD& a, const RA2floatS& b) {
  RA2floatD r;
  #include "RA2multAA.ipp"
  return r;
}

RA2floatD operator *(const RA2floatD& a, const RA2floatD& b) {
  RA2floatD r;
  #include "RA2multAA.ipp"
  return r;
}

RA2floatS operator *(const RA2floatS& a, const int& b) {
  RA2floatS r;
  #include "RA2multAP.ipp"
  return r;
}

RA2floatS operator *(const RA2floatS& a, const float& b) {
  RA2floatS r;
  #include "RA2multAP.ipp"
  return r;
}

RA2floatD operator *(const RA2floatS& a, const double& b) {
  RA2floatD r;
  #include "RA2multAP.ipp"
  return r;
}

RA2floatD operator *(const RA2floatD& a, const int& b) {
  RA2floatD r;
  #include "RA2multAP.ipp"
  return r;
}

RA2floatD operator *(const RA2floatD& a, const float& b) {
  RA2floatD r;
  #include "RA2multAP.ipp"
  return r;
}

RA2floatD operator *(const RA2floatD& a, const double& b) {
  RA2floatD r;
  #include "RA2multAP.ipp"
  return r;
}

RA2floatS operator *(const int& a, const RA2floatS& b) {
  RA2floatS r;
  #include "RA2multPA.ipp"
  return r;
}

RA2floatD operator *(const int& a, const RA2floatD& b) {
  RA2floatD r;
  #include "RA2multPA.ipp"
  return r;
}

RA2floatS operator *(const float& a, const RA2floatS& b) {
  RA2floatS r;
  #include "RA2multPA.ipp"
  return r;
}

RA2floatD operator *(const float& a, const RA2floatD& b) {
  RA2floatD r;
  #include "RA2multPA.ipp"
  return r;
}

RA2floatD operator *(const double& a, const RA2floatS& b) {
  RA2floatD r;
  #include "RA2multPA.ipp"
  return r;
}

RA2floatD operator *(const double& a, const RA2floatD& b) {
  RA2floatD r;
  #include "RA2multPA.ipp"
  return r;
}
