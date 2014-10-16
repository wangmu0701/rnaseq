// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "RA1types.hpp"

const unsigned int RA1floatS::arrSz;
RA1floatS::RA1floatS() {
  RA1floatS& l = *this;
  const float r = 0.0F;
  #include "RA1asgnP.ipp"
}

void RA1floatS::set(const int& direction, const int& degree, const float& passive) {
  #include "RA1set.ipp"
}

float RA1floatS::get(const int& direction, const int& degree) {
  float passive;
  #include "RA1get.ipp"
  return passive;
}

const RA1floatS& RA1floatS::operator =(const int& r) {
  const RA1floatS& ret = *this;
  RA1floatS& l = *this;
  #include "RA1asgnP.ipp"
  return ret;
}

const RA1floatS& RA1floatS::operator =(const float& r) {
  const RA1floatS& ret = *this;
  RA1floatS& l = *this;
  #include "RA1asgnP.ipp"
  return ret;
}

RA1floatS::RA1floatS(const int& r) {
  const RA1floatS& ret = *this;
  RA1floatS& l = *this;
  #include "RA1asgnP.ipp"
}

RA1floatS::RA1floatS(const float& r) {
  const RA1floatS& ret = *this;
  RA1floatS& l = *this;
  #include "RA1asgnP.ipp"
}

void RA1floatS::toArray(float arr[RA1floatS::arrSz]) {
  RA1floatS& l = *this;
  arr[0] = l.v;
  #include "RA1toArray.ipp"
}

void RA1floatS::fromArray(float arr[RA1floatS::arrSz]) {
  RA1floatS& l = *this;
  l.v = arr[0];
  #include "RA1fromArray.ipp"
}

const unsigned int RA1floatD::arrSz;
RA1floatD::RA1floatD() {
  RA1floatD& l = *this;
  const double r = 0.0;
  #include "RA1asgnP.ipp"
}

void RA1floatD::set(const int& direction, const int& degree, const double& passive) {
  #include "RA1set.ipp"
}

double RA1floatD::get(const int& direction, const int& degree) {
  double passive;
  #include "RA1get.ipp"
  return passive;
}

const RA1floatD& RA1floatD::operator =(const RA1floatS& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnA.ipp"
  return ret;
}

const RA1floatD& RA1floatD::operator =(const int& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnP.ipp"
  return ret;
}

const RA1floatD& RA1floatD::operator =(const float& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnP.ipp"
  return ret;
}

const RA1floatD& RA1floatD::operator =(const double& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnP.ipp"
  return ret;
}

RA1floatD::RA1floatD(const RA1floatS& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnA.ipp"
}

RA1floatD::RA1floatD(const int& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnP.ipp"
}

RA1floatD::RA1floatD(const float& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnP.ipp"
}

RA1floatD::RA1floatD(const double& r) {
  const RA1floatD& ret = *this;
  RA1floatD& l = *this;
  #include "RA1asgnP.ipp"
}

void RA1floatD::toArray(double arr[RA1floatD::arrSz]) {
  RA1floatD& l = *this;
  arr[0] = l.v;
  #include "RA1toArray.ipp"
}

void RA1floatD::fromArray(double arr[RA1floatD::arrSz]) {
  RA1floatD& l = *this;
  l.v = arr[0];
  #include "RA1fromArray.ipp"
}

float makeFPE1(const float& n, const float& d) {
  float r;
  r = n / d;
  return r;
}
