// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#ifndef _RA1types_INCLUDE_
#define _RA1types_INCLUDE_
  #include "RA1prec.hpp"
  // RA1floatS
  class RA1floatS {
public:
    float v;
    float d1_1;
    float d2_1;
    float d3_1;
    float d4_1;
    float d5_1;
    float d6_1;
    float d7_1;
    RA1floatS();
    float get(const int& direction, const int& degree);
    void set(const int& direction, const int& degree, const float& passive);
    const RA1floatS& operator =(const int& r);
    const RA1floatS& operator =(const float& r);
    RA1floatS(const int& r);
    RA1floatS(const float& r);
    static const unsigned int arrSz = 8;
    void toArray(float arr[arrSz]);
    void fromArray(float arr[arrSz]);
  };
  // RA1floatD
  class RA1floatD {
public:
    double v;
    double d1_1;
    double d2_1;
    double d3_1;
    double d4_1;
    double d5_1;
    double d6_1;
    double d7_1;
    RA1floatD();
    double get(const int& direction, const int& degree);
    void set(const int& direction, const int& degree, const double& passive);
    const RA1floatD& operator =(const RA1floatS& r);
    const RA1floatD& operator =(const int& r);
    const RA1floatD& operator =(const float& r);
    const RA1floatD& operator =(const double& r);
    RA1floatD(const RA1floatS& r);
    RA1floatD(const int& r);
    RA1floatD(const float& r);
    RA1floatD(const double& r);
    static const unsigned int arrSz = 8;
    void toArray(double arr[arrSz]);
    void fromArray(double arr[arrSz]);
  };
  float makeFPE1(const float& n, const float& d);
#endif
