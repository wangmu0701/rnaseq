// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
#ifndef _RA2types_INCLUDE_
#define _RA2types_INCLUDE_
  #include "RA2prec.hpp"
  // RA2floatS
  class RA2floatS {
public:
    float v;
    float d1_1;
    float d1_2;
    float d2_1;
    float d2_2;
    float d3_1;
    float d3_2;
    float d4_1;
    float d4_2;
    float d5_1;
    float d5_2;
    float d6_1;
    float d6_2;
    float d7_1;
    float d7_2;
    float d8_1;
    float d8_2;
    float d9_1;
    float d9_2;
    float d10_1;
    float d10_2;
    float d11_1;
    float d11_2;
    float d12_1;
    float d12_2;
    float d13_1;
    float d13_2;
    float d14_1;
    float d14_2;
    float d15_1;
    float d15_2;
    float d16_1;
    float d16_2;
    float d17_1;
    float d17_2;
    float d18_1;
    float d18_2;
    float d19_1;
    float d19_2;
    float d20_1;
    float d20_2;
    float d21_1;
    float d21_2;
    float d22_1;
    float d22_2;
    float d23_1;
    float d23_2;
    float d24_1;
    float d24_2;
    float d25_1;
    float d25_2;
    float d26_1;
    float d26_2;
    float d27_1;
    float d27_2;
    float d28_1;
    float d28_2;
    RA2floatS();
    float get(const int& direction, const int& degree);
    void set(const int& direction, const int& degree, const float& passive);
    const RA2floatS& operator =(const int& r);
    const RA2floatS& operator =(const float& r);
    RA2floatS(const int& r);
    RA2floatS(const float& r);
    static const unsigned int arrSz = 57;
    void toArray(float arr[arrSz]);
    void fromArray(float arr[arrSz]);
  };
  // RA2floatD
  class RA2floatD {
public:
    double v;
    double d1_1;
    double d1_2;
    double d2_1;
    double d2_2;
    double d3_1;
    double d3_2;
    double d4_1;
    double d4_2;
    double d5_1;
    double d5_2;
    double d6_1;
    double d6_2;
    double d7_1;
    double d7_2;
    double d8_1;
    double d8_2;
    double d9_1;
    double d9_2;
    double d10_1;
    double d10_2;
    double d11_1;
    double d11_2;
    double d12_1;
    double d12_2;
    double d13_1;
    double d13_2;
    double d14_1;
    double d14_2;
    double d15_1;
    double d15_2;
    double d16_1;
    double d16_2;
    double d17_1;
    double d17_2;
    double d18_1;
    double d18_2;
    double d19_1;
    double d19_2;
    double d20_1;
    double d20_2;
    double d21_1;
    double d21_2;
    double d22_1;
    double d22_2;
    double d23_1;
    double d23_2;
    double d24_1;
    double d24_2;
    double d25_1;
    double d25_2;
    double d26_1;
    double d26_2;
    double d27_1;
    double d27_2;
    double d28_1;
    double d28_2;
    RA2floatD();
    double get(const int& direction, const int& degree);
    void set(const int& direction, const int& degree, const double& passive);
    const RA2floatD& operator =(const RA2floatS& r);
    const RA2floatD& operator =(const int& r);
    const RA2floatD& operator =(const float& r);
    const RA2floatD& operator =(const double& r);
    RA2floatD(const RA2floatS& r);
    RA2floatD(const int& r);
    RA2floatD(const float& r);
    RA2floatD(const double& r);
    static const unsigned int arrSz = 57;
    void toArray(double arr[arrSz]);
    void fromArray(double arr[arrSz]);
  };
  float makeFPE(const float& n, const float& d);
#endif
