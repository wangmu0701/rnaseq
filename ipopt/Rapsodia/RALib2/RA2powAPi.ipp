// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
r.v = pow(double(a.v), double(b));
if (a.v == 0.0)
  {
    if (b <= 0)
      {
        r.d1_1 = makeFPE(0.0, 0.0);
        r.d1_2 = makeFPE(0.0, 0.0);
        r.d2_1 = makeFPE(0.0, 0.0);
        r.d2_2 = makeFPE(0.0, 0.0);
        r.d3_1 = makeFPE(0.0, 0.0);
        r.d3_2 = makeFPE(0.0, 0.0);
        r.d4_1 = makeFPE(0.0, 0.0);
        r.d4_2 = makeFPE(0.0, 0.0);
        r.d5_1 = makeFPE(0.0, 0.0);
        r.d5_2 = makeFPE(0.0, 0.0);
        r.d6_1 = makeFPE(0.0, 0.0);
        r.d6_2 = makeFPE(0.0, 0.0);
        r.d7_1 = makeFPE(0.0, 0.0);
        r.d7_2 = makeFPE(0.0, 0.0);
        r.d8_1 = makeFPE(0.0, 0.0);
        r.d8_2 = makeFPE(0.0, 0.0);
        r.d9_1 = makeFPE(0.0, 0.0);
        r.d9_2 = makeFPE(0.0, 0.0);
        r.d10_1 = makeFPE(0.0, 0.0);
        r.d10_2 = makeFPE(0.0, 0.0);
        r.d11_1 = makeFPE(0.0, 0.0);
        r.d11_2 = makeFPE(0.0, 0.0);
        r.d12_1 = makeFPE(0.0, 0.0);
        r.d12_2 = makeFPE(0.0, 0.0);
        r.d13_1 = makeFPE(0.0, 0.0);
        r.d13_2 = makeFPE(0.0, 0.0);
        r.d14_1 = makeFPE(0.0, 0.0);
        r.d14_2 = makeFPE(0.0, 0.0);
        r.d15_1 = makeFPE(0.0, 0.0);
        r.d15_2 = makeFPE(0.0, 0.0);
        r.d16_1 = makeFPE(0.0, 0.0);
        r.d16_2 = makeFPE(0.0, 0.0);
        r.d17_1 = makeFPE(0.0, 0.0);
        r.d17_2 = makeFPE(0.0, 0.0);
        r.d18_1 = makeFPE(0.0, 0.0);
        r.d18_2 = makeFPE(0.0, 0.0);
        r.d19_1 = makeFPE(0.0, 0.0);
        r.d19_2 = makeFPE(0.0, 0.0);
        r.d20_1 = makeFPE(0.0, 0.0);
        r.d20_2 = makeFPE(0.0, 0.0);
        r.d21_1 = makeFPE(0.0, 0.0);
        r.d21_2 = makeFPE(0.0, 0.0);
        r.d22_1 = makeFPE(0.0, 0.0);
        r.d22_2 = makeFPE(0.0, 0.0);
        r.d23_1 = makeFPE(0.0, 0.0);
        r.d23_2 = makeFPE(0.0, 0.0);
        r.d24_1 = makeFPE(0.0, 0.0);
        r.d24_2 = makeFPE(0.0, 0.0);
        r.d25_1 = makeFPE(0.0, 0.0);
        r.d25_2 = makeFPE(0.0, 0.0);
        r.d26_1 = makeFPE(0.0, 0.0);
        r.d26_2 = makeFPE(0.0, 0.0);
        r.d27_1 = makeFPE(0.0, 0.0);
        r.d27_2 = makeFPE(0.0, 0.0);
        r.d28_1 = makeFPE(0.0, 0.0);
        r.d28_2 = makeFPE(0.0, 0.0);
      }
    else
      {
        {
          if (b == 1)
            {
              r.d1_1 = a.d1_1;
              r.d1_2 = a.d1_2;
              r.d2_1 = a.d2_1;
              r.d2_2 = a.d2_2;
              r.d3_1 = a.d3_1;
              r.d3_2 = a.d3_2;
              r.d4_1 = a.d4_1;
              r.d4_2 = a.d4_2;
              r.d5_1 = a.d5_1;
              r.d5_2 = a.d5_2;
              r.d6_1 = a.d6_1;
              r.d6_2 = a.d6_2;
              r.d7_1 = a.d7_1;
              r.d7_2 = a.d7_2;
              r.d8_1 = a.d8_1;
              r.d8_2 = a.d8_2;
              r.d9_1 = a.d9_1;
              r.d9_2 = a.d9_2;
              r.d10_1 = a.d10_1;
              r.d10_2 = a.d10_2;
              r.d11_1 = a.d11_1;
              r.d11_2 = a.d11_2;
              r.d12_1 = a.d12_1;
              r.d12_2 = a.d12_2;
              r.d13_1 = a.d13_1;
              r.d13_2 = a.d13_2;
              r.d14_1 = a.d14_1;
              r.d14_2 = a.d14_2;
              r.d15_1 = a.d15_1;
              r.d15_2 = a.d15_2;
              r.d16_1 = a.d16_1;
              r.d16_2 = a.d16_2;
              r.d17_1 = a.d17_1;
              r.d17_2 = a.d17_2;
              r.d18_1 = a.d18_1;
              r.d18_2 = a.d18_2;
              r.d19_1 = a.d19_1;
              r.d19_2 = a.d19_2;
              r.d20_1 = a.d20_1;
              r.d20_2 = a.d20_2;
              r.d21_1 = a.d21_1;
              r.d21_2 = a.d21_2;
              r.d22_1 = a.d22_1;
              r.d22_2 = a.d22_2;
              r.d23_1 = a.d23_1;
              r.d23_2 = a.d23_2;
              r.d24_1 = a.d24_1;
              r.d24_2 = a.d24_2;
              r.d25_1 = a.d25_1;
              r.d25_2 = a.d25_2;
              r.d26_1 = a.d26_1;
              r.d26_2 = a.d26_2;
              r.d27_1 = a.d27_1;
              r.d27_2 = a.d27_2;
              r.d28_1 = a.d28_1;
              r.d28_2 = a.d28_2;
            }
          else
            {
              r.d1_1 = 0.0;
              r.d1_2 = a.d1_1 * a.d1_1;
              r.d2_1 = 0.0;
              r.d2_2 = a.d2_1 * a.d2_1;
              r.d3_1 = 0.0;
              r.d3_2 = a.d3_1 * a.d3_1;
              r.d4_1 = 0.0;
              r.d4_2 = a.d4_1 * a.d4_1;
              r.d5_1 = 0.0;
              r.d5_2 = a.d5_1 * a.d5_1;
              r.d6_1 = 0.0;
              r.d6_2 = a.d6_1 * a.d6_1;
              r.d7_1 = 0.0;
              r.d7_2 = a.d7_1 * a.d7_1;
              r.d8_1 = 0.0;
              r.d8_2 = a.d8_1 * a.d8_1;
              r.d9_1 = 0.0;
              r.d9_2 = a.d9_1 * a.d9_1;
              r.d10_1 = 0.0;
              r.d10_2 = a.d10_1 * a.d10_1;
              r.d11_1 = 0.0;
              r.d11_2 = a.d11_1 * a.d11_1;
              r.d12_1 = 0.0;
              r.d12_2 = a.d12_1 * a.d12_1;
              r.d13_1 = 0.0;
              r.d13_2 = a.d13_1 * a.d13_1;
              r.d14_1 = 0.0;
              r.d14_2 = a.d14_1 * a.d14_1;
              r.d15_1 = 0.0;
              r.d15_2 = a.d15_1 * a.d15_1;
              r.d16_1 = 0.0;
              r.d16_2 = a.d16_1 * a.d16_1;
              r.d17_1 = 0.0;
              r.d17_2 = a.d17_1 * a.d17_1;
              r.d18_1 = 0.0;
              r.d18_2 = a.d18_1 * a.d18_1;
              r.d19_1 = 0.0;
              r.d19_2 = a.d19_1 * a.d19_1;
              r.d20_1 = 0.0;
              r.d20_2 = a.d20_1 * a.d20_1;
              r.d21_1 = 0.0;
              r.d21_2 = a.d21_1 * a.d21_1;
              r.d22_1 = 0.0;
              r.d22_2 = a.d22_1 * a.d22_1;
              r.d23_1 = 0.0;
              r.d23_2 = a.d23_1 * a.d23_1;
              r.d24_1 = 0.0;
              r.d24_2 = a.d24_1 * a.d24_1;
              r.d25_1 = 0.0;
              r.d25_2 = a.d25_1 * a.d25_1;
              r.d26_1 = 0.0;
              r.d26_2 = a.d26_1 * a.d26_1;
              r.d27_1 = 0.0;
              r.d27_2 = a.d27_1 * a.d27_1;
              r.d28_1 = 0.0;
              r.d28_2 = a.d28_1 * a.d28_1;
              for(j=3;j<=int(b);j+=1)
                {
                  r.d1_2 = r.v * a.d1_1 + r.d1_1 * a.v;
                  r.d2_2 = r.v * a.d2_1 + r.d2_1 * a.v;
                  r.d3_2 = r.v * a.d3_1 + r.d3_1 * a.v;
                  r.d4_2 = r.v * a.d4_1 + r.d4_1 * a.v;
                  r.d5_2 = r.v * a.d5_1 + r.d5_1 * a.v;
                  r.d6_2 = r.v * a.d6_1 + r.d6_1 * a.v;
                  r.d7_2 = r.v * a.d7_1 + r.d7_1 * a.v;
                  r.d8_2 = r.v * a.d8_1 + r.d8_1 * a.v;
                  r.d9_2 = r.v * a.d9_1 + r.d9_1 * a.v;
                  r.d10_2 = r.v * a.d10_1 + r.d10_1 * a.v;
                  r.d11_2 = r.v * a.d11_1 + r.d11_1 * a.v;
                  r.d12_2 = r.v * a.d12_1 + r.d12_1 * a.v;
                  r.d13_2 = r.v * a.d13_1 + r.d13_1 * a.v;
                  r.d14_2 = r.v * a.d14_1 + r.d14_1 * a.v;
                  r.d15_2 = r.v * a.d15_1 + r.d15_1 * a.v;
                  r.d16_2 = r.v * a.d16_1 + r.d16_1 * a.v;
                  r.d17_2 = r.v * a.d17_1 + r.d17_1 * a.v;
                  r.d18_2 = r.v * a.d18_1 + r.d18_1 * a.v;
                  r.d19_2 = r.v * a.d19_1 + r.d19_1 * a.v;
                  r.d20_2 = r.v * a.d20_1 + r.d20_1 * a.v;
                  r.d21_2 = r.v * a.d21_1 + r.d21_1 * a.v;
                  r.d22_2 = r.v * a.d22_1 + r.d22_1 * a.v;
                  r.d23_2 = r.v * a.d23_1 + r.d23_1 * a.v;
                  r.d24_2 = r.v * a.d24_1 + r.d24_1 * a.v;
                  r.d25_2 = r.v * a.d25_1 + r.d25_1 * a.v;
                  r.d26_2 = r.v * a.d26_1 + r.d26_1 * a.v;
                  r.d27_2 = r.v * a.d27_1 + r.d27_1 * a.v;
                  r.d28_2 = r.v * a.d28_1 + r.d28_1 * a.v;
                  r.d1_1 = r.v * a.v;
                  r.d2_1 = r.v * a.v;
                  r.d3_1 = r.v * a.v;
                  r.d4_1 = r.v * a.v;
                  r.d5_1 = r.v * a.v;
                  r.d6_1 = r.v * a.v;
                  r.d7_1 = r.v * a.v;
                  r.d8_1 = r.v * a.v;
                  r.d9_1 = r.v * a.v;
                  r.d10_1 = r.v * a.v;
                  r.d11_1 = r.v * a.v;
                  r.d12_1 = r.v * a.v;
                  r.d13_1 = r.v * a.v;
                  r.d14_1 = r.v * a.v;
                  r.d15_1 = r.v * a.v;
                  r.d16_1 = r.v * a.v;
                  r.d17_1 = r.v * a.v;
                  r.d18_1 = r.v * a.v;
                  r.d19_1 = r.v * a.v;
                  r.d20_1 = r.v * a.v;
                  r.d21_1 = r.v * a.v;
                  r.d22_1 = r.v * a.v;
                  r.d23_1 = r.v * a.v;
                  r.d24_1 = r.v * a.v;
                  r.d25_1 = r.v * a.v;
                  r.d26_1 = r.v * a.v;
                  r.d27_1 = r.v * a.v;
                  r.d28_1 = r.v * a.v;
                }
            }
        }
      }
  }
else
  {
    recip = 1.0 / a.v;
    s.d1_1 = 1 * a.d1_1;
    t.d1_1 = recip * (b * (r.v * s.d1_1) - (0.0));
    r.d1_1 = t.d1_1 / 1;
    s.d1_2 = 2 * a.d1_2;
    t.d1_2 = recip * (b * (r.v * s.d1_2 + r.d1_1 * s.d1_1) - (a.d1_1 * t.d1_1));
    r.d1_2 = t.d1_2 / 2;
    s.d2_1 = 1 * a.d2_1;
    t.d2_1 = recip * (b * (r.v * s.d2_1) - (0.0));
    r.d2_1 = t.d2_1 / 1;
    s.d2_2 = 2 * a.d2_2;
    t.d2_2 = recip * (b * (r.v * s.d2_2 + r.d2_1 * s.d2_1) - (a.d2_1 * t.d2_1));
    r.d2_2 = t.d2_2 / 2;
    s.d3_1 = 1 * a.d3_1;
    t.d3_1 = recip * (b * (r.v * s.d3_1) - (0.0));
    r.d3_1 = t.d3_1 / 1;
    s.d3_2 = 2 * a.d3_2;
    t.d3_2 = recip * (b * (r.v * s.d3_2 + r.d3_1 * s.d3_1) - (a.d3_1 * t.d3_1));
    r.d3_2 = t.d3_2 / 2;
    s.d4_1 = 1 * a.d4_1;
    t.d4_1 = recip * (b * (r.v * s.d4_1) - (0.0));
    r.d4_1 = t.d4_1 / 1;
    s.d4_2 = 2 * a.d4_2;
    t.d4_2 = recip * (b * (r.v * s.d4_2 + r.d4_1 * s.d4_1) - (a.d4_1 * t.d4_1));
    r.d4_2 = t.d4_2 / 2;
    s.d5_1 = 1 * a.d5_1;
    t.d5_1 = recip * (b * (r.v * s.d5_1) - (0.0));
    r.d5_1 = t.d5_1 / 1;
    s.d5_2 = 2 * a.d5_2;
    t.d5_2 = recip * (b * (r.v * s.d5_2 + r.d5_1 * s.d5_1) - (a.d5_1 * t.d5_1));
    r.d5_2 = t.d5_2 / 2;
    s.d6_1 = 1 * a.d6_1;
    t.d6_1 = recip * (b * (r.v * s.d6_1) - (0.0));
    r.d6_1 = t.d6_1 / 1;
    s.d6_2 = 2 * a.d6_2;
    t.d6_2 = recip * (b * (r.v * s.d6_2 + r.d6_1 * s.d6_1) - (a.d6_1 * t.d6_1));
    r.d6_2 = t.d6_2 / 2;
    s.d7_1 = 1 * a.d7_1;
    t.d7_1 = recip * (b * (r.v * s.d7_1) - (0.0));
    r.d7_1 = t.d7_1 / 1;
    s.d7_2 = 2 * a.d7_2;
    t.d7_2 = recip * (b * (r.v * s.d7_2 + r.d7_1 * s.d7_1) - (a.d7_1 * t.d7_1));
    r.d7_2 = t.d7_2 / 2;
    s.d8_1 = 1 * a.d8_1;
    t.d8_1 = recip * (b * (r.v * s.d8_1) - (0.0));
    r.d8_1 = t.d8_1 / 1;
    s.d8_2 = 2 * a.d8_2;
    t.d8_2 = recip * (b * (r.v * s.d8_2 + r.d8_1 * s.d8_1) - (a.d8_1 * t.d8_1));
    r.d8_2 = t.d8_2 / 2;
    s.d9_1 = 1 * a.d9_1;
    t.d9_1 = recip * (b * (r.v * s.d9_1) - (0.0));
    r.d9_1 = t.d9_1 / 1;
    s.d9_2 = 2 * a.d9_2;
    t.d9_2 = recip * (b * (r.v * s.d9_2 + r.d9_1 * s.d9_1) - (a.d9_1 * t.d9_1));
    r.d9_2 = t.d9_2 / 2;
    s.d10_1 = 1 * a.d10_1;
    t.d10_1 = recip * (b * (r.v * s.d10_1) - (0.0));
    r.d10_1 = t.d10_1 / 1;
    s.d10_2 = 2 * a.d10_2;
    t.d10_2 = recip * (b * (r.v * s.d10_2 + r.d10_1 * s.d10_1) - (a.d10_1 * t.d10_1));
    r.d10_2 = t.d10_2 / 2;
    s.d11_1 = 1 * a.d11_1;
    t.d11_1 = recip * (b * (r.v * s.d11_1) - (0.0));
    r.d11_1 = t.d11_1 / 1;
    s.d11_2 = 2 * a.d11_2;
    t.d11_2 = recip * (b * (r.v * s.d11_2 + r.d11_1 * s.d11_1) - (a.d11_1 * t.d11_1));
    r.d11_2 = t.d11_2 / 2;
    s.d12_1 = 1 * a.d12_1;
    t.d12_1 = recip * (b * (r.v * s.d12_1) - (0.0));
    r.d12_1 = t.d12_1 / 1;
    s.d12_2 = 2 * a.d12_2;
    t.d12_2 = recip * (b * (r.v * s.d12_2 + r.d12_1 * s.d12_1) - (a.d12_1 * t.d12_1));
    r.d12_2 = t.d12_2 / 2;
    s.d13_1 = 1 * a.d13_1;
    t.d13_1 = recip * (b * (r.v * s.d13_1) - (0.0));
    r.d13_1 = t.d13_1 / 1;
    s.d13_2 = 2 * a.d13_2;
    t.d13_2 = recip * (b * (r.v * s.d13_2 + r.d13_1 * s.d13_1) - (a.d13_1 * t.d13_1));
    r.d13_2 = t.d13_2 / 2;
    s.d14_1 = 1 * a.d14_1;
    t.d14_1 = recip * (b * (r.v * s.d14_1) - (0.0));
    r.d14_1 = t.d14_1 / 1;
    s.d14_2 = 2 * a.d14_2;
    t.d14_2 = recip * (b * (r.v * s.d14_2 + r.d14_1 * s.d14_1) - (a.d14_1 * t.d14_1));
    r.d14_2 = t.d14_2 / 2;
    s.d15_1 = 1 * a.d15_1;
    t.d15_1 = recip * (b * (r.v * s.d15_1) - (0.0));
    r.d15_1 = t.d15_1 / 1;
    s.d15_2 = 2 * a.d15_2;
    t.d15_2 = recip * (b * (r.v * s.d15_2 + r.d15_1 * s.d15_1) - (a.d15_1 * t.d15_1));
    r.d15_2 = t.d15_2 / 2;
    s.d16_1 = 1 * a.d16_1;
    t.d16_1 = recip * (b * (r.v * s.d16_1) - (0.0));
    r.d16_1 = t.d16_1 / 1;
    s.d16_2 = 2 * a.d16_2;
    t.d16_2 = recip * (b * (r.v * s.d16_2 + r.d16_1 * s.d16_1) - (a.d16_1 * t.d16_1));
    r.d16_2 = t.d16_2 / 2;
    s.d17_1 = 1 * a.d17_1;
    t.d17_1 = recip * (b * (r.v * s.d17_1) - (0.0));
    r.d17_1 = t.d17_1 / 1;
    s.d17_2 = 2 * a.d17_2;
    t.d17_2 = recip * (b * (r.v * s.d17_2 + r.d17_1 * s.d17_1) - (a.d17_1 * t.d17_1));
    r.d17_2 = t.d17_2 / 2;
    s.d18_1 = 1 * a.d18_1;
    t.d18_1 = recip * (b * (r.v * s.d18_1) - (0.0));
    r.d18_1 = t.d18_1 / 1;
    s.d18_2 = 2 * a.d18_2;
    t.d18_2 = recip * (b * (r.v * s.d18_2 + r.d18_1 * s.d18_1) - (a.d18_1 * t.d18_1));
    r.d18_2 = t.d18_2 / 2;
    s.d19_1 = 1 * a.d19_1;
    t.d19_1 = recip * (b * (r.v * s.d19_1) - (0.0));
    r.d19_1 = t.d19_1 / 1;
    s.d19_2 = 2 * a.d19_2;
    t.d19_2 = recip * (b * (r.v * s.d19_2 + r.d19_1 * s.d19_1) - (a.d19_1 * t.d19_1));
    r.d19_2 = t.d19_2 / 2;
    s.d20_1 = 1 * a.d20_1;
    t.d20_1 = recip * (b * (r.v * s.d20_1) - (0.0));
    r.d20_1 = t.d20_1 / 1;
    s.d20_2 = 2 * a.d20_2;
    t.d20_2 = recip * (b * (r.v * s.d20_2 + r.d20_1 * s.d20_1) - (a.d20_1 * t.d20_1));
    r.d20_2 = t.d20_2 / 2;
    s.d21_1 = 1 * a.d21_1;
    t.d21_1 = recip * (b * (r.v * s.d21_1) - (0.0));
    r.d21_1 = t.d21_1 / 1;
    s.d21_2 = 2 * a.d21_2;
    t.d21_2 = recip * (b * (r.v * s.d21_2 + r.d21_1 * s.d21_1) - (a.d21_1 * t.d21_1));
    r.d21_2 = t.d21_2 / 2;
    s.d22_1 = 1 * a.d22_1;
    t.d22_1 = recip * (b * (r.v * s.d22_1) - (0.0));
    r.d22_1 = t.d22_1 / 1;
    s.d22_2 = 2 * a.d22_2;
    t.d22_2 = recip * (b * (r.v * s.d22_2 + r.d22_1 * s.d22_1) - (a.d22_1 * t.d22_1));
    r.d22_2 = t.d22_2 / 2;
    s.d23_1 = 1 * a.d23_1;
    t.d23_1 = recip * (b * (r.v * s.d23_1) - (0.0));
    r.d23_1 = t.d23_1 / 1;
    s.d23_2 = 2 * a.d23_2;
    t.d23_2 = recip * (b * (r.v * s.d23_2 + r.d23_1 * s.d23_1) - (a.d23_1 * t.d23_1));
    r.d23_2 = t.d23_2 / 2;
    s.d24_1 = 1 * a.d24_1;
    t.d24_1 = recip * (b * (r.v * s.d24_1) - (0.0));
    r.d24_1 = t.d24_1 / 1;
    s.d24_2 = 2 * a.d24_2;
    t.d24_2 = recip * (b * (r.v * s.d24_2 + r.d24_1 * s.d24_1) - (a.d24_1 * t.d24_1));
    r.d24_2 = t.d24_2 / 2;
    s.d25_1 = 1 * a.d25_1;
    t.d25_1 = recip * (b * (r.v * s.d25_1) - (0.0));
    r.d25_1 = t.d25_1 / 1;
    s.d25_2 = 2 * a.d25_2;
    t.d25_2 = recip * (b * (r.v * s.d25_2 + r.d25_1 * s.d25_1) - (a.d25_1 * t.d25_1));
    r.d25_2 = t.d25_2 / 2;
    s.d26_1 = 1 * a.d26_1;
    t.d26_1 = recip * (b * (r.v * s.d26_1) - (0.0));
    r.d26_1 = t.d26_1 / 1;
    s.d26_2 = 2 * a.d26_2;
    t.d26_2 = recip * (b * (r.v * s.d26_2 + r.d26_1 * s.d26_1) - (a.d26_1 * t.d26_1));
    r.d26_2 = t.d26_2 / 2;
    s.d27_1 = 1 * a.d27_1;
    t.d27_1 = recip * (b * (r.v * s.d27_1) - (0.0));
    r.d27_1 = t.d27_1 / 1;
    s.d27_2 = 2 * a.d27_2;
    t.d27_2 = recip * (b * (r.v * s.d27_2 + r.d27_1 * s.d27_1) - (a.d27_1 * t.d27_1));
    r.d27_2 = t.d27_2 / 2;
    s.d28_1 = 1 * a.d28_1;
    t.d28_1 = recip * (b * (r.v * s.d28_1) - (0.0));
    r.d28_1 = t.d28_1 / 1;
    s.d28_2 = 2 * a.d28_2;
    t.d28_2 = recip * (b * (r.v * s.d28_2 + r.d28_1 * s.d28_1) - (a.d28_1 * t.d28_1));
    r.d28_2 = t.d28_2 / 2;
  }
