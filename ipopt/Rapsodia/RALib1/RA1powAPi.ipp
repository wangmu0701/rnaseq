// This file was generated by Rapsodia (see www.mcs.anl.gov/Rapsodia)
r.v = pow(double(a.v), double(b));
if (a.v == 0.0)
  {
    if (b <= 0)
      {
        r.d1_1 = makeFPE1(0.0, 0.0);
        r.d2_1 = makeFPE1(0.0, 0.0);
        r.d3_1 = makeFPE1(0.0, 0.0);
        r.d4_1 = makeFPE1(0.0, 0.0);
        r.d5_1 = makeFPE1(0.0, 0.0);
        r.d6_1 = makeFPE1(0.0, 0.0);
        r.d7_1 = makeFPE1(0.0, 0.0);
      }
    else
      {
        {
          if (b == 1)
            {
              r.d1_1 = a.d1_1;
              r.d2_1 = a.d2_1;
              r.d3_1 = a.d3_1;
              r.d4_1 = a.d4_1;
              r.d5_1 = a.d5_1;
              r.d6_1 = a.d6_1;
              r.d7_1 = a.d7_1;
            }
          else
            {
              r.d1_1 = 0.0;
              r.d2_1 = 0.0;
              r.d3_1 = 0.0;
              r.d4_1 = 0.0;
              r.d5_1 = 0.0;
              r.d6_1 = 0.0;
              r.d7_1 = 0.0;
              for(j=3;j<=int(b);j+=1)
                {
                  r.d1_1 = r.v * a.v;
                  r.d2_1 = r.v * a.v;
                  r.d3_1 = r.v * a.v;
                  r.d4_1 = r.v * a.v;
                  r.d5_1 = r.v * a.v;
                  r.d6_1 = r.v * a.v;
                  r.d7_1 = r.v * a.v;
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
    s.d2_1 = 1 * a.d2_1;
    t.d2_1 = recip * (b * (r.v * s.d2_1) - (0.0));
    r.d2_1 = t.d2_1 / 1;
    s.d3_1 = 1 * a.d3_1;
    t.d3_1 = recip * (b * (r.v * s.d3_1) - (0.0));
    r.d3_1 = t.d3_1 / 1;
    s.d4_1 = 1 * a.d4_1;
    t.d4_1 = recip * (b * (r.v * s.d4_1) - (0.0));
    r.d4_1 = t.d4_1 / 1;
    s.d5_1 = 1 * a.d5_1;
    t.d5_1 = recip * (b * (r.v * s.d5_1) - (0.0));
    r.d5_1 = t.d5_1 / 1;
    s.d6_1 = 1 * a.d6_1;
    t.d6_1 = recip * (b * (r.v * s.d6_1) - (0.0));
    r.d6_1 = t.d6_1 / 1;
    s.d7_1 = 1 * a.d7_1;
    t.d7_1 = recip * (b * (r.v * s.d7_1) - (0.0));
    r.d7_1 = t.d7_1 / 1;
  }
