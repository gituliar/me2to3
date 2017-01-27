Off FinalStats;
Off Statistics;

#include me2to3.h # Declarations

* Group A
#ifdef `qp2qpg'
  Local ME2 = A(p1,p2,-p3,-p4,-p1-p2+p3+p4)/4/N^2;
  #include me2to3.h # A
#endif
#ifdef `qP2qPg'
  Local ME2 = A(p1,-p4,-p3,p2,-p1-p2+p3+p4)/4/N^2;
  #include me2to3.h # A
#endif
#ifdef `qQ2Ppg'
  Local ME2 = A(p1,p3,-p2,-p4,-p1+p2-p3+p4)/4/N^2;
  #include me2to3.h # A
#endif
#ifdef `qg2qpP'
  Local ME2 = -A(p1,-p1-p2+p3+p4,-p3,-p4,p2)/4/N^2;
  #include me2to3.h # A
#endif

* Group B
#ifdef `qq2qqg'
  Local ME2 = B(p1,p2,-p3,-p4,-p1-p2+p3+p4)/4/N^2;
  #include me2to3.h # B
  #include me2to3.h # A
#endif
#ifdef `qQ2qQg'
  Local ME2 = B(p1,-p4,-p3,p2,-p1-p2+p3+p4)/4/N^2;
  #include me2to3.h # B
  #include me2to3.h # A
#endif
#ifdef `qg2qqQ'
  Local ME2 = -B(p1,-p1-p2+p3+p4,-p3,-p4,p2)/4/(1-ep)/N/V;
  #include me2to3.h # B
  #include me2to3.h # A
#endif

* Group C
#ifdef `qQ2ggg'
  Local ME2 = C(p1,p2,-p3,-p4,-p1-p2+p3+p4)/4/N^2;
  #include me2to3.h # C
#endif
#ifdef `qg2qgg'
  Local ME2 = -C(p1,-p3,p2,-p4,-p1-p2+p3+p4)/4/(1-ep)/N/V;
  #include me2to3.h # C
#endif
#ifdef `Qg2Qgg'
  Local ME2 = -C(-p3,p1,p2,-p4,-p1-p2+p3+p4)/4/(1-ep)/N/V;
  #include me2to3.h # C
#endif
#ifdef `gg2qQg'
  Local ME2 = C(-p3,-p4,p1,p2,-p1-p2+p3+p4)/4/(1-ep)^2/V^2;
  #include me2to3.h # C
#endif

* Group D
#ifdef `gg2ggg'
  Local ME2 = D(p1,p2,-p3,-p4,-p5);
  #include me2to3.h # D
#endif

#include me2to3.h # Expand

Bracket den,sp,s;
.sort
Format Mathematica;
#write <me2to3.temp.m> "(%E)", ME2;
.end
