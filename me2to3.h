* This file contains 2 -> 3 squared matrix elements in massles QCD taken
* from R.K.Ellis and J.C.Sexton "QCD Radiative Corrections to Parton Parton
* Scattering" 1985.

*--#[ Declarations:
CF A,B,C,D, Fc,Fd, den, sp, f,g;
V p,p1,...,p5;
S s,N,V,V1,...,V4,ep,[1-ep],[3+ep];

S d1,d2,d3;

#define s  "sp(p1+p2)"
#define t  "sp(p1+p3)"
#define u  "sp(p1+p4)"
#define sr "sp(p3+p4)"
#define tr "sp(p2+p4)"
#define ur "sp(p2+p3)"

#define sl "((sp(p1+p2)+sp(p3+p4))/2)"
#define sm "((sp(p1+p2)-sp(p3+p4))/2)"
#define tl "((sp(p1+p3)+sp(p2+p4))/2)"
#define tm "((sp(p1+p3)-sp(p2+p4))/2)"
#define ul "((sp(p1+p4)+sp(p2+p3))/2)"
#define um "((sp(p1+p4)-sp(p2+p3))/2)"

#define ap "((sp(p1+p3)+sp(p2+p3))/2)"
#define bp "((sp(p1+p4)+sp(p2+p4))/2)"
#define cp "((sp(p1+p5)+sp(p2+p5))/2)"
#define am "((sp(p1+p3)-sp(p2+p3))/2)"
#define bm "((sp(p1+p4)-sp(p2+p4))/2)"
#define cm "((sp(p1+p5)-sp(p2+p5))/2)"

#define s12 "(sp(p1+p2)/2)"
#define s13 "(sp(p1+p3)/2)"
#define s14 "(sp(p1+p4)/2)"
#define s15 "(sp(p1+p5)/2)"
#define s23 "(sp(p2+p3)/2)"
#define s24 "(sp(p2+p4)/2)"
#define s25 "(sp(p2+p5)/2)"
#define s34 "(sp(p3+p4)/2)"
#define s35 "(sp(p3+p5)/2)"
#define s45 "(sp(p4+p5)/2)"
*--#] Declarations:

*--#[ A:
id A(p1?,p2?,p3?,p4?,p5?) = 1/`t'/`tr'/`s15'/`s25'/`s35'/`s45'*(
   + V1*(`sl'^2+`sm'^2+`ul'^2+`um'^2)*(1/2*`ul'*(`s'*`sr'+`t'*`tr'-`u'*`ur')+1/4*`u'*(`s'*`t'+`sr'*`tr')+1/4*`ur'*(`s'*`tr'+`sr'*`t'))
   - V2*(`sl'^2+`sm'^2+`ul'^2+`um'^2)*(1/2*`sl'*(`s'*`sr'-`t'*`tr'-`u'*`ur')+`ul'*`t'*`tr'+`tl'*`u'*`ur')
   + ep*V1*((`ul'*(`sm'^2+`um'^2+1/2*`tl'^2)+1/4*`tl'*(1-2*ep)*(`sm'^2-`tm'^2-`um'^2))*(`sm'^2+`tm'^2-`um'^2)-`sm'*`tm'*`um'*(2*`sm'^2+2*`um'^2+`tl'^2)-`tm'^3*`sm'*`um'*(1-2*ep))
   + ep*V2*(-1/2*`tl'^2*(`sl'*(`sm'^2-`tm'^2-`um'^2)+2*`ul'*`tm'^2+2*`tl'*`um'^2)+1/4*`tl'*(1-2*ep)*(`sm'^4+`tm'^4+`um'^4-2*`sm'^2*`tm'^2-2*`sm'^2*`um'^2-2*`tm'^2*`um'^2)-`sl'*(`sm'^4-`um'^4)-2*`tl'*`um'^4+(`sl'-2*`ul')*`tm'^2*(`sm'^2+`um'^2)-2*`tl'*`sm'^2*`um'^2)
);
*--#] A:

*--#[ B:
id B(p1?,p2?,p3?,p4?,p5?) =
   + A(p1,p2,p3,p4,p5)
   + A(p1,p2,p4,p3,p5)
   + 1/`s15'/`s25'/`s35'/`s45'/`t'/`tr'/`u'/`ur'*(
     + (`s'*`sr'-`t'*`tr'-`u'*`ur')*(`s'^2+`sr'^2)*V3*(
       + 1/4*`sl'*(`s'*`sr'-`t'*`tr'-`u'*`ur')
       + 1/2*`ul'*`t'*`tr'+1/2*`tl'*`u'*`ur')
     + (`s'*`sr'-`t'*`tr'-`u'*`ur')*(`s'^2+`sr'^2)*V4*(
       + 1/4*`sl'*(`s'*`sr'-`t'*`tr'-`u'*`ur')
       - 1/2*`ul'*`t'*`tr'
       - 1/2*`tl'*`u'*`ur'
       - 1/4*`s'*(`t'*`u'+`tr'*`ur')
       - 1/4*`sr'*(`t'*`ur'+`u'*`tr'))
     - ep*(`s'*`sr'-`t'*`tr'-`u'*`ur')*(
       - 1/4*(5+ep)*(V3+V4)*`sl'*(`sm'^2-`tm'^2-`um'^2)^2
       + (`sm'^2-`tm'^2-`um'^2)*(
         + (1/2*(1+ep)*(1-2*ep)*(V3+V4)-1/2*(5+ep)*(V3-V4))*(`tm'^2*`ul'+`um'^2*`tl')
         - (1-ep)*(V3+V4)*`sl'*(`tm'^2+`um'^2)
         - 1/2*(V3+V4)*`sl'*(`tl'^2+`ul'^2+`tm'^2+`um'^2)
         - (5+ep)*V4*`sm'*`tm'*`um')
       - 2*V4*`sm'*`tm'*`um'*(
         + `tl'^2+`ul'^2+`tm'^2+`um'^2
         + 2*(1-ep)*(`tm'^2+`um'^2)
         + (1+ep)*(1-2*ep)*`tl'*`ul')
       + (1+ep)*(V3+V4)*`sl'*`t'*`tr'*`u'*`ur'
       - (V3-V4)*(`tm'^2*`ul'+`um'^2*`tl')*(`tl'^2+`ul'^2+`tm'^2+`um'^2)
       - 2*(1-ep)*(V3-V4)*(`tm'^2*`ul'+`um'^2*`tl')*(`tm'^2+`um'^2)
       + (1+ep)*(1-2*ep)*(V3-V4)*`sl'*`tm'^2*`um'^2)
     + ep*(1+ep)*`t'*`tr'*`u'*`ur'*(
       + 2*V3*(`sl'*`tl'*`ul'+`tl'*`um'^2+`ul'*`tm'^2)
       + 2*V4*(`sl'*`tl'*`ul'-`tl'*`um'^2-`ul'*`tm'^2)
       + 8*ep*V4*`sm'*`tm'*`um'))
;
*--#] B:

*--#[ C:
id C(p1?,p2?,p3?,p4?,p5?) =
   + Fc(p1,p2, p3,p4,p5)
   + Fc(p1,p2, p3,p5,p4)
   + Fc(p1,p2, p4,p3,p5)
   + Fc(p1,p2, p4,p5,p3)
   + Fc(p1,p2, p5,p3,p4)
   + Fc(p1,p2, p5,p4,p3)

   + Fc(p2,p1, p3,p4,p5)
   + Fc(p2,p1, p3,p5,p4)
   + Fc(p2,p1, p4,p3,p5)
   + Fc(p2,p1, p4,p5,p3)
   + Fc(p2,p1, p5,p3,p4)
   + Fc(p2,p1, p5,p4,p3)
;

id Fc(p1?,p2?,p3?,p4?,p5?) = V*N^-2*(
   - N^2*[1-ep]^2*(`ap'*`bp'-`am'*`bm')/`s34'*
       (2*(`ap'^2+`am'^2)/`s14'/`s24'/`s15'/`s25' + (`cp'^2+`cm'^2)/`s13'/`s23'/`s14'/`s24')
   + N^4*[1-ep]*(`ap'*`bp'-`am'*`bm')/`s12'/`s45'/`s35'*
       (2*(`bp'^2+`bm'^2)/`s13'/`s23' + (`cp'^4-`cm'^4)/`s14'/`s24'/`s13'/`s23')
   + (N^2+1)/4*`s12'/`s13'/`s23'/`s14'/`s24'*
       (2*[1-ep]^2*(`cp'^2+`cm'^2) - d1*(`ap'^2+`am'^2+`bp'^2+`bm'^2-`cp'^2-`cm'^2)
       - 2*d2*(`cp'^2-`cm'^2) - d3*(`ap'^2-`am'^2+`bp'^2-`bm'^2-`cp'^2+`cm'^2))
   + N^2*(-ep*(1+3*ep)*(`ap'^2-`am'^2+`bp'^2-`bm'^2)*
       (`s13'*`s14'*`s15'-`s23'*`s24'*`s25')*`cm'/2)/`s34'/`s13'/`s23'/`s14'/`s24'/`s15'/`s25'
   + 2*N^2*d2*(`ap'*`bp'-`am'*`bm')/`s15'/`s25'/`s34'
   + N^2*(1/2*(d1-2*d2+d3)*(`ap'+`bp'-`cp')*(`ap'+`bp'+`cp')*
       ((`ap'+`bp')*(`ap'+`bp'-`cp')-`ap'*`bp'-3*`am'*`bm')+2*d2*`ap'*`bp'*(`ap'^2+`bp'^2)
       + (2*d3-6*d2)*`am'*`bm'*(`ap'^2+`bp'^2) + (d3-4*d2)*`ap'*`bp'*`cm'^2
       + d3*`am'*`bm'*(`am'^2+`bm'^2) + (2*d2-d3)*`cp'*(`ap'*`bm'^2+`am'^2*`bp')
       + (6*d2-3*d3)*`am'*`bm'*(`ap'+`bp')*`cp' + (6*d3-8*d2)*`ap'*`bp'*`am'*`bm')/`s34'/`s13'/`s23'/`s14'/`s24'
   + N^4*(-2*d2*(`ap'*`bp'-`am'*`bm')*`cp'^4 - d2*(`ap'+`bp')^3*`cp'*`ap'*`bp'
       + 2*d2*(`ap'*`bp'-`am'*`bm')*`cp'^2*`cm'^2 + d2*(`ap'*`bp'+`am'*`bm')*(`ap'+`bp')*`cp'^3
       + d2*(`ap'*`bm'+`am'*`bp')*`cp'^3*`cm' - d2*`ap'*`bp'*(`ap'^2-`am'^2+`bp'^2-`bm'^2)*`cp'^2
       + 4*d2*(`ap'^2-`am'^2)*(`bp'^2-`bm'^2)*`cp'^2 + 2*d2*(`ap'*`bp'+`am'*`bm')*(`ap'*`bp'-`am'*`bm')*`cp'^2
       + d2*`cp'*(`ap'+`bp')*(`ap'*`bp'-`am'*`bm')^2
       - 2*d2*(`ap'*`bp'-`am'*`bm')*(`ap'*`bm'-`am'*`bp')*`cp'*`cm' - d2*(`ap'*`bp'+`am'*`bm')*(`ap'+`bp')*`cp'*`am'*`bm'
       - d2*(`ap'*`bm'+`am'*`bp')*`ap'*`bp'*`cp'*`cm' + d2*`ap'*`bp'*(`ap'+`bp')*`cp'*`cm'^2
       - 2*d2*(`ap'^2-`am'^2+`bp'^2-`bm'^2)*`ap'*`bp'*`am'*`bm' + 2*d2*(`ap'*`bp'-`am'*`bm')*(`ap'^2+`bp'^2)*(`am'^2+`bm'^2)
       - 2*d2*(`ap'*`bp'-`am'*`bm')*(`ap'^2+`bp'^2)^2 + d2*`ap'^4*(`bp'^2-`bm'^2)
       + d2*`bp'^4*(`ap'^2-`am'^2) + d2*`ap'*`bp'*`ap'^2*(`ap'^2-`am'^2) + d2*`ap'*`bp'*`bp'^2*(`bp'^2-`bm'^2)
       + 1/2*(d1+d3)*(`ap'*`bp'-`am'*`bm')*(`ap'^4+`bp'^4+`cp'^4-2*`ap'^2*`bp'^2-2*`ap'^2*`cp'^2-2*`bp'^2*`cp'^2)
       - d3*(`ap'*`bp'-`am'*`bm')*(`ap'^2*`am'^2+`bp'^2*`bm'^2+`cp'^2*`cm'^2)
       + d3*(`ap'*`bp'-`am'*`bm')*(`ap'^2*`bm'^2+`am'^2*`bp'^2+`ap'^2*`cm'^2+`am'^2*`cp'^2+`bp'^2*`cm'^2+`bm'^2*`cp'^2)
       )/`s12'/`s35'/`s45'/`s13'/`s23'/`s14'/`s24'
);
*--#] C:

*--#[ D:
id D(p1?,p2?,p3?,p4?,p5?) = 1/10*(
     f(p1,p2,p3,p4,p5)
);
id f(?a) = distrib_(1,1,Fd,f,?a);
#do i=1,4
  id f(?a) = distrib_(1,1,g,f,?a);
  id Fd(?a)*g(p1?) = Fd(?a, p1);
#enddo
id f() = 1;

id Fd(p1?,p2?,p3?,p4?,p5?) = V*N^3/`s12'/`s23'/`s34'/`s45'/`s15'*(
  + 1/2*[1-ep]^2*(`s12'^4+`s13'^4+`s14'^4+`s15'^4+`s23'^4+`s24'^4+`s25'^4+`s34'^4+`s35'^4+`s45'^4)
  - 3*d3*(`s12'^2*`s23'^2+`s23'^2*`s34'^2+`s34'^2*`s45'^2+`s45'^2*`s15'^2+`s15'^2*`s12'^2)
  + 6*d3*(`s12'*`s23'^2*`s34'+`s23'*`s34'^2*`s45'+`s34'*`s45'^2*`s15'+`s45'*`s15'^2*`s12'+`s15'*`s12'^2*`s23')
  - 6*d3*(`s12'*`s23'*`s34'*`s45'+`s23'*`s34'*`s45'*`s15'+`s34'*`s45'*`s15'*`s12'+`s45'*`s15'*`s12'*`s23'+`s15'*`s12'*`s23'*`s34')
);
*--#] D:

*--#[ Expand:
*id V = N^2 - 1;
id V1 = V^2/N;
id V2 = V/N;
id V3 = V*(N^2+1)/(2*N^2);
id V4 = V^2/(2*N^2);

id d1 = ep*[1-ep]^2;
id d2 = ep*[1-ep];
id d3 = ep*[3+ep];

id [1-ep] = 1-ep;
id [1-ep]^-1 = 1/(1-ep);
id [3+ep] = 3+ep;

id   sp(p1+p2) =   s;
id 1/sp(p1+p2) = 1/s;
#do i=1,5
  #do j=`i'+1,5
    id   sp(-p`i'+p`j') =   sp(p`i'-p`j');
    id   sp(-p`i'-p`j') =   sp(p`i'+p`j');
    id 1/sp(-p`i'-p`j') = 1/sp(p`i'+p`j');
    #do k=`j'+1,5
      id sp(-p`i'+p`j'+p`k')    =   sp(p`i'-p`j'-p`k');
      id sp(-p`i'+p`j'+p`k')^-1 = 1/sp(p`i'-p`j'-p`k');
      id sp(-p`i'+p`j'-p`k')    =   sp(p`i'-p`j'+p`k');
      id sp(-p`i'+p`j'-p`k')^-1 = 1/sp(p`i'-p`j'+p`k');
      id sp(-p`i'-p`j'+p`k')    =   sp(p`i'+p`j'-p`k');
      id sp(-p`i'-p`j'+p`k')^-1 = 1/sp(p`i'+p`j'-p`k');
    #enddo
  #enddo
#enddo
.sort

repeat;
id 1/sp(?args)*sp(?args) = 1;
id sp(?args)*1/sp(?args) = 1;
endrepeat;
.sort
id sp(?args)^-1 = den(sp(?args));
.sort
*--#] Expand:
