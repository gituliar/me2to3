Author: Oleksandr Gituliar <oleksandr@gituliar.net>
Last modified: 2017-02-02

This directory contains QCD squared matrix elements for 2 → 3 processes.
It is assumed that
  * pᵢ² = 0, i.e., all external particles are massless
  * n-dimensional space-time, n=4-2ε
  * physical kinematics p₁+p₂=p₃+p₄+p₅
  * s = (p₁+p₂)²

Initial expressions were taken from [1] and processed with the help of FORM.
In total we have 12 physical processes which according to [1, Table 7] are
divided into 4 groups:

    Group A        Group B        Group C        Group D
    -------        -------        -------        -------
    qp2qpg         qq2qqg         qQ2ggg         gg2ggg
    qP2qPg         qQ2qQg         qg2qgg
    qQ2Ppg         qg2qqQ         Qg2Qgg
    qg2qpP                        gg2qQg

where we use the notation:
  * q is a j-type quark
  * p is a k-type quark
  * Q is a j-type qnti-quark
  * P is a k-type qnti-quark
  * g is a gluon

[1] R.K.Ellis and J.C.Sexton "QCD Radiative Corrections to Parton Parton
Scattering" 1985
