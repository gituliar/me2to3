#!/usr/bin/env math-script

(* This script verifies if a given matrix element squared (ME2) is free of
errors. To do this, it reads a list of phase-space points and numerical values
of the ME2 generated by KaTie and compares them with values of the analytical
ME2 in these points. *)

If[
  Length[$CommandLine] < 5
  ,
  Print["Usage:"];
  Print["  comapre.ms <expr> <pts>"];
  Abort[];
];

$me2 = Get[$CommandLine[[4]]];
$pts = Get[$CommandLine[[5]]];

sp[{E_, px_, py_, pz_}] := E^2 - px^2 - py^2 - pz^2;

$N[i_] := $me2 /. {
            ep -> 0,
            N -> 3,
            s -> sp[p1 + p2], 
            den[ex_] :> 1/ex
          } /. {
            p1 -> -$pts[[i]][[1]], 
            p2 -> -$pts[[i]][[2]],
            p3 ->  $pts[[i]][[3]], 
            p4 ->  $pts[[i]][[4]]};

$n = $pts[[1]][[6]]/$N[1];

For[$i=1, $i<=Length[$pts], $i++,
  $real     = $n $N[$i];
  $expected = $pts[[$i]][[6]];
  If[
    TrueQ[Abs[($real - $expected)/$real] > 10^-9]
    ,
    Print["Wrong value at point "<>ToString[$i]];
    Print["  Got:      "<>ToString[$real, InputForm]];
    Print["  Expected: "<>ToString[$expected, InputForm]];
    Abort[];
  ];
];
Print["OK"];