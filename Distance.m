(* ::Package:: *)

(* Input directory *)
dirRoot = FileNameJoin[{DirectoryName[$InputFileName], "Data"}];


(* ::Item:: *)
(*Weights*)


Unprotect[\[Rho]t, \[Rho]x];
{\[Rho]t, \[Rho]x} = {1 + 10^-20, 11/10};
Protect[\[Rho]t, \[Rho]x];


(* ::Item:: *)
(*vHat (Fourier coefficients of the approximate solution v)*)


MapIndexed[(vHat[First[#2]] = Import[FileNameJoin[{dirRoot, #1, "vHat.m"}]])&, 
  {"1st", "2nd", "3rd"}
];


(* ::Item:: *)
(*Bound on norm of A*)


bndNormA[1] = 19906846310/1894822989;
bndNormA[2] = 70711511471/1600927163;
bndNormA[3] = 67312911313/1520322491;


(* ::Item:: *)
(*\[Delta]*)


delta[1] = 13081330001308133/209065086978900000000000;
delta[2] = 263430000026343/141217277216243600000000;
delta[3] = 9160870000916087/571871816392180000000000;


(* ::Item:: *)
(*\[Epsilon]*)


Do[epsilon[i] = delta[i]*bndNormA[i], {i, 1, 3}];


(* ::Item:: *)
(*Distance*)


Print[
  Function[{i, j},
    Module[{M, N, vplus, vminus},
      M = Max[First /@ Dimensions /@ {vHat[i], vHat[j]}];
      N = Max[Last /@ Dimensions /@ {vHat[i], vHat[j]}];
      vplus = PadRight[vHat[i], {M, N}] + PadRight[vHat[j], {M, N}];
      vminus = PadRight[vHat[i], {M, N}] - PadRight[vHat[j], {M, N}];

      Sum[Abs[#[[m + 1, n + 1]]] * \[Rho]t^(2 m + 1) * \[Rho]x^(2 n + 1),
          {m, 0, M - 1},
          {n, 0, N - 1}
        ] > epsilon[i] + epsilon[j]& /@ {vplus, vminus}
    ]
  ] @@@ {{1, 2}, {1, 3}, {2, 3}}
];
