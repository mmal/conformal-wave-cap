(* ::Package:: *)

(* ::Subsection:: *)
(*Code*)


(* ::Subsubsection:: *)
(*Basic*)


(* ::Item::Closed:: *)
(*coeffsSquareP*)


Clear[coeffsSquareP];

coeffsSquareP[uhat_List, {MM_Integer, NN_Integer}] := 
  Module[{u = uhat, vHat},
   vHat = 1/4 * 
     ParallelTable[
      (2 - KroneckerDelta[p, 0]) * (
        Sum[u[[m + p + 1, n + q + 1]] u[[m + 1, n + 1]], 
          {m, 0, MM - 1 - p}, 
          {n, 0, NN - 1 - q}
        ] + 
        (1 - KroneckerDelta[q, 0]) * 
          Sum[u[[m + p + 1, n - q + 1]] u[[m + 1, n + 1]], 
            {m, 0, MM - 1 - p}, 
            {n, q, NN - 1}
          ]
      ) + 
      (2 - KroneckerDelta[p, 0]) * 
        Sum[u[[m + p + 1, n + 1]] u[[m + 1, q - 1 - n + 1]], 
          {m, 0, MM - 1 - p}, 
          {n, Max[0, q - NN], Min[NN - 1, q - 1]}
        ] + 
      (2 - KroneckerDelta[q, 0]) * 
        Sum[u[[m + 1, n + q + 1]] u[[p - 1 - m + 1, n + 1]], 
          {m, Max[0, p - MM], Min[MM - 1, p - 1]}, 
          {n, 0, NN - 1 - q}
        ] + 
      Sum[u[[m + 1, n + 1]] u[[p - 1 - m + 1, q - 1 - n + 1]], 
        {m, Max[0, p - MM], Min[MM - 1, p - 1]}, 
        {n, Max[0, q - NN], Min[NN - 1, q - 1]}
      ],
      {p, 0, 2 MM - 1}, 
      {q, 0, 2 NN - 1},
      Method -> "FinestGrained"
    ];

   vHat
  ];


(* ::Item::Closed:: *)
(*coeffsCubeP*)


Clear[coeffsCubeP];

coeffsCubeP[uhat_List, vHat_List, {MM_Integer, NN_Integer}] := 
  Module[{u = uhat, v = vHat, what},
   what = 1/4 * 
     ParallelTable[
      Sum[
        Sum[
          v[[m - k + 1, n - l + 1]] u[[m + 1, n + 1]], 
          {m, k, MM - 1}
        ] + 
        Sum[
          v[[m + k + 1 + 1, n - l + 1]] u[[m + 1, n + 1]], 
          {m, 0, Min[MM - 1, 2 MM - 2 - k]}
        ],
        {n, l, NN - 1}
      ] + 
      Sum[
        Sum[
          v[[m - k + 1, n + l + 1 + 1]] u[[m + 1, n + 1]], 
          {m, k, MM - 1}
        ] + 
        Sum[
          v[[m + k + 1 + 1, n + l + 1 + 1]] u[[m + 1, n + 1]], 
          {m, 0, Min[MM - 1, 2 MM - 2 - k]}
        ],
        {n, 0, Min[NN - 1, 2 NN - 2 - l]}
      ] + 
      Sum[
        Sum[
          v[[m - k + 1, l - n + 1]] u[[m + 1, n + 1]], 
          {m, k, MM - 1}
        ] + 
        Sum[
          v[[m + k + 1 + 1, l - n + 1]] u[[m + 1, n + 1]], 
          {m, 0, Min[MM - 1, 2 MM - 2 - k]}
        ],
        {n, Max[0, l - 2 NN + 1], Min[NN - 1, l]}
      ] + 
      Sum[
        Sum[
          v[[k - m + 1, n - l + 1]] u[[m + 1, n + 1]], 
          {n, l, NN - 1}
        ] + 
        Sum[
          v[[k - m + 1, n + l + 1 + 1]] u[[m + 1, n + 1]], 
          {n, 0, Min[NN - 1, 2 NN - 2 - l]}
        ],
        {m, Max[0, k - 2 MM + 1], Min[MM - 1, k]}
      ] + 
      Sum[
        v[[k - m + 1, l - n + 1]] u[[m + 1, n + 1]], 
        {m, Max[0, k - 2 MM + 1], Min[MM - 1, k]}, 
        {n, Max[0, l - 2 NN + 1], Min[NN - 1, l]}
      ],
      {k, 0, 3 MM - 2}, 
      {l, 0, 3 NN - 2}, 
      Method -> "FinestGrained"
    ];

   what
  ];


(* ::Item::Closed:: *)
(*coeffC*)


Clear[coeffC];

coeffC[\[Mu]_Integer, \[Nu]_Integer, m_Integer, n_Integer, a2_?MatrixQ, MM_Integer, NN_Integer] := 
  1/4 * (
    If[0 <= m - \[Mu] <= -1 + 2 MM && 0 <= n - \[Nu] <= -1 + 2 NN, 
      a2[[ (m - \[Mu]) + 1, (n - \[Nu]) + 1 ]], 
      0] + 
    If[0 <= m - \[Mu] <= -1 + 2 MM && 0 <= -n + \[Nu] <= -1 + 2 NN, 
      a2[[ (m - \[Mu]) + 1, (-n + \[Nu]) + 1 ]], 
      0] + 
    If[0 <= m - \[Mu] <= -1 + 2 MM && 0 <= 1 + n + \[Nu] <= -1 + 2 NN, 
      a2[[ (m - \[Mu]) + 1, (1 + n + \[Nu]) + 1 ]], 
      0] + 
    If[0 <= -m + \[Mu] <= -1 + 2 MM && 0 <= n - \[Nu] <= -1 + 2 NN, 
      a2[[ (-m + \[Mu]) + 1, (n - \[Nu]) + 1 ]], 
      0] + 
    If[0 <= -m + \[Mu] <= -1 + 2 MM && 0 <= -n + \[Nu] <= -1 + 2 NN, 
      a2[[ (-m + \[Mu]) + 1, (-n + \[Nu]) + 1 ]], 
      0] + 
    If[0 <= -m + \[Mu] <= -1 + 2 MM && 0 <= 1 + n + \[Nu] <= -1 + 2 NN, 
      a2[[ (-m + \[Mu]) + 1, (1 + n + \[Nu]) + 1 ]], 
      0] + 
    If[0 <= 1 + m + \[Mu] <= -1 + 2 MM && 0 <= n - \[Nu] <= -1 + 2 NN, 
      a2[[ (1 + m + \[Mu]) + 1, (n - \[Nu]) + 1 ]], 
      0] + 
    If[0 <= 1 + m + \[Mu] <= -1 + 2 MM && 0 <= -n + \[Nu] <= -1 + 2 NN, 
      a2[[ (1 + m + \[Mu]) + 1, (-n + \[Nu]) + 1 ]], 
      0] + 
    If[0 <= 1 + m + \[Mu] <= -1 + 2 MM && 0 <= 1 + n + \[Nu] <= -1 + 2 NN, 
      a2[[ (1 + m + \[Mu]) + 1, (1 + n + \[Nu]) + 1 ]], 
      0]
  );


(* ::Item::Closed:: *)
(*cmpNormH0Pmn*)


Clear[cmpNormH0Pmn];

cmpNormH0Pmn[{m_Integer, n_Integer}, {MM_Integer, NN_Integer}, a2_?MatrixQ, \[CapitalOmega]_?NumericQ, {\[Rho]t_, \[Rho]x_}] := 
  Module[{\[Alpha], \[Beta], c, d, z},
   \[Alpha][mm_, nn_] := -\[CapitalOmega]^2 (2 mm + 1)^2 + 4 (nn + 1)^2;
   
   \[Beta][mm_, nn_] := -((4 (2 nn + 1)) / \[Alpha][mm, nn - 1]);
   
   c[\[Mu]_, \[Nu]_] := coeffC[\[Mu], \[Nu], m, n, a2, MM, NN] / \[Alpha][\[Mu], \[Nu]];
   
   d[\[Mu]_, n + 2 NN - 1] := d[\[Mu], n + 2 NN - 1] = 0;
   d[\[Mu]_, n + 2 NN - 2] := d[\[Mu], n + 2 NN - 2] = c[\[Mu], n + 2 NN - 1] * \[Beta][\[Mu], n + 2 NN - 1];
   d[\[Mu]_, \[Nu]_] := d[\[Mu], \[Nu]] = c[\[Mu], \[Nu] + 1] * \[Beta][\[Mu], \[Nu] + 1] + d[\[Mu], \[Nu] + 1];
   
   z[\[Mu]_, \[Nu]_] := c[\[Mu], \[Nu]] + d[\[Mu], \[Nu]];
   
   Sum[Abs[-3 z[\[Mu], \[Nu]]] * \[Rho]t^(2 \[Mu] + 1) * \[Rho]x^(2 \[Nu] + 1), 
     {\[Mu], 0, m + 2 MM - 1}, 
     {\[Nu], n + 2 NN - 1, 0, -1}] / (\[Rho]t^(2 m + 1) * \[Rho]x^(2 n + 1))
  ];


(* ::Item::Closed:: *)
(*Coefficients z_{\mu,\nu}^{m,n}*)


Clear[coeffZ];

coeffZ[{m_Integer, n_Integer}, {MM_Integer, NN_Integer}, a2_?MatrixQ, \[CapitalOmega]_?NumericQ] :=
  Module[{\[Alpha], \[Beta], c, d, z},
   \[Alpha][mm_, nn_] := -\[CapitalOmega]^2 (2mm+1)^2 + 4(nn+1)^2;
   \[Beta][mm_, nn_] := -((4(2nn+1))/\[Alpha][mm,nn-1]);
   c[\[Mu]_, \[Nu]_] := coeffC[\[Mu], \[Nu], m, n, a2, MM, NN]/\[Alpha][\[Mu], \[Nu]];

   d[\[Mu]_, n+2NN-1] := d[\[Mu], n+2NN-1] = 0;
   d[\[Mu]_, n+2NN-2] := d[\[Mu], n+2NN-2] = c[\[Mu], n+2NN-1]\[Beta][\[Mu], n+2NN-1];
   d[\[Mu]_, \[Nu]_] := d[\[Mu],\[Nu]] = c[\[Mu], \[Nu]+1]\[Beta][\[Mu], \[Nu]+1] + d[\[Mu], \[Nu]+1];

   z[\[Mu]_, \[Nu]_] := c[\[Mu],\[Nu]] + d[\[Mu],\[Nu]];

   Reverse/@Table[z[\[Mu], \[Nu]], {\[Mu], 0, m+2MM-1}, {\[Nu], n+2NN-1, 0, -1}]
];


(* ::Item::Closed:: *)
(*Function g*)


Clear[funG];

funG[{m_Integer, n_Integer}, c_?MatrixQ, {Mhat_Integer, Nhat_Integer}] := Module[{\[Alpha]},
  \[Alpha][mm_, nn_] := -\[CapitalOmega]0^2 (2 mm + 1)^2 + 4 (nn + 1)^2;

  3*Sum[
    Abs[c[[mhat + 1, nhat + 1]]]*(
      Abs[1/\[Alpha][mhat + m - Mhat, nhat + n - Nhat]] +
      Abs[
        (4 (2 (nhat + n - Nhat) + 1))/(
          \[Alpha][mhat + m - Mhat, nhat + n - Nhat - 1]
          \[Alpha][mhat + m - Mhat, nhat + n - Nhat]
        )
      ]*(1 - \[Rho]x^(-2 (nhat + n - Nhat)))/(\[Rho]x^2 - 1)
    )*\[Rho]t^(2 (mhat - Mhat)) \[Rho]x^(2 (nhat - Nhat)),
    {mhat, 0, 2 Mhat},
    {nhat, 0, 2 Nhat}
  ]
];


(* ::Subsubsection:: *)
(*Tools*)


(* ::Item::Closed:: *)
(*Timing*)


$scriptStartTime = SessionTime[];


(* ::Item::Closed:: *)
(*Sets Sk*)


Clear[pointsInSk];

pointsInSk[-1, {A_, B_, CC_}, {M_Integer, NN_Integer}, {dM_Integer, dN_Integer}] := 
  Join[
    Flatten[Table[{m, n}, {m, 0, M - 1}, {n, 0, A * m + B}], 1],
    Flatten[Table[{m, n}, {m, M, (NN + CC)/A}, {n, Ceiling[Max[0, A * m - CC]], NN - 1}], 1]
  ];

pointsInSk[k_Integer, {A_, B_, CC_}, {M_Integer, NN_Integer}, {dM_Integer, dN_Integer}] := 
  Join[
    Flatten[Table[{m, n}, {m, M + k * dM, M + (k + 1) * dM - 1}, {n, NN + k * dN, A * m + B}], 1],
    Flatten[Table[{m, n}, {m, M + (k + 1) * dM, (NN + (k + 1) * dN + CC)/A}, {n, Ceiling[Max[NN + k * dN, A * m - CC]], NN + (k + 1) * dN - 1}], 1]
  ];


(* ::Item::Closed:: *)
(*Printing functions*)


(* Print a numbered tree with one branch per string *)
printTree[lines : {__String}] := Module[{branches},
  branches = MapIndexed[
    If[First[#2] == Length[lines], "\:2514\[HorizontalLine] ", "\:251c\[HorizontalLine] "] <> #1 &,
    lines
  ];
  Print[
    "#", ++$printCount, "  ",
    DateString[],
    "\n",
    StringRiffle[branches, "\n"],
    "\n"
  ]
];


printMessage[message_String] := printTree[{message}];


SetAttributes[printTest, HoldAll];
printTest[
   expr_Less | expr_Greater | expr_Equal |
   expr_LessEqual | expr_GreaterEqual
] := Module[{text, result, status},
  text = ToString[
    Unevaluated[expr], InputForm, PageWidth -> Infinity
  ];
  result = expr;
  status = Which[
    result === True,  "True",
    result === False, "False",
    True, "UNRESOLVED: " <> ToString[result, InputForm]
  ];
  printTree[{
    "condition  " <> text,
    "result     " <> status
  }]
];


SetAttributes[printVariable, HoldAll];
Options[printVariable] = {"Numeric" -> True};

printVariable[sym_Symbol, OptionsPattern[]] := Module[
  {name = SymbolName[Unevaluated[sym]],
  exact = ToString[sym, InputForm]},

  printTree[
    If[OptionValue["Numeric"],
      {
        "exact    " <> name <> " = " <> exact,
        "numeric  " <> name <> " = " <>
          ToString[N[sym, 64], InputForm]
      },
      {
        name <> " = " <> exact
      }
    ]
  ]
];


$printCount = 0;


(* ::Subsection:: *)
(*Inputs*)


(* ::Subsubsection:: *)
(*Configuration*)


(* Input directory *)
dirRoot = FileNameJoin[{DirectoryName[$InputFileName], "Data/1st"}];


(* Number of parallel kernels *)
nproc = 16;


(* ::Subsubsection:: *)
(*Data*)


(* ::Item:: *)
(*vHat (Fourier coefficients of the approximate solution v)*)


vHat = Import[FileNameJoin[{dirRoot, "vHat.m"}]];


Unprotect[MM, NN];
{MM, NN} = Dimensions[vHat];
Protect[MM, NN];


(* ::Item:: *)
(*aCal (the non-trivial part of the isomorphism A)*)


{aCal, setABar} = Import[FileNameJoin[{dirRoot, "aCal.m"}]];


(* ::Item:: *)
(*LambdaTilde (quadratic operator)*)


vSquaredTilde = Import[FileNameJoin[{dirRoot, "vSquaredTilde.m"}]];


{MMM, NNN}=Dimensions[vSquaredTilde]/2;


Mtilde = 2*MMM-1;
Ntilde = 2*NNN-1;


cTildeTab = Table[
  coeffC[m, n, Mtilde, Ntilde, vSquaredTilde, MMM, NNN],
  {m, 0, 2 Mtilde},
  {n, 0, 2 Ntilde}
];


(* ::Subsubsection:: *)
(*Parameters and bounds*)


(* ::Item:: *)
(*Weights*)


Unprotect[\[Rho]t, \[Rho]x];


{\[Rho]t, \[Rho]x} = {1 + 10^-20, 11/10};


Protect[\[Rho]t, \[Rho]x];


(* ::Item:: *)
(*Solution frequency \[CapitalOmega]0=(2p+1)/q*)


Unprotect[\[CapitalOmega]0, p, q];


p = 11;
q = 8;
\[CapitalOmega]0 = (2*p + 1) / q;


Protect[\[CapitalOmega]0, p, q];


(* ::Item:: *)
(*Threshold for H_0P_{m,n} elements*)


Unprotect[sigma];


sigma = 9/10;


Protect[sigma];


(* ::Item:: *)
(*Bounds*)


Unprotect[bndNormV, bndNormN0, bndNormA, bndNormH, bndNormLinv];


bndNormV = 11932577993/2618142364;
bndNormN0 = 1251640/200388335065877;
bndNormA = 19906846310/1894822989;
bndNormH = 5959447941/6621608675;


bndNormLinv = 1720890135138878791416128/1850768595563353538700837;


Protect[bndNormV, bndNormN0, bndNormA, bndNormH, bndNormLinv];


(* ::Item:: *)
(*Constants in the theorem*)


Unprotect[cnstDelta, cnstK];


cnstDelta = 13081330001308133/209065086978900000000000;
cnstK = 6817223577/7573215052;


Protect[cnstDelta, cnstK];


(* ::Subsection:: *)
(*Computations*)


(* ::Subsubsection:: *)
(*Initialisation*)


(* ::Item:: *)
(*System report*)


 printMessage["System report: \n" <>
  "Mathematica version: " <> ToString[$Version] <> "\n" <>
  "Machine name: " <> ToString[$MachineName] <> "\n" <>
  "Machine type: " <> ToString[$MachineType] <> "\n" <>
  "Operating system: " <> ToString[$OperatingSystem] <> "\n" <>
  "System: " <> ToString[$System] <> "\n" <>
  "Processor type: " <> ToString[$ProcessorType] <> "\n" <>
  "Processor cores: " <> ToString[$ProcessorCount] <> "\n" <>
  "System memory (GiB): " <> ToString[N[$SystemMemory/2^30]] <> "\n" <>
  "Script command line: " <> ToString[$ScriptCommandLine, InputForm] <> "\n"
]


(* ::Item:: *)
(*Launching kernels*)


LaunchKernels[nproc];


printVariable[$KernelCount, "Numeric" -> False];


(* ::Subsubsection:: *)
(*Norm of v*)


(* ::Item:: *)
(*Computing norm of v*)


normV = Sum[
  Abs[vHat[[m + 1, n + 1]]] * \[Rho]t^(2 m + 1) * \[Rho]x^(2 n + 1),
  {m, 0, MM - 1},
  {n, 0, NN - 1}
];


printVariable[normV];


(* ::Subsubsection:: *)
(*Residuum N(0)*)


(* ::Item:: *)
(*Generating vSquaredHat (coefficients of v^2)*)


printMessage["Computing vSquaredHat"];


vSquaredHat = coeffsSquareP[vHat, {MM, NN}];//AbsoluteTiming


(* ::Item:: *)
(*Generating vCubedHat (coefficients of v^3)*)


printMessage["Computing vCubedHat"];


vCubedHat = coeffsCubeP[vHat, vSquaredHat, {MM, NN}];//AbsoluteTiming


(* ::Item:: *)
(*Computing norm of the residuum N(0)*)


printMessage["Computing normN0"];


normN0 = Module[{\[Alpha], \[Beta]},
  \[Alpha][mm_, nn_] := -\[CapitalOmega]0^2 (2 mm + 1)^2 + 4 (nn + 1)^2;
  \[Beta][mm_, nn_] := -((4 (2 nn + 1)) / (\[Alpha][mm, nn - 1] * \[Alpha][mm, nn]));
  Total @ Flatten @ ParallelTable[
    Sum[
      Abs[
        vCubedHat[[m + 1, n + 1]] / \[Alpha][m, n] + 
        Sum[
          vCubedHat[[m + 1, j + 1]] * \[Beta][m, j],
          {j, n + 1, 3 NN - 2}
        ] + 
        If[And[m < MM, n < NN], vHat[[m + 1, n + 1]], 0]
      ] * \[Rho]t^(2 m + 1) * \[Rho]x^(2 n + 1),
      {n, 0, 3 NN - 2}
    ],
    {m, 0, 3 MM - 2},
    Method -> "FinestGrained"
  ]
];


printVariable[normN0];


(* ::Subsubsection:: *)
(*Norm of \[CapitalLambda]tilde and the difference \[CapitalLambda]-\[CapitalLambda]tilde*)


(* ::Item:: *)
(*Constructing \[CapitalLambda]*)


Mhat = 2*MM-1;
Nhat = 2*NN-1;


cTab = Table[
  coeffC[m, n, Mhat, Nhat, vSquaredHat, MM, NN],
  {m, 0, 2 Mhat},
  {n, 0, 2 Nhat}
];


(* ::Item:: *)
(*Computing ||\[CapitalLambda]tilde||*)


normLambdaTilde = Sum[
  \[Rho]t^(2 Abs[mtilde - Mtilde]) * \[Rho]x^(2 Abs[ntilde - Ntilde]) * Abs[cTildeTab[[mtilde + 1, ntilde + 1]]],
  {mtilde, 0, 2 Mtilde},
  {ntilde, 0, 2 Ntilde}
];


printVariable[normLambdaTilde];


(* ::Item:: *)
(*Computing ||\[CapitalLambda]tilde-\[CapitalLambda]||*)


normLambdaDiff = With[
  {cTildeTab2 = CenterArray[cTildeTab, {2 Mhat + 1, 2 Nhat + 1}]},
  Sum[
    Abs[cTab[[m + 1, n + 1]] - cTildeTab2[[m + 1, n + 1]]] * \[Rho]t^(2 Abs[m - Mhat]) * \[Rho]x^(2 Abs[n - Nhat]),
    {m, 0, 2 Mhat},
    {n, 0, 2 Nhat}
  ]
];


printVariable[normLambdaDiff];


(* ::Subsubsection:: *)
(*Norm of H0*)


(* ::Item:: *)
(*Determining band \[ScriptCapitalB]0*)


(* ::Subitem:: *)
(*Finding \[Eta] (from Lemma 2, Eq. (26))*)


beta = (sigma / 3) * (normLambdaTilde ^ (-1));


\[Eta]0 = Ceiling[
  Max[
    \[Eta] /. Solve[beta == (4\[Eta](3\[Eta]+1)-((2p+1)/q)^2)/(16(\[Eta]-(2p+1)/(2q))(\[Eta]+1-(2p+1)/(2q))(\[Eta]+(2p+1)/(2q))(\[Eta]+1+(2p+1)/(2q))), \[Eta]]
  ]
];


printVariable[\[Eta]0, "Numeric" -> False];


(* ::Subitem:: *)
(*Finding \[Chi] (Lemma 3, Eq. (28))*)


\[Chi]0 = Ceiling[
  Max[
    \[Chi] /. Solve[beta == 1/((((2p+1)/q)(2\[Chi]+1)-2)(((2p+1)/q)(2\[Chi]+1)+2)), \[Chi]]
  ]
];


printVariable[\[Chi]0, "Numeric" -> False];


(* ::Item:: *)
(*Determining the band \[ScriptCapitalB]*)


\[Eta] = Ceiling[\[Eta]0 + ((2p+1)/q*Mtilde + Ntilde)];
printVariable[\[Eta], "Numeric" -> False];


\[Chi] = Ceiling[\[Chi]0 + (Mtilde+q/(2p+1)*Ntilde)];
printVariable[\[Chi], "Numeric" -> False];


(* ::Item:: *)
(*Defining set \[ScriptCapitalS]_g*)


printMessage["Defining setG"];


setSm1 = pointsInSk[-1, {\[CapitalOmega]0, \[Eta], \[CapitalOmega]0*\[Chi]}, {Mtilde, Ntilde}, {q, 2p+1}];
setS0 = pointsInSk[0, {\[CapitalOmega]0, \[Eta], \[CapitalOmega]0*\[Chi]}, {Mtilde, Ntilde}, {q, 2p+1}];


setG = NestWhileList[
  With[{dataG = ParallelTable[{mn, funG[mn, cTildeTab, {Mtilde, Ntilde}]}, {mn, # + {q, 2p+1}&/@#}, Method->"FinestGrained"]},
    Cases[dataG, {{m_, n_}, g_?(#>=sigma&)} :> {m, n}]
]&, # - {q, 2p+1}&/@setS0, Length[#] > 0&][[2;;-2]];


(* ::Item:: *)
(*Defining set \[ScriptCapitalS]_{H_0}*)


printMessage["Defining setH0"];


setH0 = Flatten[Prepend[setG, setSm1], 1];


(* ::Item:: *)
(*Computing ||H0P_{m,n}||/||P_{m,n}|| on \[ScriptCapitalS]_{H_0}*)


printMessage["Computing dataH0Pmn"];


(* Reverse for performance reasons and better utilisation of the ParallelTable *)
dataH0Pmn = Reverse@ParallelTable[
  {mn, cmpNormH0Pmn[mn, {MMM, NNN}, vSquaredTilde, \[CapitalOmega]0, {\[Rho]t, \[Rho]x}]},
  {mn, Reverse[setH0]},
  Method->"FinestGrained"
];


(* ::Item:: *)
(*Selecting sites which need nontrivial A , set \[ScriptCapitalS]_A*)


printMessage["Defining setA"];


setA = Cases[dataH0Pmn, {{m_, n_}, H0_?(# >= sigma&)} :> {m, n}];


If[Not[And@@(MemberQ[setABar, #]&/@setA)], 
  printMessage["ERROR: setA is not a subset of setABar. Aborting ..."];
  Abort[]
];


(* ::Item:: *)
(*Computing norm of A*)


Kmax = Length[setABar];


normA = Max[{1, 
  Max[
    Table[
      Module[{mQ, nQ, Qbar},
        {mQ, nQ} = setA[[Q]];
        Qbar = Position[setABar, {mQ, nQ}][[1, 1]];
        Sum[
          Module[{mK, nK},
            {mK, nK} = setABar[[K]];
            Abs[aCal[[Qbar, K]]] * \[Rho]t^(2 mK + 1) * \[Rho]x^(2 nK + 1)
          ],
          {K, Kmax}
        ] / (\[Rho]t^(2 mQ + 1) * \[Rho]x^(2 nQ + 1))
      ],
      {Q, Length[setA]}
    ]
  ]}
];


printVariable[normA];


(* ::Item:: *)
(*Application of A on \[ScriptCapitalS]_A*)


mQup = Max[setABar[[All,1]]];
nQup = Max[setABar[[All,2]]];


printMessage["Computing dataZtab"];


dataZtab = ParallelTable[
  PadRight[
    coeffZ[setABar[[Q]], {MMM, NNN}, vSquaredTilde, \[CapitalOmega]0],
    {mQup + 2 MMM, nQup + 2 NNN}
  ],
  {Q, Kmax},
  Method -> "FinestGrained"
];


printMessage["Computing dataHPmn"];


dataHPmn = ParallelTable[
  Module[{mk, nk, K, dataSumZ, mq, nq},
    printMessage[ToString[{"$KernelID" -> $KernelID, "II" -> II}]];
    {mk, nk} = setA[[II]];
    K = Position[setABar, setA[[II]]][[1, 1]];
    dataSumZ = -3 * Sum[aCal[[K, JJ]] * dataZtab[[JJ]], {JJ, Kmax}];
    Do[
      {mq, nq} = setABar[[Q]];
      dataSumZ[[mq + 1, nq + 1]] -= aCal[[K, Q]],
      {Q, Kmax}
    ];
    dataSumZ[[mk + 1, nk + 1]] += 1;
    
    1 / (\[Rho]t^(2 mk + 1) * \[Rho]x^(2 nk + 1)) * 
      Sum[
        Abs[dataSumZ[[\[Mu] + 1, \[Nu] + 1]]] * 
        \[Rho]t^(2 \[Mu] + 1) * 
        \[Rho]x^(2 \[Nu] + 1),
        {\[Mu], 0, mQup + 2 MMM - 1},
        {\[Nu], 0, nQup + 2 NNN - 1}
      ]
  ],
  {II, Length[setA]},
  Method -> "FinestGrained"
];


maxDataHPmn = Max[dataHPmn];


printVariable[maxDataHPmn];


(* ::Item:: *)
(*Computing norm of the 3*L^{-1}(\[CapitalLambda]tilde-\[CapitalLambda])A term*)


errLambda = 3 * bndNormLinv * normA * normLambdaDiff;


printVariable[errLambda];


(* ::Item:: *)
(*Computing ||H||*)


normH = Max[maxDataHPmn, sigma] + errLambda;


(* ::Subsubsection:: *)
(*Verification of bounds*)


printTest[normV <= bndNormV];


printTest[normN0 <= bndNormN0];


printTest[normA <= bndNormA];


printTest[normH <= bndNormH];


(* ::Subsubsection:: *)
(*Proof of the theorem*)


Block[{\[Delta] = cnstDelta, K = cnstK},
  printTest[
    bndNormH + 6 bndNormLinv*bndNormV*bndNormA^2 \[Delta] +
      3 bndNormLinv*bndNormA^3 \[Delta]^2 < K
  ];
  printTest[bndNormN0 < (1 - K)*\[Delta]];
];


(* ::Subsubsection:: *)
(*Finalisation*)


(* ::Item:: *)
(*Reporting memory usage*)


Module[{mainPeak, subkernelPeaks, totalPeak},
  mainPeak = MaxMemoryUsed[];
  subkernelPeaks = ParallelEvaluate[MaxMemoryUsed[]];
  totalPeak = mainPeak + Total[subkernelPeaks];

  printMessage["Memory usage:\n"<>
  "Main kernel peak: "<>ToString[N[mainPeak/2^30]]<>" GiB\n"
  <>"Subkernel peaks: "<>ToString[N[subkernelPeaks/2^30]]<>" GiB\n"
  <>"Sum of kernel peaks: "<>ToString[N[totalPeak/2^30]]<>" GiB"];
];


(* ::Item:: *)
(*Reporting running time*)


$scriptElapsedTime = SessionTime[] - $scriptStartTime;


printMessage["Running time (seconds): " <> ToString[N@$scriptElapsedTime, InputForm] <> "\n"
];


(* ::Item:: *)
(*Closing kernels*)


CloseKernels[];


printMessage["Exiting ..."];
