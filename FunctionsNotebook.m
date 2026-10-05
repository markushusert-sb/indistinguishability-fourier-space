BeginPackage["FunctionsNotebook`"]
Get["Fouriergroups`"]
genpointg::usage=""
integratenumericallyoverpixel::usage=""
gaugeerrors::usage=""
plotsymmetryerrors::usage=""
combinederrors::usage=""
amplitudeerrors::usage=""
plotbasevectors::usage=""
setbasevectors::usage=""
getidcestotest::usage=""
findlocalmaximum::usage=""
fouriercoeffpixel::usage=""
calcfft::usage=""
plotfft::usage=""
plotfftone::usage=""
plotselectedvectors::usage=""
Begin["`Private`"]
gencosets[gens_] := Module[{donecosets, lastcosets, donetrafo},
        lastcosets = donecosets = {identitytrafo[Length[gens[[1,1]]]]};
        donetrafo = {donecosets[[1]][[2]]};
        While[Length[lastcosets] > 0,
        donecosets = 
      Join[donecosets, genlayercosets[gens, lastcosets, donetrafo]];
        ];
        donecosets
  ]
genpointg[gens_] := gencosets[({0*Range[Length[#]],#}&)/@ gens][[All,2]]
identitytrafo[n_:2]:={(0&)/@Range[n], IdentityMatrix[n]}
cosetop[x_, y_] := {x[[1]] + x[[2]] . y[[1]],
  x[[2]] . y[[2]]} (*coset x multiplied with coset y*)
genlayercosets[gens_, lastcosets_, donetrafo_] :=
 (*creates new layer of potential cosets, by right and left multiplying gens to lastcosets and und updating lastcosets and donetrafo after*)
 Module[{potentcoset, cosets2add},
  cosets2add = {};
  Do [(
    potentcoset = cosetop[gens[[i]], lastcosets[[j]]];
    If [! MemberQ[donetrafo, potentcoset[[2]]],
     AppendTo[cosets2add, potentcoset];
     AppendTo[donetrafo, potentcoset[[2]]];
     ];
    potentcoset = cosetop[lastcosets[[j]], gens[[i]]];
    If [! MemberQ[donetrafo, potentcoset[[2]]],
     AppendTo[cosets2add, potentcoset];
     AppendTo[donetrafo, potentcoset[[2]]];
     ]
    )
   , {i, 1, Length[gens]}, {j, 1, Length[lastcosets]}];
  lastcosets = cosets2add
  ]
SetAttributes[genlayercosets, HoldAll]
findlocalmaximum[imagefile_String,vector:_?VectorQ]:=(
	rasterpoints= (# + vector) & /@ Tuples[Subdivide[-0.25, 0.25, 16], 2];
	amplitudes=integratenumericallyoverpixel[imagefile,rasterpoints][[;;,3]];
	maxPoint=rasterpoints[[First@Ordering[amplitudes, -1]]];
	Print["local maximum at ",maxPoint];
	{maxPoint,ListContourPlot[
		MapThread[Append, {rasterpoints, amplitudes}],
		Contours -> 20,
		Epilog -> {
				Red,
				PointSize[0.02],
				Point[maxPoint],
				Black,
				Text[
					Style["Local maximum", 10, Bold],
					maxPoint,
					{0, -1.5}
				]
			},
		PlotLegends ->BarLegend[Automatic]
	]}
	)
Nextint[x_] := Floor[x] + 1
bilinearidx[array_,idx : _?VectorQ] :=(*see \
https://en.wikipedia.org/wiki/Bilinear_interpolation*)
(
	 {Nextint[idx[[1]]] - idx[[1]],idx[[1]] -Floor[idx[[1]]]} . (
	{{array[[Floor[idx[[1]]], Floor[idx[[2]]]]],array[[Floor[idx[[1]]], Nextint[idx[[2]]]]]},
{array[[Nextint[idx[[1]]], Floor[idx[[2]]]]],array[[Nextint[idx[[1]]], Nextint[idx[[2]]]]]}}.
{Nextint[idx[[2]]] - idx[[2]], idx[[2]] - Floor[idx[[2]]]})
	 )
setbasevectors[basevectors:_?MatrixQ,numberoffundfreqs:_?IntegerQ,intervallsize:_?IntegerQ]:=(
Column[Function[{vecidx},
   Row[(
      Column[{
         Subscript["k",
          StringJoin[ToString[vecidx],
           ",", # /. {1 -> "x", 2 -> "y"}]],

         Manipulator[
          Dynamic[basevectors[[vecidx, #]]], {-intervallsize/2, intervallsize/2},
          AppearanceElements -> {"InputField"}, Appearance -> "Open"]
         }] &) /@ Range[2]]] /@ Range[numberoffundfreqs]]
)
SetAttributes[setbasevectors, HoldAll]
plotbasevectors[basevectors_?MatrixQ,numberoffundfreqs_?IntegerQ]:=(
GraphicsRow[{
	ListPlot[plotfftsettings[],
		Epilog -> {
     Join @@ (({{Red, PointSize[Large], Point[basevectors[[#]]]},
           Text[Style[Subscript["k", "i"], Red, 14, Bold],
            basevectors[[#]], {-1, -1.5}]} &) /@
        Range[numberoffundfreqs])
     }
	],
	diffractionlegend[FFTamplitudes,sizes]		
	}]
)
calcfft[pixels:_?MatrixQ,diffractiondiagramthreshold_,intervallsize_]:=(Module[{FFTcenter,FFTcoords,FFTamplitudes},
	nrpixel = Reverse[Dimensions[pixels]];(*reversing so that first component goes in x-direction*)
	FFTcenter=(Conjugate[(*conjugate because we have negative exponent in fourier trafo, but mathematica a positive one*)
	Transpose[Reverse[(*reverse and transpose, since in the indixes, the first element is the rows (vertical downwards) and the second element is the cols (horizontal rightward), while for Listplot, it is the x coordinates (horizontal rightwards) and the y coordinates (vertical upwards)*)RotateLeft[
	Fourier[pixels], {-Floor[(nrpixel[[1]] - 1)/2], -Floor[nrpixel[[2]]/2]}]]]]/Sqrt[Times @@ nrpixel])[[Floor[nrpixel[[2]]/2] + 1 - interval ;; Floor[nrpixel[[2]]/2] + 1 + interval, Floor[nrpixel[[1]]/2] + 1 - interval ;; Floor[nrpixel[[1]]/2] + 1 + interval]];
FFTfiltered =
  SparseArray[
   Map[(If[Abs[#] > diffractiondiagramthreshold, #, 0] &),
    FFTcenter, {2}]];
FFTfiltered[[intervallsize/2+1,intervallsize/2+1]]=0;
(*Extract coordinates and values*)
rules = Most@
    ArrayRules[Abs[FFTfiltered]];  (*drop default zero rule*);
idces = (# - {intervallsize/2 + 1, intervallsize/2 + 1} &) /@ Keys[rules];
FFTcoords =N[Transpose[{{Sqrt[Times@@nrpixel]/nrpixel[[1]],0},{0,Sqrt[Times@@nrpixel]/nrpixel[[2]]}}.Transpose[idces]]];
(* N[Transpose[{{Sqrt[Times @@ nrpixel]/nrpixel[[1]], 0}, {0,
      Sqrt[Times @@ nrpixel]/nrpixel[[2]]}} . Transpose[idces]]];*)
FFTamplitudes = Values[rules];
{FFTcenter,FFTcoords,FFTamplitudes}
])
(* Define a circular marker graphic *)
hollowCircle[size_]:= Graphics[
   {
     EdgeForm[Black],
     FaceForm[None],
     Disk[{0,0},Offset[size]]
   },
	ImageSize->size*2.4
]
integratenumericallyoverpixel[imagefile_String,wavevectors:_?MatrixQ]:=Module[{pythonresult},Export["wavevectors.csv",wavevectors];
(*carry out calculation in python, because too expensive in mathematica *)RunProcess[{"python","integrate_image_numerically.py",imagefile,"wavevectors.csv"}];
pythonresult=Import["wavevectors.csv"][[2;;]](*Returns list of 4-tuples with each tuple: (k_1,k_2,abs(fourier(rho)),arg(fourier(rho)))*);
pythonresult[[All,3;;4]]
]
plotsymmetryerrors[fouriercoeffs:_?VectorQ,idcestotest:{{_?IntegerQ..}..},symmetries:{{_?StringQ,_?MatrixQ}..}]:=Module[{idxtofourier,baseidces,otheridces},(*fouriercoeffs: list of pairs of amplitude and phase*)
idxtofourier=AssociationThread[idcestotest,fouriercoeffs];
baseidces=IdentityMatrix[Length[idcestotest[[1]]]];
otheridces=Complement[idcestotest,baseidces];
Labeled[GraphicsRow[
	Function[{symmetry},Labeled[Histogram[combinederrors[symmetry[[2]],baseidces,otheridces,idxtofourier]],symmetry[[1]],Bottom]]/@symmetries
],Superscript["\[CapitalDelta]","sym"],Top]
]
combinederrors[symmetry:_?MatrixQ,baseidces:{{_?IntegerQ..}..},otheridces:{{_?IntegerQ..}..},idxtofourier_Association]:=(1-(1-2*Abs[gaugeerrors[symmetry,baseidces,otheridces,idxtofourier][[2,All]]])*(1-Abs[intensityerrors[symmetry,Join[baseidces,otheridces],idxtofourier]]))

diffractionlegend[FFTamplitudes:_?VectorQ,sizes:_?VectorQ]:=GraphicsColumn[
	(Labeled[hollowCircle[#[[1]]],Style[StringJoin[{"|\[Rho]|=",ToString[#[[2]]]}], FontSize -> 18],Right]&)/@{{Max[sizes]/3,Max[FFTamplitudes]/3},{2*Max[sizes]/3,2*Max[FFTamplitudes]/3},{Max[sizes],Max[FFTamplitudes]}},
	Spacings -> 5, Method -> {"ShrinkWrap" -> True},ImageSize->Small,Frame->True
	]
plotfftsettings[FFTcoords:_?MatrixQ,sizes:_?VectorQ]:=Sequence[
		({#}&)/@FFTcoords,
		PlotStyle -> Black,
		PlotMarkers -> ( hollowCircle[#] &) /@ sizes ,
		Frame -> True,
		AxesLabel -> {Subscript["k", "x"], Subscript["k", "y"]},
		PlotLabel -> "Diffraction Diagram",
		AspectRatio -> Automatic,
		GridLines -> Automatic
]
plotfft[FFTcoords:_?MatrixQ,FFTamplitudes:_?VectorQ]:=(
	sizes = FFTamplitudes/(Max[FFTamplitudes])*10;
	GraphicsRow[{ListPlot[plotfftsettings[FFTcoords,sizes]
		],
		diffractionlegend[FFTamplitudes,sizes]		
	}]
)
plotselectedvectors[FFTcoords:_?MatrixQ,FFTamplitudes:_?VectorQ,consideredfrequencies:_?MatrixQ]:=(
sizes = FFTamplitudes/(Max[FFTamplitudes])*10;
GraphicsRow[{
	ListPlot[plotfftsettings[FFTcoords,sizes],
	 Epilog -> {
		 ({{Blue, PointSize[Large], Point[#]}} &) /@ consideredfrequencies
		 },
	 AspectRatio -> "Automatic", GridLines -> "Automatic",
	 AxesLabel -> {Subscript["k", "x"], Subscript["k", "y"]},
	 PlotLabel -> "Diffraction Diagram"],
	 GraphicsColumn[{
			diffractionlegend[FFTamplitudes, sizes],
			GraphicsColumn[{Labeled[
	 			Graphics[{FaceForm[Blue], Disk[{0, 0}, Offset[3]]},ImageSize -> 8], 
		 		Style[TraditionalForm[Symbol["\[ScriptCapitalM]"]^Style["s", FontSlant -> "Plain"]], FontSize -> 20],
		 Right]},Frame->True,Spacings->0,Method -> {"ShrinkWrap" -> True},ImageSize->Small]
	}]
}]
)
End[]
EndPackage[]
