#!/usr/bin/env wolframscript
Clear["Global`*"]
Needs["Fouriergroups`"]
(*read arguments*)
argv = Echo[Rest @ $ScriptCommandLine,"argv"];
file=Echo[argv[[1]],"file="]
regiontoplot=type/.{"1"->2,"2"->1};(*we are plotting (in black) the polygons whose fourier trafo we do NOT want to take*)
intervall=Echo[ToExpression[argv[[2]]],"intervall="]
midpoint = {intervall + 1, intervall + 1};
diffractiondiagramthreshold=Echo[ToExpression[argv[[3]]],"threshold="]
multiples=Echo[ToExpression[argv[[4]]],"multiples"]
gaugethreshold=Echo[ToExpression[argv[[5]]],"gaugethreshold="]
symmetriesfiles=Echo[argv[[6;;]],"symmetriesfiles="]
symmetriestocheck=({StringReplace[#, ".csv" -> ""],Import[#]}&)/@symmetriesfiles
basevectorsfile="basevectors.csv"
holohedrygeneratorsfiles={"generator_holohedry_1.csv","generator_holohedry_2.csv"}

scriptDirectory = DirectoryName[ExpandFileName[First[$ScriptCommandLine]]];
fftfiles=Echo[(StringJoin[{StringReplace[file, ".png" -> ""],#}]&)/@FFTcsvfilenames[intervall],"fftfiles"]
fouriercoefffile=Echo[StringReplace[file, ".png" -> StringJoin[{"fourier_coefficients_",".txt"}]],"fourier coefficients integration"];

image=Import[file];
data=ImageData[image]
Print["Dimensions=",Dimensions[data]];
nrpixel=Dimensions[data]
If[Length[Dimensions[data]]>2,Throw["used image is not grayscale"]]

Print["Calculating FFT"];
FFTcenter=(Conjugate[(*conjugate because we have negative exponent in fourier trafo, but mathematica a positive one*)
Transpose[Reverse[(*reverse and transpose, since in the indixes, the first element is the rows (vertical downwards) and the second element is the cols (horizontal rightward), while for Listplot, it is the x coordinates (horizontal rightwards) and the y coordinates (vertical upwards)*)RotateLeft[
Fourier[data], {-Floor[nrpixel[[2]]/2], -Floor[nrpixel[[1]]/2]}]]]]/Sqrt[Times @@ nrpixel])[[Floor[nrpixel[[1]]/2] + 1 - intervall ;; Floor[nrpixel[[1]]/2] + 1 + intervall, Floor[nrpixel[[2]]/2] + 1 - intervall ;; Floor[nrpixel[[2]]/2] + 1 + intervall]]
Clear[data]
Clear[image]
Export[StringJoin[{"FFT_heatmap_",StringReplace[file, ".png" -> ""],".csv"}],Abs[FFTcenter]]
(*print diffraction image*)
FFTfiltered=SparseArray[Map[(If[Abs[#]>diffractiondiagramthreshold,#,0]&),FFTcenter,{2}]];
peakpositionaspixels=FFTfiltered["NonzeroPositions"];
Print["peakpositionaspixels=",(#-midpoint&)/@peakpositionaspixels];
peakpositionaspixelsstretched=N[Transpose[{{Sqrt[Times@@nrpixel]/nrpixel[[1]],0},{0,Sqrt[Times@@nrpixel]/nrpixel[[2]]}}.Transpose[(#-midpoint&)/@peakpositionaspixels]]];(*apply a stretch to the wave vectors in case of a rectangular image, the so that the basis of the fourier space remains orthonormal*)
peakintensities=(Abs[FFTfiltered[[#[[1]], #[[2]]]]] &) /@peakpositionaspixels;
peakphases=(Arg[FFTfiltered[[#[[1]], #[[2]]]]] &) /@peakpositionaspixels;
exportDataToCSV[StringJoin[{"diffractiondiagram_",StringReplace[file, ".png" -> ""],".csv"}],Transpose[Join[Transpose[peakpositionaspixelsstretched], {peakintensities,peakphases}]],"#k_1,k_2,absolute value of peak,phase of peak"]

(*Determine indices for which to verify gauge*)
(*Do not modify existing sets of basevectors, just add a new counter*)
basevectors=Import[basevectorsfile];
holohedry=generategroup[Import/@holohedrygeneratorsfiles]
Print["Number of elements in holohedry=",Length[holohedry]]
idcestotestorbit=(DeleteDuplicates[Function[{trafo}, trafo . #] /@ holohedry] &) /@ 
 Select[
	Tuples[Range[-multiples, multiples], Length[basevectors]],
	(0<Total[Abs[#]] <= multiples &)]
idcestotest=Flatten[Select[idcestotestorbit,(AllTrue[#,Function[{idx},Norm[#.basevectors]<intervall]]&)],1] (*remove long range frequencies beyond intervall*)
idxtoabsfourierfft=AssociationThread[idcestotest,(bilinearidx[Abs[FFTcenter], midpoint + ({{Sqrt[Times@@nrpixel]/nrpixel[[2]],0},{0,Sqrt[Times@@nrpixel]/nrpixel[[1]]}}.(#.basevectors))]&)/@idcestotest](*estimate amplitude of fourier coefficient using FFT*)
idcestotestexisting=Select[idcestotest,(idxtoabsfourierfft[#]>gaugethreshold&)](*only those with non-negligeble amplitude*)
idcestocalc=Echo[DeleteDuplicates[Flatten[(Function[{trafo}, trafo . #] /@ holohedry &) /@idcestotestexisting,1]],"idces where fourier transform is needed"](*orbit of all existing frequencies*)

Export[StringJoin[{"Indices_",StringReplace[file, ".png" -> ""],".csv"}],Transpose[Join[Sequence @@ Transpose/@{idcestocalc,idcestocalc.basevectors}]]]
baseidces=IdentityMatrix[Length[basevectors]]
otheridces=Complement[idcestocalc,baseidces]
baseidcespos=baseidces.basevectors
otheridcespos=otheridces.basevectors

(*calculate exact fourier transform of pixel image*)
Print["Calculating gauges by Fourier transformation of pixels"]
If[FileExistsQ[fouriercoefffile]&& FileDate[fouriercoefffile] >(*basis vectors file has not been modified*) FileDate[basevectorsfile],
{idxtofourier,idceslefttocalc}=readidcesfromdict[fouriercoefffile,idcestocalc],
(idxtofourier=<||>;idceslefttocalc=idcestocalc)
]
Print["calculating fourier transform for ",Length[idceslefttocalc]," idces"]
pixelintegrationresult=Echo[If[Length[idceslefttocalc]>0,integratenumericallyoverpixel[file,idceslefttocalc.basevectors],{}],"pixelintres"]
fouriercoeffs=MapThread[(#1*Exp[I*#2]&),{pixelintegrationresult[[All,3]],pixelintegrationresult[[All,4]]}]
idxtofourier=Merge[
	{idxtofourier, AssociationThread[idceslefttocalc,fouriercoeffs]},
	First]
Export[fouriercoefffile,idxtofourier]
baseidcesfourier=(idxtofourier[#]&)/@baseidces
otheridcesfourier=(idxtofourier[#]&)/@otheridces
exportDataToCSV[StringJoin["wavevectors_base_",StringReplace[file, ".png" -> ""],".csv"],Transpose[Join[Transpose[baseidcespos], {Abs/@baseidcesfourier,Arg/@baseidcesfourier}]],"#k_1,k_2,absolute value of peak,phase of peak"]
exportDataToCSV[StringJoin["wavevectors_other_",StringReplace[file, ".png" -> ""],".csv"],Transpose[Join[Transpose[otheridcespos], {Abs/@otheridcesfourier,Arg/@otheridcesfourier}]],"#k_1k_2absolute value of peakphase of peak"]

(*calculate gauges*)
((Print["calculating gaugerrors of ",#[[1]]];writegaugeerrors[#,basevectors,baseidces,otheridces,idxtofourier,StringJoin["gaugeerrors_",StringReplace[file, ".png" -> ""],"_"]])&)/@symmetriestocheck

