#!/usr/bin/env wolframscript
Clear["Global`*"]
Needs["Fouriergroups`"]
(*read arguments*)
argv = Echo[Rest @ $ScriptCommandLine,"argv"];
file=Echo[argv[[1]],"file="]
regiontoplot=type/.{"1"->2,"2"->1};(*we are plotting (in black) the polygons whose fourier trafo we do NOT want to take*)
interval=Echo[ToExpression[argv[[2]]],"interval="]
midpoint = {interval + 1, interval + 1};
diffractiondiagramthreshold=Echo[ToExpression[argv[[3]]],"threshold="]
multiples=Echo[ToExpression[argv[[4]]],"multiples"]
gaugethreshold=Echo[ToExpression[argv[[5]]],"gaugethreshold="]
symmetriesfiles=Echo[argv[[6;;]],"symmetriesfiles="]
symmetriestocheck=({StringReplace[#, ".csv" -> ""],Import[#]}&)/@symmetriesfiles
basisvectorsfile="basisvectors.csv"

scriptDirectory = DirectoryName[ExpandFileName[First[$ScriptCommandLine]]];
fftfiles=Echo[(StringJoin[{StringReplace[file, ".png" -> ""],#}]&)/@FFTcsvfilenames[interval],"fftfiles"]
fouriercoefffile=Echo[StringReplace[file, ".png" -> StringJoin[{"fourier_coefficients_",".txt"}]],"fourier coefficients integration"];

image=Import[file];
data=ImageData[image]
Print["Dimensions=",Dimensions[data]];
nrpixel=Dimensions[data](*rows than columns*)
If[Length[Dimensions[data]]>2,Throw["used image is not grayscale"]]

Print["Calculating FFT"];
{FFTcenter,positions,intensities,phases}=calcfft[data,interval,diffractiondiagramthreshold];
exportDataToCSV[StringJoin[{"diffractiondiagram_",StringReplace[file, ".png" -> ""],".csv"}],Transpose[Join[Transpose[positions], {intensities,phases}]],"#k_1,k_2,absolute value of peak,phase of peak"]
Clear[data]
Clear[image]

If[!FileExistsQ[basisvectorsfile],Throw["No basis vectors specified, calculation cannot proceed"]]

basisvectors=Import[basisvectorsfile];
holohedry=generategroup[Import/@symmetriesfiles]
Print["Number of elements in holohedry=",Length[holohedry]]
idcestocalc=getidcestocalculate[FFTcenter,nrpixel,gaugethreshold,holohedry,multiples,basisvectors]

Export[StringJoin[{"Indices_",StringReplace[file, ".png" -> ""],".csv"}],Transpose[Join[Sequence @@ Transpose/@{idcestocalc,idcestocalc.basisvectors}]]]
baseidces=IdentityMatrix[Length[basisvectors]]
otheridces=Complement[idcestocalc,baseidces]
baseidcespos=baseidces.basisvectors
otheridcespos=otheridces.basisvectors

(*calculate exact fourier transform of pixel image*)
Print["Calculating gauges by Fourier transformation of pixels"]
If[FileExistsQ[fouriercoefffile]&& FileDate[fouriercoefffile] >(*basis vectors file has not been modified*) FileDate[basisvectorsfile],
{idxtofourier,idceslefttocalc}=readidcesfromdict[fouriercoefffile,idcestocalc],
(idxtofourier=<||>;idceslefttocalc=idcestocalc)
]
Print["calculating fourier transform for ",Length[idceslefttocalc]," idces"]
pixelintegrationresult=Echo[If[Length[idceslefttocalc]>0,integratenumericallyoverpixel[file,idceslefttocalc.basisvectors],{}],"pixelintres"]
fouriercoeffs=MapThread[(#1*Exp[I*#2]&),{pixelintegrationresult[[All,3]],pixelintegrationresult[[All,4]]}]
idxtofourier=Merge[
	{idxtofourier, AssociationThread[idceslefttocalc,fouriercoeffs]},
	First]
Export[fouriercoefffile,idxtofourier]
baseidcesfourier=(idxtofourier[#]&)/@baseidces
otheridcesfourier=(idxtofourier[#]&)/@otheridces
exportDataToCSV[StringJoin["fourier_coefficients_base_",StringReplace[file, ".png" -> ""],".csv"],Transpose[Join[Transpose[baseidcespos], {Abs/@baseidcesfourier,Arg/@baseidcesfourier}]],"#k_1,k_2,absolute value of peak,phase of peak"]
exportDataToCSV[StringJoin["fourier_coefficients_other_",StringReplace[file, ".png" -> ""],".csv"],Transpose[Join[Transpose[otheridcespos], {Abs/@otheridcesfourier,Arg/@otheridcesfourier}]],"#k_1k_2absolute value of peakphase of peak"]

(*calculate gauges*)
((Print["calculating gaugerrors of ",#[[1]]];writegaugeerrors[#,basisvectors,baseidces,otheridces,idxtofourier,StringJoin["gaugeerrors_",StringReplace[file, ".png" -> ""],"_"]])&)/@symmetriestocheck

