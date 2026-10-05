(* ::Package:: *)
BeginPackage["Fouriergroups`"]
extinctionQ::usage="extinctionQ[vector,fouriergroup] returns True if the wavevector is an extinction" 
diffmod1::usage="distance (unsigned) of 2 numbers modulo 1"
diffmod1signed::usage="difference (signed) of 2 numbers modulo 1"
fourierspectrumviagroupaction::usage="fourierspectrumviagroupaction[polygons,n,group] returns fourier point spectrum of polygons under fourier group as list of pairs of wave vectors and coefficients"
fouriercoeff::usage="fouriercoeff[polygon(s),wavevector] returns fouriercoefficient of polygon "
fouriercoeffpixel::usage="returns fourier coefficient of image specified by pixel contained in data. The image is interpreted as having its upper left corner in the origin"
fourierspectrumtocsv::usage="transform a fourier point spectrum to" 
latticepointsfundreg::usage="latticepointsfundreg[n,group] gives integral coefficients of all lattice vectors in the fundamental region of lattice holohedry"
latticepointsfundregphase::usage="latticepointsfundreg[n,group] gives integral coefficients of all lattice vectors in the fundamental region of point group"
latticepoints::usage="latticepoints[n,group] gives integral coefficients of all lattice vectors for given n"
fourierspectrumlegend::usage="legend for csv file containing fourier spectrum"
gauge::usage="calculates gauge function for given transformation"
intensityerrors::usage="intensityerrors[symmetry,idces,idces_to_fourier] calculates first deviation measure for given symmetry"
gaugeerrors::usage="gaugeerrors[symmetry,idces_base,idces_other,idces_to_fourier] returns gauge function and second deviation measure for given symmmetry"
writegaugeerrors::usage="calculates gauge errors and exports gaugeerrors as function of amplitudes to a file"
generategroup::usage="provide list of generators described as elements of GL(integers,n), returns all elements of generated subgroup of GL(integers,n)"
integratenumericallyoverpixel::usage="calculate fourier transform for a list of wavevectors over a grayscale image"
symmetriesreplacementrule::usage="pattern/.symmetriesreplacementrule yields orthogonal transformations to be investigated for a given pattern"
FFTcsvfilenames::usage="FFTcsvfilename[intervall,pixel] yields filename containing FFT-values for the frequencies in the center square of length 2*intervall+1, calculated using <pixel> pixel in each dimension"
bilinearidx::usage="bilinearidx[array_,idx : _?VectorQ] indexes into array using floating point index idx by bilinearly interpolating between arrayva&lues around idx"
readidcesfromdict::usage=""
exportDataToCSV::usage="export to csv while prepending header"
relativediff::usage=""
getidcestocalculate::usage"get idces of studied frequencies with regards to basis vectors"
calcfft::usage"calcfft[data,interval,threshold] calculates FFT in data for frequencies contained in interval and filters by amplitude using threshold"
Begin["`Private`"]
calcfft[data:_?MatrixQ,interval:_?IntegerQ,threshold:_?NumericQ]:=Module[{FFTcenter,FFTfiltered,peakpositionaspixels,positions,intensities,phases},
	nrpixel=Dimensions[data];(*rows than columns*)
	midpoint = {interval + 1, interval + 1};
	FFTcenter=(Conjugate[(*conjugate because we have negative exponent in fourier trafo, but mathematica a positive one*)
	Transpose[Reverse[(*reverse and transpose, since in the indixes, the first element is the rows (vertical downwards) and the second element is the cols (horizontal rightward), while for Listplot, it is the x coordinates (horizontal rightwards) and the y coordinates (vertical upwards)*)
		RotateLeft[
			Fourier[data],
			{-Floor[(nrpixel[[1]] - 1)/2], -Floor[nrpixel[[2]]/2]}]]]
]/Sqrt[Times @@ nrpixel])
	[[Floor[nrpixel[[2]]/2] + 1 - interval ;; Floor[nrpixel[[2]]/2] + 1 + interval, Floor[nrpixel[[1]]/2] + 1 - interval ;; Floor[nrpixel[[1]]/2] + 1 + interval]];
	FFTfiltered=SparseArray[Map[(If[Abs[#]>threshold,#,0]&),FFTcenter,{2}]];
	peakpositionaspixels=FFTfiltered["NonzeroPositions"];
	positions=N[Transpose[{{Sqrt[Times@@nrpixel]/nrpixel[[1]],0},{0,Sqrt[Times@@nrpixel]/nrpixel[[2]]}}.Transpose[(#-midpoint&)/@peakpositionaspixels]]];(*apply a stretch to the wave vectors in case of a rectangular image, the so that the basis of the fourier space remains orthonormal*)
	intensities=(Abs[FFTfiltered[[#[[1]], #[[2]]]]] &) /@peakpositionaspixels;
	phases=(Arg[FFTfiltered[[#[[1]], #[[2]]]]] &) /@peakpositionaspixels;
	{FFTcenter,positions,intensities,phases}
]
getidcestocalculate[FFTcenter:_?MatrixQ,nrpixel:{_Integer,_Integer},gaugethreshold:_?NumericQ, holohedry:{_?MatrixQ..}, multiples:_?IntegerQ, basisvectors:_?MatrixQ]:=Module[{idcestotestorbit,midpoint,idcestotest,idxtoabsfourierfft,idcestotestexisting,interval},
	interval=(Length[FFTcenter]-1)/2;
	midpoint = {interval + 1, interval + 1};
	idcestotestorbit=(DeleteDuplicates[Function[{trafo}, trafo . #] /@ holohedry] &) /@ 
	 Select[
		Tuples[Range[-multiples, multiples], Length[basisvectors]],
		(0<Total[Abs[#]] <= multiples &)];
	idcestotest=Flatten[Select[idcestotestorbit,(AllTrue[#,Function[{idx},Norm[#.basisvectors]<interval]]&)],1]; (*remove long range frequencies beyond interval*)
	idxtoabsfourierfft=AssociationThread[idcestotest,(bilinearidx[Abs[FFTcenter], midpoint + ({{Sqrt[Times@@nrpixel]/nrpixel[[2]],0},{0,Sqrt[Times@@nrpixel]/nrpixel[[1]]}}.(#.basisvectors))]&)/@idcestotest];(*estimate amplitude of fourier coefficient using FFT*)
	idcestotestexisting=Select[idcestotest,(idxtoabsfourierfft[#]>gaugethreshold&)];(*only those with non-negligeble amplitude*)
	Echo[DeleteDuplicates[Flatten[(Function[{trafo}, trafo . #] /@ holohedry &) /@idcestotestexisting,1]],"idces where fourier transform is needed"](*orbit of all existing frequencies*)
]
readidcesfromdict[dict_String,idces:{_?VectorQ..}]:=readidcesfromdict[ToExpression[Import[dict]],idces]
readidcesfromdict[dict_Association,idces:{_?VectorQ..}]:={dict,Select[idces,(!KeyExistsQ[dict,#]&)]}
exportDataToCSV[fileName_, data_, header_String] := Module[
  {tempFile, stream, csvContent},

  (* Export the data without the header *)
  Export[fileName, data, "CSV"];

  (* Read the CSV content as plain text *)
  csvContent = Import[fileName, "Text"];

  (* Create a temporary file and write the header *)
  tempFile = CreateTemporary[];
  stream = OpenWrite[tempFile];
  WriteLine[stream, header];

  (* Write the original CSV content after the header *)
  WriteString[stream, csvContent];

  (* Close the stream *)
  Close[stream];

  (* Overwrite the original file with the updated content *)
  CopyFile[tempFile, fileName, OverwriteTarget -> True];
  DeleteFile[tempFile];
]
relativediff[x_?NumericQ,y_?NumericQ]:=(x-y)/Max[x,y]

Nextint[x_] := Floor[x] + 1
FFTcsvfilenames[intervall_Integer]:={StringJoin["intervall_",ToString[intervall],"real.csv"],StringJoin["intervall_",ToString[intervall],"imaginary.csv"]}
bilinearidx[array_,idx :_?VectorQ] :=(*see \
https://en.wikipedia.org/wiki/Bilinear_interpolation*)
(
	 {Nextint[idx[[1]]] - idx[[1]],idx[[1]] -Floor[idx[[1]]]} . (
	{{array[[Floor[idx[[1]]], Floor[idx[[2]]]]],array[[Floor[idx[[1]]], Nextint[idx[[2]]]]]},
{array[[Nextint[idx[[1]]], Floor[idx[[2]]]]],array[[Nextint[idx[[1]]], Nextint[idx[[2]]]]]}}.
{Nextint[idx[[2]]] - idx[[2]], idx[[2]] - Floor[idx[[2]]]})
	 )
generategroup[generators:{{{__Integer} ..}...}]:=NestWhile[
		(DeleteDuplicates[Join[#, Flatten[Outer[Dot, #, generators, 1], 1]]] &),
		{IdentityMatrix[Length[generators[[1,1]]]]}
		, Length[#2] > Length[#1] &, 2]
integratenumericallyoverpixel[imagefile_String,wavevectors:{_?VectorQ..}]:=(Export["wavevectors.csv",wavevectors];
(*carry out calculation in python, because too expensive in mathematica *)Print[RunProcess[{"integrate_numerically_over_pixel.py",imagefile,"wavevectors.csv"}]];
Import["wavevectors.csv"][[2;;]](*Returns list of 4-tuples with each tuple: (k_1,k_2,abs(fourier(rho)),arg(fourier(rho)))*)
)
intensityerrors[symmetry:{_?VectorQ ..},otheridces:{_?VectorQ..},idxtofourier_Association]:=(relativediff[Abs[idxtofourier[#]],Abs[idxtofourier[symmetry.#]]]&)/@otheridces
gaugeerrors[symmetry:{_?VectorQ ..},baseidces:{_?VectorQ..},otheridces:{_?VectorQ..},idxtofourier_Association]:=Module[{gaugeatbase,actualgauges,gaugeerrors},
		gaugeatbase=(gauge[symmetry,#,idxtofourier]&)/@baseidces;
		actualgauges=(gauge[symmetry,#,idxtofourier]&)/@otheridces;
		{Join[gaugeatbase,actualgauges],Join[(0&)(*by definition we have no error for the gauge function at the base*)/@gaugeatbase,MapThread[(diffmod1signed[#1,gaugeatbase.#2(*propagates gauge function from basevectors to considered vector*)]&),{actualgauges,otheridces}]]}
	]
diffmod1signed[x_?NumericQ, y_?NumericQ] :=MinimalBy[{Mod[x - y, 1], -Mod[y - x, 1]}, Abs][[1]]
diffmod1[x_?NumericQ,
  y_?NumericQ] := (((*Print["difference of ",x," and ",y,
     " is ",#];*)#) &)[
  Min[1 - Abs[Mod[x, 1] - Mod[y, 1]], Abs[Mod[x, 1] - Mod[y, 1]]]]

writegaugeerrors[symmetry:{_?StringQ, {_?VectorQ ..}},basevectors:{_?VectorQ..},baseidces:{_?VectorQ..},otheridces:{_?VectorQ..},idxtofourier_Association,name_String]:=exportDataToCSV[StringJoin[name,StringReplace[symmetry[[1]]," "->"_"],".csv"],Transpose[
Join[Transpose[Join[baseidces,otheridces].basevectors],
{Sequence@@gaugeerrors[symmetry[[2]],baseidces,otheridces,idxtofourier],
(Abs[idxtofourier[#]]&)/@Join[baseidces,otheridces],
intensityerrors[symmetry[[2]],Join[baseidces,otheridces],idxtofourier]}]],"#k_1,k_2,gauge function,gauge error,amplitudes,relative amplitudedifference"]
gauge[symmetry : _?MatrixQ, idx : _?VectorQ,idxtofourier_Association] :=
 Module[{newidx, oldval, newval},
	newidx = symmetry.idx;
  oldval = idxtofourier[idx];
  newval = idxtofourier[newidx];
  (*Print["relative difference in amplitudes: ",(Abs[newval]-Abs[oldval])/Abs[oldval]];*)
  Mod[(Arg[newval] - Arg[oldval])/2/Pi, 1]
 ]
phasefunatbasevecs[realspaceele:affineelepattern]:=((-realspaceele[[1]].realspaceele[[2]].# &)/@ {{1,0},{0,1}});
End[]

EndPackage[]
