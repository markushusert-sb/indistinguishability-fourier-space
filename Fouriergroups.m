(* ::Package:: *)
BeginPackage["Fouriergroups`"]
Get["Settings`"]
extinctionQ::usage="extinctionQ[vector,fouriergroup] returns True if the wavevector is an extinction" 
fourierspectrumviagroupaction::usage="fourierspectrumviagroupaction[polygons,n,group] returns fourier point spectrum of polygons under fourier group as list of pairs of wave vectors and coefficients"
fouriercoeff::usage="fouriercoeff[polygon(s),wavevector] returns fouriercoefficient of polygon "
fouriercoeffpixel::usage="returns fourier coefficient of image specified by pixel contained in data. The image is interpreted as having its upper left corner in the origin"
fourierspectrumtocsv::usage="transform a fourier point spectrum to" 
latticepointsfundreg::usage="latticepointsfundreg[n,group] gives integral coefficients of all lattice vectors in the fundamental region of lattice holohedry"
latticepointsfundregphase::usage="latticepointsfundreg[n,group] gives integral coefficients of all lattice vectors in the fundamental region of point group"
latticepoints::usage="latticepoints[n,group] gives integral coefficients of all lattice vectors for given n"
fourierspectrumlegend::usage="legend for csv file containing fourier spectrum"
gauge::usage="calculates gauge function for given transformation"
writegaugeerrors::usage="calculates gauge errors and exports gaugeerrors as function of amplitudes to a file"
generategroup::usage="provide list of generators described as elements of GL(integers,n), returns all elements of generated subgroup of GL(integers,n)"
integratenumericallyoverpixel::usage="calculate fourier transform for a list of wavevectors over a grayscale image"
symmetriesreplacementrule::usage="pattern/.symmetriesreplacementrule yields orthogonal transformations to be investigated for a given pattern"
FFTcsvfilenames::usage="FFTcsvfilename[intervall,pixel] yields filename containing FFT-values for the frequencies in the center square of length 2*intervall+1, calculated using <pixel> pixel in each dimension"
bilinearidx::usage="bilinearidx[array_,idx : vecpattern] indexes into array using floating point index idx by bilinearly interpolating between arrayva&lues around idx"
readidcesfromdict::usage=""
exportDataToCSV::usage="export to csv while prepending header"
relativediff::usage=""
Begin["`Private`"]
readidcesfromdict[dict_String,idces:{vecpattern..}]:=readidcesfromdict[ToExpression[Import[dict]],idces]
readidcesfromdict[dict_Association,idces:{vecpattern..}]:={dict,Select[idces,(!KeyExistsQ[dict,#]&)]}
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
bilinearidx[array_,idx : vecpattern] :=(*see \
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
integratenumericallyoverpixel[imagefile_String,wavevectors:{vecpattern..}]:=(Export["wavevectors.csv",wavevectors];
(*carry out calculation in python, because too expensive in mathematica *)Print[RunProcess[{"integrate_numerically_over_pixel.py",imagefile,"wavevectors.csv"}]];
Import["wavevectors.csv"][[2;;]](*Returns list of 4-tuples with each tuple: (k_1,k_2,abs(fourier(rho)),arg(fourier(rho)))*)
)
intensityerrors[symmetry:{vecpattern ..},otheridces:{vecpattern..},idxtofourier_Association]:=(relativediff[Abs[idxtofourier[#]],Abs[idxtofourier[symmetry.#]]]&)/@otheridces
gaugeerrors[symmetry:{vecpattern ..},baseidces:{vecpattern..},otheridces:{vecpattern..},idxtofourier_Association]:=Module[{gaugeatbase,actualgauges,gaugeerrors},
		gaugeatbase=(gauge[symmetry,#,idxtofourier]&)/@baseidces;
		actualgauges=(gauge[symmetry,#,idxtofourier]&)/@otheridces;
		{Join[gaugeatbase,actualgauges],Join[(0&)(*by definition we have no error for the gauge function at the base*)/@gaugeatbase,MapThread[(diffmod1signed[#1,gaugeatbase.#2(*propagates gauge function from basevectors to considered vector*)]&),{actualgauges,otheridces}]]}
	]
writegaugeerrors[symmetry:{_?StringQ, {vecpattern ..}},basevectors:{vecpattern..},baseidces:{vecpattern..},otheridces:{vecpattern..},idxtofourier_Association,name_String]:=exportDataToCSV[StringJoin[name,StringReplace[symmetry[[1]]," "->"_"],".csv"],Transpose[
Join[Transpose[Join[baseidces,otheridces].basevectors],
{Sequence@@gaugeerrors[symmetry[[2]],baseidces,otheridces,idxtofourier],
(Abs[idxtofourier[#]]&)/@Join[baseidces,otheridces],
intensityerrors[symmetry[[2]],Join[baseidces,otheridces],idxtofourier]}]],"#k_1,k_2,gauge function,gauge error,amplitudes,relative amplitudedifference"]
gauge[symmetry : {vecpattern ..}, idx : vecpattern,idxtofourier_Association] :=
 Module[{newidx, oldval, newval},
	newidx = symmetry.idx;
  oldval = idxtofourier[idx];
  newval = idxtofourier[newidx];
  (*Print["relative difference in amplitudes: ",(Abs[newval]-Abs[oldval])/Abs[oldval]];*)
  Mod[(Arg[newval] - Arg[oldval])/2/Pi, 1]
 ]
phasefunatbasevecs[realspaceele:affineelepattern]:=((-realspaceele[[1]].realspaceele[[2]].# &)/@ {{1,0},{0,1}});
checklatticevecfourier[a :vecpattern, b :vecpattern, name_] :=checklatticevec[a,b,name]
fouriercoeffpixel[pixel:{vecpattern..}, k : {_Integer ..}] :=If[Length[DeleteDuplicates[Dimensions[pixel]]]>1,Throw["pixeldata is not square"],
Total[Flatten[MapIndexed[
	(#1*Exp[-2 Pi*I*(k[[1]]*(#2[[2]]-1)(*2nd index, columns, goes in plus x direction*) - k[[2]]*(#2[[1]]-1)(*1st index, rows, goes in minus y direction*))/Length[pixel]] &), 
    Length[pixel], {2}]]]/Length[pixel]^2]
End[]

EndPackage[]
