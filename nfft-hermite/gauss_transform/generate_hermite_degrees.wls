#!/usr/bin/env wolframscript
(* ::Package:: *)

Now


samples = 512


degrees = Array[# &, 16, {32, 512}]


intervalFactor = 1.05


h[n_, x_] := N[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2]]


Xr = ConstantArray[0, {Dimensions[degrees][[1]], samples}]
Do[
	n = degrees[[i]];
	X = RandomReal[{N[-Sqrt[2n + 1]] * intervalFactor, N[Sqrt[2n + 1]] * intervalFactor}, {samples}];
	Xr[[i]] = X
, {i, Dimensions[degrees][[1]]}]


Export["Xr.mat", Xr, "Table"]


Do[
	H = Parallelize[Table[h[n, x], {n, degrees[[i]]}, {x, Xr[[i]]}]];
	Export["H" <> IntegerString[degrees[[i]], 10, 3] <> ".mat", H, "Table"]
, {i, Dimensions[degrees][[1]]}]


Now
