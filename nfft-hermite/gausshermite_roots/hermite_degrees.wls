#!/usr/bin/env wolframscript
(* ::Package:: *)

Now


SetDirectory["/home/adam/Documents/surp-2019/nfft-hermite/gausshermite_roots/data/nodes"]
F = FileNames[]
SetDirectory["/home/adam/Documents/surp-2019/nfft-hermite/gausshermite_roots/data"]


h[n_, x_] := SetPrecision[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2], 50]


psi[n_, x_] := SetPrecision[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!], 50]


Do[
	X = Import[Directory[] <> "/nodes/" <> f, "List"];
	degree = Dimensions[X][[1]];
	H = Parallelize[Table[psi[n, x], {n, degree}, {x, X}]];
	Export[Directory[] <> "/hermites/H_" <> IntegerString[degree, 10, 4] <> ".mat", H, "Table"];
	Print[degree],
{f, F}]


Now
