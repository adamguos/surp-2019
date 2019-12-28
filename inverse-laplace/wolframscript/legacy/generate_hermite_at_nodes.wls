#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Generates Hermite function values at given nodes *)
(* Expects 2 arguments: nodesDir, hermitesDir *)

start = Now

nodesDir = $ScriptCommandLine[[2]]
hermitesDir = $ScriptCommandLine[[3]]

F = FileNames[All, nodesDir]

(* Arbitrarily set precision to 50 digits *)
h[n_, x_] := SetPrecision[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2], 50]	(* Hermite function *)
psi[n_, x_] := SetPrecision[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!], 50]	(* Hermite polynomial *)

Do[
	X = Import[f, "List"];
	degree = Dimensions[X][[1]];
	H = Parallelize[Table[h[n, x], {n, degree}, {x, X}]];
	Export[hermitesDir <> "/H_" <> IntegerString[degree, 10, 4] <> ".mat", H, "Table"];
	Print["n = " <> ToString[degree] <> " saved"],
{f, F}]

Print["Time taken: " <> TextString[Now - start]]