#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Generates Hermite function values at equally spaced points *)
(* Expects 2 arguments: degree (n), number of points *)
(* Writes 2 files: hermite_t.mat, t.mat *)

start = Now

degree = ToExpression[$ScriptCommandLine[[2]]]
points = ToExpression[$ScriptCommandLine[[3]]]
low = -degree * Sqrt[2]
up = degree * Sqrt[2]

X = SetPrecision[Array[# &, points, {low, up}], 50]

(* Arbitrarily set precision to 50 digits *)
h[n_, x_] := SetPrecision[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2], 50]	(* Hermite function *)
psi[n_, x_] := SetPrecision[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!], 50]	(* Hermite polynomial *)

H = Parallelize[Table[h[n, x], {n, degree}, {x, X}]]
Export["data/hermite_t.mat", H, "Table"]
Export["data/t.mat", X, "List"]

Print["Time taken: " <> TextString[Now - start]]