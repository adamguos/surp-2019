#!/usr/bin/env wolframscript
(* ::Package:: *)

Now


degree = 192


samples = 512


trials = 1


intervalFactor = 1.05


interval = N[Sqrt[2 * degree + 1]] * intervalFactor
-interval


X = ConstantArray[0, {trials, samples}]
For[i = 1, i <= Dimensions[X][[1]], i++, X[[i]] = RandomReal[{-interval, interval}, {Dimensions[X][[2]]}]]


Export["X" <> IntegerString[degree, 10, 3] <> ".mat", X, "Table"]


h[n_, x_] := N[HermiteH[n, x] / N[Sqrt[Sqrt[Pi] * 2^n * n!], 300] * Exp[-x^2 / 2], 300]


Do[
	H = Parallelize[Table[N[h[n, x]], {n, degree}, {x, X[[i]]}]];
	Export["H" <> IntegerString[degree, 10, 3] <> "_" <> IntegerString[i, 10, 4] <> ".mat", H, "Table"],
{i, Dimensions[X][[1]]}]


Now
