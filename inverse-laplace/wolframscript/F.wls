#!/usr/bin/env wolframscript
(* ::Package:: *)

precision = 100;
interval = {-3, 3};

n = 64;
points = 2^16;
w = (Sqrt[5] - 1) / 2;

A = {1, 3, 1};
T = {-1, 0, 0.75};

F[t_] := (Pi * w)^(-1/2) * Exp[ -(Sqrt[5] / 2) * SetPrecision[t, precision]^2 ] * (Total[Table[ A[[j]] * Exp[ 2 * SetPrecision[t, precision] * T[[j]] ], {j, Dimensions[A][[1]]} ]]);

t = Array[# &, points, interval];
ft = ParallelTable[F[ti], {ti, t}];
Export["temp/t.mat", SetPrecision[t, precision], "List"];
Export["temp/ft.mat", SetPrecision[ft, precision], "List"];
