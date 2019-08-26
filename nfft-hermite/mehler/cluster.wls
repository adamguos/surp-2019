#!/usr/bin/env wolframscript
(* ::Package:: *)

(* DEFINE CONSTANTS/PARAMETERS *)

(* SetDirectory["~/Documents/surp-2019/nfft-hermite/mehler"] *)
precision = 200;
interval = {-2, 2};
(* interval = {-degree * Sqrt[2], degree * Sqrt[2]}; *)

n = 80;
numOfPoints = 2048;
w = (Sqrt[5] - 1) / 2;

A = {1, 2, 1};
T = {-1, 0, 0.75};

(* err[t_] := RandomReal[{-0.0005, 0.0005}]; *)
err[t_] := 0



(* PRE-PROCESSING FUNCTIONS *)

(* Compute gausshermite nodes and weights using Julia and import *)
computeGaussHermite[n_] :=
	Module[{nodes, lambda}, 
		Run["/usr/local/bin/julia export_gausshermite_nodes_weights.jl " <> ToString[n]];
		nodes = Import["data/nodes.mat", "List"];
		lambda = Import["data/weights.mat", "List"];
		Return[{nodes, lambda}]
	]
importGaussHermite[] := {Import["data/nodes.mat", "List"], Import["data/weights.mat", "List"]}

(* Multiply gausshermite weights by factor of exp(node_i^2) *)
computeWeights[nodes_, lambda_] :=
	Module[{W},
		W = ConstantArray[0, { Dimensions[lambda][[1]] }];
		Do[
			W[[i]] = lambda[[i]] * Exp[ nodes[[i]] ^ 2 ],
		{ i, Dimensions[lambda][[1]] }];
		Return[W]
	]

(* Compute or import Hermite function values for gausshermite nodes *)
hermiteFunc[n_, x_] := HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2]
computePsi[degree_, nodes_] :=
	Module[{psi},
		psi = ParallelTable[hermiteFunc[n, x], {n, degree}, {x, nodes}];
		Export["data/hermite_at_nodes.mat", SetPrecision[psi, precision], "Table"];
		Return[psi]
	]
importPsi[] := Import["data/hermite_at_nodes.mat", "Table"]

(* Compute or import Hermite function values for many points on [-n\sqrt{2}, n\sqrt{2}] *)
generateHermiteT[points_, degree_] :=
	Module[{t, psi},
		t = Array[# &, points, interval];
		psi = ParallelTable[hermiteFunc[n, x], {n, degree}, {x, t}];
		Export["data/t.mat", SetPrecision[t, precision], "List"];
		Export["data/hermite_at_t.mat", SetPrecision[psi, precision], "Table"];
		Return[{t, psi}]
	]
importHermiteT[] := {Import["data/t.mat", "List"], Import["data/hermite_at_t.mat", "Table"]}



(* COMPUTE/IMPORT DATA *)

(* Re-compute data *)
Print["Computing Gauss Hermite nodes and weights..."]; {nodes, lambda} = computeGaussHermite[2n];
W = computeWeights[nodes, lambda];
Print["Computing Hermite function values at nodes..."]; psi = computePsi[n, nodes];
Print["Computing Hermite function values at many points t..."]; {t, psit} = generateHermiteT[numOfPoints, n];

(* Import previously-computed data *)
(* Print["Importing Gauss Hermite nodes and weights..."]; {nodes, lambda} = importGaussHermite[];
W = computeWeights[nodes, lambda];
Print["Importing Hermite function values at nodes..."]; psi = importPsi[];
Print["Importing Hermite function values at many points t..."]; {t, psit} = importHermiteT[]; *)



(* DEFINE COMPUTATIONAL FUNCTIONS *)

(* functions depend on constants/parameters defined above *)

F[t_] := (Pi * w)^(-1/2) * Exp[ -(Sqrt[5] / 2) * SetPrecision[t, precision]^2 ] * (Total[Table[ A[[j]] * Exp[ 2 * SetPrecision[t, precision] * T[[j]] ], {j, Dimensions[A][[1]]} ]]);

mu[k_] := Total[Table[ W[[l]] * F[ nodes[[l]] ] * psi[[k, l]], {l, 2n} ]] * w^(-k);

H[x_] := Piecewise[{ {1, Abs[x] <= 1/2}, {Exp[-Exp[2 / (1 - 2 * Abs[x])] / (1 - Abs[x])], 1/2 < Abs[x] < 1}, {0, Abs[x] >= 1} }];

Tn[xi_] := Total[Table[ H[k/n] * mu[k] * psit[[k, xi]], {k, n} ]];

gram[nodes_, degree_] :=
	Module[{psi},
		psi = Table[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2], {n, degree / 2}, {x, nodes}];
		Return[Table[Total[Table[W[[j]] * psi[[m, j]] * psi[[l, j]], {j, Dimensions[W][[1]]}]], {m, degree / 2}, {l, degree / 2}]];
	]

(* create tables of values computed by ifft_hermite *)

Fl = ParallelTable[F[t], {t, nodes}];
Export["data/fl.mat", SetPrecision[Fl, precision], "List"];
Run["octave-cli ifft_run.m"];	(* compute inverse FFT from sampled values of F (at zeros of Hermite polynomials) *)

Fhat = Import["data/f_hat.mat", "List"];
mul = Table[Fhat[[k]] / w^(k-1), {k, Dimensions[Fhat][[1]]}];
Tnl[xi_] := Total[Table[ H[k/n] * mul[[k]] * psit[[k, xi]], {k, n} ]];

(* test plot mu *)
mut = Range[n];
mutn = Table[mu[t], {t, mut}];
Export["temp/mut.mat", SetPrecision[mut, precision], "List"];
Export["temp/mutn.mat", SetPrecision[mutn, precision], "List"];
(* mut = Range[Dimensions[mul][[1]]];
Export["temp/mut.mat", SetPrecision[mut, precision], "List"];
Export["temp/mutn.mat", SetPrecision[mul, precision], "List"]; *)



(* EXECUTE PROGRAM *)

(* Compute Gram residual *)
G = gram[nodes, n];
Image[G]
Print["Gram residual norm: ", Norm[G - IdentityMatrix[Dimensions[G][[1]]], Infinity]]

(* Compute Tn *)
time = Now;
tn = ConstantArray[0, { Dimensions[t][[1]] }]
tn = ParallelTable[
	Print[ToString[xi], " out of ", ToString[Dimensions[t][[1]]] ". Time elapsed: ", ToString[Now - time]];
	time = Now;
	Tn[xi],
{xi, Dimensions[t][[1]]}];

(* Compute Tnl using ifft *)
time = Now;
tnl = ParallelTable[
	Print[ToString[xi], " out of ", ToString[Dimensions[t][[1]]] ". Time elapsed: ", ToString[Now - time]];
	time = Now;
	Tnl[xi],
{xi, Dimensions[t][[1]]}];

(* Export t and tn *)
dirName = "output/" <> DateString[{"Year", "Month", "Day", "-", "Hour", "Minute", "Second"}];
CreateDirectory[dirName]
SetDirectory[dirName];
Export["t.mat", SetPrecision[t, precision], "List"];
Export["tn.mat", SetPrecision[tn, precision], "List"];
Export["tnl.mat", SetPrecision[tnl, precision], "List"];
Export["mut.mat", SetPrecision[mut, precision], "List"];
Export["mutn.mat", SetPrecision[mutn, precision], "List"];
Export["A.mat", SetPrecision[A, precision], "List"];
Export["T.mat", SetPrecision[T, precision], "List"];