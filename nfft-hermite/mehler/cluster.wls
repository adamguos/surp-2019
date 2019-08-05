#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory["/home/adam/Documents/surp-2019/nfft-hermite/mehler"];
precision = 100;

(* PRE-PROCESSING FUNCTIONS *)

(* Compute gausshermite nodes and weights using Julia and import *)
computeGaussHermite[n_] :=
	Module[{nodes, lambda}, 
		Run["julia export_gausshermite_nodes_weights.jl " <> ToString[n]];
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
			W[[i]] = SetPrecision[lambda[[i]] * Exp[ nodes[[i]] ^ 2 ], precision],
		{ i, Dimensions[lambda][[1]] }];
		Return[W]
	]

(* Compute or import Hermite function values for gausshermite nodes *)
hermiteFunc[n_, x_] := SetPrecision[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2], precision]
computePsi[degree_, nodes_] :=
	Module[{psi},
		psi = Parallelize[Table[hermiteFunc[n, x], {n, degree}, {x, nodes}]];
		Export["data/hermite_at_nodes.mat", psi, "Table"];
		Return[psi]
	]
importPsi[] := Import["data/hermite_at_nodes.mat", "Table"]

(* Compute or import Hermite function values for many points on [-n\sqrt{2}, n\sqrt{2}] *)
generateHermiteT[points_, degree_] :=
	Module[{t, psi},
		(* t = SetPrecision[Array[# &, points, {-degree * Sqrt[2], degree * Sqrt[2]}], precision]; *)
		t = SetPrecision[Array[# &, points, {-30, 30}], precision];
		psi = Parallelize[Table[hermiteFunc[n, x], {n, degree}, {x, t}]];
		Export["data/t.mat", t, "List"];
		Export["data/hermite_at_t.mat", psi, "Table"];
		Return[{t, psi}]
	]
importHermiteT[] := {Import["data/t.mat", "List"], Import["data/hermite_at_t.mat", "Table"]}

(* DEFINE CONSTANTS/PARAMETERS *)

n = 128;

A = {3, 2, 1, 4, 5};
T = {-5, -2, 2, 5, 8};
K = Dimensions[A][[1]];

w = (Sqrt[5] - 1) / 2;

numOfPoints = 2048;

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
F[t_] := SetPrecision[(Pi * w)^(-1/2) * Exp[ -(Sqrt[5] / 2) * t^2 ] * Total[Table[ A[[j]] * Exp[ 2 * t * T[[j]] ], {j, 3} ]] , precision];

mu[k_] := SetPrecision[Total[Table[ W[[l]] * F[ nodes[[l]] ] * psi[[k, l]], {l, 2n} ]] * w^(-k) , precision];

H[x_] := Piecewise[{ {1, Abs[x] <= 1/2}, {Exp[-Exp[2 / (1 - 2 * Abs[x])] / (1 - Abs[x])], 1/2 < Abs[x] < 1}, {0, Abs[x] >= 1} }];

(* Tn[xi_] := Total[Table[ H[k/n] * mu[k] * psit[[k, xi]], {k, n} ]]; *)
Tn[xi_] := Total[Table[ H[k/n] * mu[k] * psit[[k, xi]], {k, n/2} ]];

gram[nodes_, degree_] :=
	Module[{psi},
		psi = Table[HermiteH[n, x] / Sqrt[Sqrt[Pi] * 2^n * n!] * Exp[-x^2 / 2], {n, degree / 2}, {x, nodes}];
		Return[Table[Total[Table[W[[j]] * psi[[m, j]] * psi[[l, j]], {j, Dimensions[W][[1]]}]], {m, degree / 2}, {l, degree / 2}]];
	]

(* EXECUTE PROGRAM *)

(* Compute Gram residual *)
G = gram[nodes, n];
Image[G]
Print["Gram residual norm: ", Norm[G - IdentityMatrix[Dimensions[G][[1]]], Infinity]]

(* Plot mu *)
mut = Range[n];
mutn = Table[mu[xi], {xi, mut}];
ListLinePlot[Transpose@{mut, mutn}]

(* Distribute definitions of symbols to kernels for parallelisation *)
DistributeDefinitions[Tn, t, K];

(* Compute Tn *)
time = Now;
tn = ParallelTable[
	Print[ToString[xi-1], " out of ", ToString[Dimensions[t][[1]]] ". Time elapsed: ", ToString[Now - time]];
	time = Now;
	Tn[xi],
{xi, Dimensions[t][[1]]}];

(* Export t and tn *)
dirName = "output/" <> DateString[{"Year", "Month", "Day", "-", "Hour", "Minute", "Second"}];
CreateDirectory[dirName];
Export[dirName <> "/t.mat", t, "List"];
Export[dirName <> "/tn.mat", tn, "List"];

(* Import Tn *)
(* tn = Import["output/tn.mat"] *)

ListLinePlot[Transpose@{t, tn}]