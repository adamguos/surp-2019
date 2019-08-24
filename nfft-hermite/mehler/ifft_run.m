Fl = load("data/Fl.mat");
addpath ../zhuang
Fhat = ifft_hermite(Fl);
save -ascii "data/Fhat.mat" Fhat;