fl = load("data/fl.mat");
addpath ../zhuang
f_hat = ifft_hermite(fl);
save -ascii "data/f_hat.mat" f_hat;