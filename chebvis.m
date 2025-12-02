clear all; close all; clc; 

Nz = 10;
N = Nz;
%N = Nz + 1;

k_indices = N:-1:0;

cheb_11 = cos(k_indices * pi / N);
y = zeros(1,N+1);

plot(cheb_11, y,'o');