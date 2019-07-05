clear; close all; clc
alpha_vec = logspace(-3,1,20);
tau_vec = logspace(-3,-0.01,20);
regularizationTerm = @(v,tau) 0.0;
N = 100;
dt = 1e-3;
Nblock = 400;
filename = 'filename.mat';
integrateSystem(alpha_vec,tau_vec,N,filename,Nblock,dt,regularizationTerm)