% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 2nd April, 2019.

clear all; clc; close all;

%% Problem Data

% System Matrices
N  = 100;
M  = 200;
A  = [0.84  0.23
      -0.47 0.12];
B  = [0.07 -0.32
      0.23 0.58];
C  = [1 0
      2 1];
K  = [1.404 -1.042
      1.842 1.008];
L  = [0.0276   0.0448
      -0.01998 -0.0290];
n  = size(A,1);
m  = size(L,1);
Sigma_w      = [0.045  -0.011               
                -0.011 0.02];
mu_noise     = zeros(n,1);
residual_var = [2.086 0.134
                0.134 2.230];               
cov_attack   = residual_var;                
mu_attack    = [0.0 0.0]';                  
alarm_rate   = 0.05; 
B_original   = B;




A_delta = B*K*L*sqrtm(residual_var);
alfa    = 40;


cvx_begin
    variable attack_input(m,1)      
    maximize -trace((B*K*L*sqrtm(residual_var)*attack_input)'*(B*K*L*sqrtm(residual_var)*attack_input))
    subject to
        attack_input'*attack_input <= alfa;
cvx_end


