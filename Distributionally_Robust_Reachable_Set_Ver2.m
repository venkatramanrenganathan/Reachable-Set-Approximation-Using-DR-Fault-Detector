% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 27th February, 2019.

clear all; clc; close all;

%% Problem Data

% System Matrices
A  = [0.84  0.23
      -0.47 0.12];
B  = [0.07 -0.32
      0.23 0.58];
C  = [1 0
      2 1];
G  = [1.404 -1.042
      1.842 1.008];
L  = [0.0276   0.0448
      -0.01998 -0.0290];
n  = size(A,1);
m  = size(L,1);
N = 100;
M = 200;
Sigma_w      = [0.045  -0.011               % System Noise Covariance Matrix
                -0.011 0.02];
mu_noise     = zeros(n,1);
residual_var = [2.086 0.134
                0.134 2.230];               % Residual covariance
dummy        = rand(n,n);
attack_cov   = 10*Sigma_w; % dummy'*dummy;
mu_attack    = zeros(n,1);                  % Zero Mean for Attack Input
alarm_rate   = 0.05; 
B_original   = B;

% Box the function input parameters for attack
attack_input_param.A            = A;
attack_input_param.B            = B;
attack_input_param.L            = L;
attack_input_param.G            = G;
attack_input_param.mu_attack    = mu_attack;
attack_input_param.attack_cov   = attack_cov;
attack_input_param.residual_var = residual_var;
attack_input_param.alarm_rate   = alarm_rate;

% Box the function input parameters
noise_input_param.A          = A;
noise_input_param.B          = B;
noise_input_param.mu_noise   = mu_noise;
noise_input_param.Sigma_w    = Sigma_w;
noise_input_param.alarm_rate = alarm_rate;

%%% it isnt right to define a threshold for noise because in the absence of the attack 
%%% the relation between z and noise isnt the same with the relation
%%% between z and attack input in presence of attack. but any way the most
%%% critical bounding ellipsoid for the noise is corresponding to a noise
%%% inside an ellipsoid design by its covariance that has a threshold and worst case distribution based
%%% on its known mean and covariance.



%% Compute worst case attack Reachable Set

% Get Worst Case Attack
attack_details     = worst_attack_distribution(attack_input_param);
attack_probability = attack_details.probabilities;
attack_support     = attack_details.support;

for j=1:M
    x_attack = zeros(n,1);
    e_attack = zeros(n,1);
    for i=1:N
        attack_input(:,i)       = attack_support(:,sample(attack_probability));               
        index                   = (j-1)*N + i;        
        x_attack                = (A + B_original*G)* x_attack - B_original*G*e_attack;        
        e_attack                = A*e_attack - L*sqrtm(residual_var)*attack_input(:,i);        
        x_worst_attack(:,index) = x_attack;    
    end        
end
%% Compute Gaussian Noise and Uniform Attack Reachable Set

% With Gaussian Noise, we can find threshold based on inverse gamma
% function as chi-squared detector is used
attack_threshold  = ncx2inv(1-alarm_rate,size(C,1),0);

for j = 1:M    
    
    x_noise  = zeros(n,1);
    e_noise  = zeros(n,1);
    x_attack = zeros(n,1);
    e_attack = zeros(n,1);
    
    for i = 1:N    
        
        index        = (j-1)*N + i;   
        system_noise = mvnrnd(mu_noise,Sigma_w,1)'; % System Noise
        x_rand       = rand;
        
        if (x_rand < 1-alarm_rate)
            z_k = attack_threshold*x_rand/(1-alarm_rate);
        else
            z_k = 10*attack_threshold + 2*attack_threshold*(x_rand - 1 + alarm_rate)/(1 - alarm_rate);
        end
        
        attack_input = 2*rand(m,1)-1;
        attack_input = sqrt(z_k)*attack_input/norm(attack_input);
        x_noise      = (A + B_original*G)* x_noise - B_original*G*e_noise + system_noise;
        x_attack     = (A + B_original*G)* x_attack - B_original*G*e_attack;
        e_noise      = A*e_noise + system_noise;
        e_attack     = A*e_attack - L*sqrtm(residual_var)*attack_input;
        
        x_gaussian_noise(:,index) = x_noise;
        x_uniform_attack(:,index) = x_attack;    
    end        
end

%% perform the geometric sum between noise and attacke reachable sets
sum_index = index*1000;

for i=1:sum_index
    j   =  floor(rand*index)+1;
    k   =  floor(rand*index)+1;    
    jj  =  floor(rand*index)+1;
    kk  =  floor(rand*index)+1;
    x_gaussian_noise_uniform_attack(:,i) = x_gaussian_noise(:,j) + x_uniform_attack(:,k);    
    x_gaussian_noise_worst_attack(:,i)   = x_gaussian_noise(:,jj) + x_worst_attack(:,kk);
end

%% Calculate the Distributionally Robust Ellipsoid for Noise
% Call the function that will give the DR Ellipsoid
% Type = 1 - DR Case
noise_input_param.type = 1; 
DR_noise_output_param  = noise_bounding_ellipsoid(noise_input_param);

% Type = 2 - Gaussian Case
noise_input_param.type      = 2; 
Gaussian_noise_output_param = noise_bounding_ellipsoid(noise_input_param);

% Unbox the function output parameters
P_Ellipse_DR_Noise       = DR_noise_output_param.P;
P_Ellipse_Gaussian_Noise = Gaussian_noise_output_param.P;


%% Calculate the Distributionally Robust Ellipsoid for Attack
% Call the function that will give the DR Ellipsoid
% Type = 1 - DR Case
attack_input_param.type = 1; 
DR_attack_output_param  = attack_bounding_ellipsoid(attack_input_param);

% Type = 2 - Gaussian Case
attack_input_param.type      = 2; 
Gaussian_attack_output_param = attack_bounding_ellipsoid(attack_input_param);

% Unbox the function output parameters
DR_attack_ellipse       = DR_attack_output_param.P;
Gaussian_attack_ellipse = Gaussian_attack_output_param.P;


%% Plotting Cases
% Case 1: Gaussian Noise - Chi-Squared Ellipse, Worst Attack - Chi-Squared Attack Ellipse
P_case_1 = compute_minkovsky_sum(P_Ellipse_Gaussian_Noise, Gaussian_attack_ellipse);

% Case 2: Gaussian Noise - Chi-Squared Ellipse, Worst Attack - DR Attack Ellipse
P_case_2 = compute_minkovsky_sum(P_Ellipse_Gaussian_Noise, DR_attack_ellipse);

% Case 3: Gaussian Noise - DR Ellipse, Worst Attack - Chi-Squared Attack Ellipse
P_case_3 = compute_minkovsky_sum(P_Ellipse_DR_Noise, Gaussian_attack_ellipse);

% Case 4: Gaussian Noise - DR Ellipse, Worst Attack - DR Attack Ellipse
P_case_4 = compute_minkovsky_sum(P_Ellipse_DR_Noise, DR_attack_ellipse);

% Case 5: Gaussian Noise - Chi-Squared Ellipse, Uniform Attack - Chi-Squared Attack Ellipse
P_case_5 = compute_minkovsky_sum(P_Ellipse_Gaussian_Noise, Gaussian_attack_ellipse);

% Case 6: Gaussian Noise - Chi-Squared Ellipse, Uniform Attack - DR Attack Ellipse
P_case_6 = compute_minkovsky_sum(P_Ellipse_Gaussian_Noise, DR_attack_ellipse);

% Case 7: Gaussian Noise - DR Ellipse, Uniform Attack - Chi-Squared Attack Ellipse
P_case_7 = compute_minkovsky_sum(P_Ellipse_DR_Noise, Gaussian_attack_ellipse);

% Case 8: Gaussian Noise - DR Ellipse, Uniform Attack - DR Attack Ellipse
P_case_8 = compute_minkovsky_sum(P_Ellipse_DR_Noise, DR_attack_ellipse);


% Case 1
figure;
plot(x_gaussian_noise_uniform_attack(1,:),x_gaussian_noise_uniform_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_1*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off

% Case 2
figure;
plot(x_gaussian_noise_uniform_attack(1,:),x_gaussian_noise_uniform_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_2*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off

% Case 3
figure;
plot(x_gaussian_noise_uniform_attack(1,:),x_gaussian_noise_uniform_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_3*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off


% Case 4
figure;
plot(x_gaussian_noise_uniform_attack(1,:),x_gaussian_noise_uniform_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_4*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off


% Case 5
figure;
plot(x_gaussian_noise_worst_attack(1,:),x_gaussian_noise_worst_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_5*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off


% Case 6
figure;
plot(x_gaussian_noise_worst_attack(1,:),x_gaussian_noise_worst_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_6*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off

% Case 7
figure;
plot(x_gaussian_noise_worst_attack(1,:),x_gaussian_noise_worst_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_7*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off

% Case 8
figure;
plot(x_gaussian_noise_worst_attack(1,:),x_gaussian_noise_worst_attack(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_case_8*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h = fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
uistack(h(1), 'top')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off
