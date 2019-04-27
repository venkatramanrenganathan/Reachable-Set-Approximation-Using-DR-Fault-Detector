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
attack_cov   = dummy'*dummy;
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



%% Compute worst case noise and worst case attack Reachable Set

% Get Worst Case Attack
attack_details     = worst_attack_distribution(attack_input_param);
attack_probability = attack_details.probabilities;
attack_support     = attack_details.support;

% Get Worst Case Noise
noise_details     = worst_noise_distribution(noise_input_param);
noise_probability = noise_details.probabilities;
noise_support     = noise_details.support;

for j=1:M
    x_noise  = zeros(n,1);
    e_noise  = zeros(n,1);
    x_attack = zeros(n,1);
    e_attack = zeros(n,1);
    for i=1:N
        attack_input(:,i) = attack_support(:,sample(attack_probability));
        system_noise(:,i) = noise_support(:,sample(noise_probability));
    end    
    for i = 1:N    
        index                  = (j-1)*N + i;
        x_noise                = (A + B_original*G)* x_noise - B_original*G*e_noise + system_noise(:,i);
        x_attack               = (A + B_original*G)* x_attack - B_original*G*e_attack;
        e_noise                = A*e_noise + system_noise(:,i);
        e_attack               = A*e_attack - L*sqrtm(residual_var)*attack_input(:,i);
        x_plot_noise_worst(:,index)  = x_noise;
        x_plot_attack_worst(:,index) = x_attack;    
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
        
        attack_input           = 2*rand(m,1)-1;
        attack_input           = sqrt(z_k)*attack_input/norm(attack_input);
        x_noise                = (A + B_original*G)* x_noise - B_original*G*e_noise + system_noise;
        x_attack               = (A + B_original*G)* x_attack - B_original*G*e_attack;
        e_noise                = A*e_noise + system_noise;
        e_attack               = A*e_attack - L*sqrtm(residual_var)*attack_input;
        x_plot_noise_Gaussian(:,index) = x_noise;
        x_plot_attack_uniform(:,index) = x_attack;    
    end        
end

%% perform the geometric sum between noise and attacke reachable sets
sum_index = index*1000;

for i=1:sum_index
    j  =  floor(rand*index)+1;
    k  =  floor(rand*index)+1;    
    jj =  floor(rand*index)+1;
    kk =  floor(rand*index)+1;
    x_plot_total_Gaussian(:,i) = x_plot_noise_Gaussian(:,j) + x_plot_attack_uniform(:,k);
    x_plot_total_worst(:,i)    = x_plot_noise_worst(:,jj) + x_plot_attack_worst(:,kk);
end

%% Calculate the Distributionally Robust Ellipsoid for Noise
% Call the function that will give the DR Ellipsoid
noise_output_param    = cvx_noise_bounding_ellipsoid(noise_input_param);
newnoise_output_param = cvx_gaussiannoise_bounding_ellipsoid(noise_input_param);

% Unbox the function output parameters
P_noise                       = noise_output_param.P_noise;
PP_noise                      = newnoise_output_param.PP_noise;
min_volume_of_noise_ellipsoid = noise_output_param.min_volume;


%% Calculate the Distributionally Robust Ellipsoid for Attack
% Call the function that will give the DR Ellipsoid
attack_output_param    = cvx_attack_bounding_ellipsoid(attack_input_param);
newattack_output_param = cvx_gaussianattack_bounding_ellipsoid(attack_input_param);

% Unbox the function output parameters
P_attack                       = attack_output_param.P_attack;
PP_attack                      = newattack_output_param.PP_attack;
min_volume_of_attack_ellipsoid = attack_output_param.min_volume;


%% Compute the Total Bounding Ellipsoid as Minkowski Sum of both ellipsoids
Q_noise   = inv(P_noise);
QQ_noise  = inv(PP_noise);
Q_attack  = inv(P_attack);
QQ_attack = inv(PP_attack);
Q_total   = (sqrt(trace(Q_noise)) + sqrt(trace(Q_attack)))*((Q_noise/sqrt(trace(Q_noise))) +...
            (Q_attack/ sqrt(trace(Q_attack))));
QQ_total  = (sqrt(trace(QQ_noise)) + sqrt(trace(QQ_attack)))*((QQ_noise/sqrt(trace(QQ_noise))) +...
            (QQ_attack/ sqrt(trace(QQ_attack))));
P_total   = inv(Q_total);
PP_total  = inv(QQ_total);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             1                           %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the Reachable States - uniform ATTACK
figure;
plot(x_plot_attack_uniform(1,:),x_plot_attack_uniform(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Attack Ellipsoid 
syms xx yy 
ellipse_function = [xx yy]*PP_attack*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'blue','lineWidth',4);
% leg = legend('Reachable States','DR - Attack Ellipsoid','Gaussian - Attack Ellipsoid', 'Interpreter','latex');
% set(leg, 'Interpreter', 'latex')
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
title('Uniform distribution attack','FontSize', 24)
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             2                           %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the Reachable States - Gaussian NOISE
figure;
plot(x_plot_noise_Gaussian(1,:),x_plot_noise_Gaussian(2,:),'.k');
grid on
hold on;

% Plot the Distributionally Robust Noise Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*PP_noise*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'red','lineWidth',4);
axis equal;
title('Gaussian Noise','FontSize', 24)
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             3                          %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the Total reachable states and bounding ellipsoid
% Plot the Reachable States - TOTAL
figure;
plot(x_plot_total_Gaussian(1,:),x_plot_total_Gaussian(2,:),'.k');
grid on
hold on;

% Plot the Distributionally Robust Total Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*PP_total*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'red','lineWidth',4);
axis equal;
title('Gaussian Noise with Uniform Attack','FontSize', 24)
hold on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             4                           %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the Reachable States - worst-NOISE
figure;
plot(x_plot_noise_worst(1,:),x_plot_noise_worst(2,:),'.k');
grid on
hold on;
% Plot the worst Noise Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_noise*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'blue','lineWidth',4);
% leg2 = legend('Reachable States','DR - Noise Ellipsoid','Gaussian - Noise Ellipsoid');
% set(leg2, 'Interpreter', 'latex')
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
title('Worst Case Noise','FontSize', 24)
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             5                           %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the Reachable States - worst ATTACK
figure;
plot(x_plot_attack_worst(1,:),x_plot_attack_worst(2,:),'.k');
grid on
hold on;
% Plot the Distributionally Robust Attack Ellipsoid 
syms xx yy 
ellipse_function  = [xx yy]*P_attack*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'red','lineWidth',4);
axis equal;
title('Worst Case attack','FontSize', 24)
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             6                           %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the Reachable States - TOTAL
figure;
plot(x_plot_total_worst(1,:),x_plot_total_worst(2,:),'.k');
grid on
hold on;
% Plot the Gaussian Total Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*P_total*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'blue','lineWidth',4);
% leg3 = legend('Reachable States','DR - Total Ellipsoid','Gaussian - Total Ellipsoid');
% set(leg3, 'Interpreter', 'latex')
axis equal;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
title('Worst Case Noise with Worst Case Attack','FontSize', 24)
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off
