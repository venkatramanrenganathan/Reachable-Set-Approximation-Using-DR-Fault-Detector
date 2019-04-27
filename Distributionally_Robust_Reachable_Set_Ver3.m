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
attack_cov   = residual_var;                % dummy'*dummy;
mu_attack    = zeros(n,1);                  % Zero Mean for Attack Input
alarm_rate   = 0.05; 
B_original   = B;

% Box the function input parameters for attack
attack_input_param.A            = A;
attack_input_param.B            = B;
attack_input_param.L            = L;
attack_input_param.G            = G;
attack_input_param.mu_attack    = mu_attack;
attack_input_param.attack_cov   = residual_var;
attack_input_param.residual_var = residual_var;
attack_input_param.alarm_rate   = alarm_rate;

%% Compute worst case attack Reachable Set
attack_details     = worst_attack_distribution(attack_input_param);
attack_probability = attack_details.probabilities;
attack_support     = attack_details.support;
attack_threshold   = attack_details.attack_threshold;

for i = 1:N
    attack_input(:,i) = attack_support(:,sample(attack_probability));      
end
desired_attack_cov = cov(attack_input');

for j=1:M
    x_attack = zeros(n,1);
    e_attack = zeros(n,1);
    for i=1:N   
        index                   = (j-1)*N + i;        
        x_attack                = (A + B_original*G)* x_attack - B_original*G*e_attack;        
        e_attack                = A*e_attack - L*sqrtm(residual_var)*attack_input(:,i);        
        x_worst_attack(:,index) = x_attack;    
    end        
end

% % Worst Case Design - Approach - 2
% M = inv(eye(n)-A-B*G)*B*G*inv(eye(n)-A)*L*sqrtm(residual_var);
% eig_val_MM           = eig(M'*M);  
% [max_eig_MM,max_idx] = max(eig_val_MM);
% [min_eig_MM,min_idx] = min(eig_val_MM);
% [eig_vec_MM, val]    = eig(M'*M);  
% max_1_eig_vec_MM     = eig_vec_MM(:,max_idx);
% max_2_eig_vec_MM     = eig_vec_MM(:,min_idx);
% 
% 
% v_max_matrix = max_1_eig_vec_MM*max_1_eig_vec_MM';
% v_min_matrix = max_2_eig_vec_MM*max_2_eig_vec_MM';
% eig_ratio    = max_eig_MM/min_eig_MM;
% 
% 
% % Solve SDP
% cvx_begin sdp quiet    
%     variable gamma_dot
%     variable beta_dot
%     maximize   (gamma_dot)
%     subject to
%         attack_cov == gamma_dot*v_max_matrix + eig_ratio^2*v_min_matrix*(beta_dot - gamma_dot);
%         gamma_dot  <= beta_dot;
% cvx_end

%% Compute Uniform Attack Reachable Set

% Get Attack threshold based on inverse gamma function as chi-squared detector is used
chi_squared_attack_threshold = ncx2inv(1-alarm_rate,size(C,1),0);
for j = 1:M            
    x_attack = zeros(n,1);
    e_attack = zeros(n,1);    
    for i = 1:N            
        index  = (j-1)*N + i;           
        x_rand = rand;        
        if (x_rand < 1-alarm_rate)
            z_k = chi_squared_attack_threshold*x_rand/(1-alarm_rate);
        else
            z_k = 10*chi_squared_attack_threshold + 2*chi_squared_attack_threshold*...
                  (x_rand - 1 + alarm_rate)/(1 - alarm_rate);
        end        
        attack_input              = 2*rand(m,1)-1;
        attack_input              = sqrt(z_k)*attack_input/norm(attack_input);        
        e_attack                  = A*e_attack - L*sqrtm(residual_var)*attack_input;
        x_attack                  = (A + B_original*G)*x_attack - B_original*G*e_attack;                
        x_uniform_attack(:,index) = x_attack;    
    end        
end

%% Calculate the Distributionally Robust Ellipsoid for Attack
attack_input_param.attack_threshold = attack_threshold;
% Type = 1 - DR Case
attack_input_param.type = 1; 
DR_attack_output_param  = attack_bounding_ellipsoid(attack_input_param);

% Type = 2 - Chi Squared Case
attack_input_param.type      = 2; 
Gaussian_attack_output_param = attack_bounding_ellipsoid(attack_input_param);

% Unbox the function output parameters
DR_Ellipse          = DR_attack_output_param.P;
Chi_Squared_Ellipse = Gaussian_attack_output_param.P;


%% Plotting Code

% Case 1 Uniform Attack
figure;
% Plot the Chi Squared Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*Chi_Squared_Ellipse*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
hold on;
% Plot the DR Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*DR_Ellipse*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'b','lineWidth',4);
axis equal;
plot(x_uniform_attack(1,:),x_uniform_attack(2,:),'.k');
grid on
hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
legend('Chi-Squared Ellipsoid','DR Ellipsoid', 'Reachable States', 'Interpreter', 'latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off


% Case 2 Worst Attack
figure;
% Plot the Chi Squared Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*Chi_Squared_Ellipse*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'k','lineWidth',4);
hold on;
% Plot the DR Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*DR_Ellipse*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
fimplicit(ellipse_function_handle,[-7 7 -7 7],'b','lineWidth',4);
axis equal;
plot(x_worst_attack(1,:),x_worst_attack(2,:),'.k');
grid on
hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
legend('Chi-Squared Ellipsoid','DR Ellipsoid', 'Reachable States', 'Interpreter', 'latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off