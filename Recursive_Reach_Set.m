% Recursive Reachable Set Estimation
% 25 April 2019
% Venkatraman Renganathan

%% System Data

clear all; close all; clc;

dt = 1;
A  = [1 dt 0 0
      0 1  0 0 
      0 0  1 dt
      0 0  0 1];
B  = [dt^2/2 0
      dt     0
      0      dt^2/2
      0      dt];
C  = [1 0 0 0
      0 0 1 0];
Bo = [0.5 1]';
Bc = [0 0 0 0]';

%% Simulation Data
N    = 50;
h    = 3;
alfa = 13.28;
n    = size(A,1); % number of states
m    = size(B,2); % number of inputs
p    = size(C,1); % number of outputs
Q    = eye(n);
R    = eye(m);
P    = eye(n);

%%  Data Structures
X_o     = zeros(n,n,N);

% Compute predicted error covariance as solution of DARE
Sigma_e = dare(A',C',Q,R);

% Compute Steady State Kalman Gain
L = Sigma_e*C'*inv(C*Sigma_e*C' + R);

% Compute covariance of the residual
Sigma_r = C*Sigma_e*C' + R;

% Feed the initial values for recursion
gamma{1}   = alfa*(eye(n)-L*C)*Sigma_e; 
psi_o{1}   = -L*Bo;
psi_o_0    = -inv(A-L*C*A)*L*Bo;
theta_o{1} = Bo'*inv(Sigma_r)*Bo;
U_o{1}     = -(C*A*L*Bo)'*inv(Sigma_r)*[C*A*psi_o_0 Bo]; 
X_o(:,:,1) = Bo'*inv(Sigma_r)*Bo + (C*A*L*Bo)'*inv(Sigma_r)*C*A*L*Bo;

%%%%%%%% PSI_0 - BIG PROBLEM !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Psi_0 computation
for k = 2:N
    psi_o{k}   = [-((A-L*C*A)^(k-1))*L*Bo  psi_o{k-1}];    
end


% Remember Xc = Uc = Bc = theta_c = psi_c = Vc_{1,2,3} = 0
for k = 2:1:N  
    % Recursive Reachale Set Computation
    X_o(:,:,k)   = (C*A*((A-L*C*A)^(k-1))*L*Bo)'*inv(Sigma_r)*(C*A*((A-L*C*A)^(k-1))*L*Bo);
    U_o{k}       = (C*A*((A-L*C*A)^(k-1))*L*Bo)'*inv(Sigma_r)*[C*A*psi_o{k-1} Bo];
    theta_o{k}   = [X_o(:,:,k-1) U_o{k-1}; U_o{k-1}' theta_o{k-1}];
    V_o1{k-1}    = inv(X_o(:,:,k-1) - U_o{k-1}*inv(theta_o{k-1})*U_o{k-1}');
    V_o2{k-1}    = -inv(X_o(:,:,k-1))*U_o{k-1}*inv(theta_o{k-1});
    V_o3{k-1}    = inv(theta_o{k-1}) - inv(theta_o{k-1})*U_o{k-1}'*V_o1{k-1}*U_o{k-1}*inv(theta_o{k-1}); 
    theta0inv{k} = [V_o1{k-1} V_o2{k-1}; V_o2{k-1}' V_o3{k-1}];
    gamma{k}     = gamma{k-1} + (h-p)*((A-L*C*A)^{k-1}*L*Bo*V_o1{k-1}*Bo'*L'*((A-L*C*A)^{k-1})' + (A-L*C*A)^{k-1}*L*Bo*V_o1{k-1}*U_o{k-1}*inv(theta_o{k-1})*psi_o{k-1}' + ((A-L*C*A)^{k-1}*L*Bo*V_o1{k-1}*U_o{k-1}*inv(theta_o{k-1})*psi_o{k-1}')'  -  psi_o{k-1}*inv(theta_o{k-1})*U_o{k-1}'*V_o1{k-1}*U_o{k-1}*inv(theta_o{k-1})*psi_o{k-1}');

    % Reachable set examination
    mu = max(eig(gamma{k}*P));
    
    % Compute e_a* from optimization problem
    cvx_begin    
        variable e_a(n,1)        
        maximize (e_a'*P*e_a)
        subject to 
        e_a'*inv(gamma{k})*e_a <= 1;    
    cvx_end
    

end


