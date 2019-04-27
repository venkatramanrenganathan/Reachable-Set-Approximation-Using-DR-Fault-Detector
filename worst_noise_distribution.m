function sample = worst_noise_distribution(input_param)
n     = size(input_param.Sigma_w,2); % dimension
x_bar = input_param.mu_noise;
S     = input_param.Sigma_w;
A     = inv(S);
output_param = compute_noise_threshold(input_param);
z      = output_param.z;
Z      = output_param.Z;
lambda = output_param.lambda;
c      = -output_param.noise_threshold; 
x0     = 1/(1-lambda)*(x_bar - lambda*z/lambda);
S0     = 1/(1-lambda)*(S - lambda*Z/lambda);
 
%% Construction of Random variable
Z = Z/lambda;
z = z/lambda;
 
[V,E] = eig(Z - z*z');
lam   = z'*z + c;
W     = V*sqrt(E);
w1    = W(:,1);
w2    = W(:,2);
 
mu1 = w1'*w1;
mu2 = w2'*w2;
 
beta1 = roots([mu1 2*w1'*z lam]);
beta2 = roots([mu2 2*w2'*z lam]);
 
v1  = z + beta1(1)*w1;
v1r = z + beta1(2)*w1;
v2  = z + beta2(1)*w2;
v2r = z + beta2(2)*w2;
 
alpha1  = mu1/(1-beta1(1)/beta1(2))/(mu1+mu2);
alpha2  = mu2/(1-beta2(1)/beta2(2))/(mu1+mu2);
alpha1r = -alpha1*beta1(1)/beta1(2);
alpha2r = -alpha2*beta2(1)/beta2(2);


%% discrete distribution: x0 with probability 1-lambda, v2 with probability alpha2*lambda, v2r with probability alpha2r*lambda
sample.probabilities = [1-lambda alpha1*lambda alpha1r*lambda alpha2*lambda alpha2r*lambda];
sample.support       = [x0 v1 v1r v2 v2r];
 
end