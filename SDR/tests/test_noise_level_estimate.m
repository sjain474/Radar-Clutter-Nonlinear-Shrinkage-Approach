function test_noise_level_estimate

sigma = 7.2
n = 200;
p = 100;
Y = sigma * randn(n,p) / sqrt(n);
S = Y'*Y;
[V D] = eig(S);
d = diag(D);
[d sigmahat] = optimal_shrinkage(d , p/n ,'F_1');  
Shat = V * diag(d) * V';

rel_err = abs(sigma-sigmahat) / sigma
