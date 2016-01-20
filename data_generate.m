% Input:
%     K        a scalar. The number of Gaussians.
%     N        a scalar. The number of observations in each cluster.
%     sigma    a scalar. The standard deviation of each Gaussian.
% Output:
%     data     a row vector of length K * N.
%     z        a row vector of the same length with data, the indicators.
function [data, mu, prob] = data_generate(K, N, sigma)
if nargin < 1
    K = 5;
end
if nargin < 2
    N = 500;
end
if nargin < 3
    sigma = 1;
end

%--------------------------------------------------------------------------
% STEP 1: Generate the mixing measure
%--------------------------------------------------------------------------

% positions
mu = rand(1, K) + 0 : 4*sigma : 4*sigma*K;
mu = mu - mean(mu);

% weights
prob = rand(1,K);
prob = prob / sum(prob);

%--------------------------------------------------------------------------
% STEP 2: Generate latent variable z
%--------------------------------------------------------------------------
u = rand(1, K * N);
[~, ~, z] = histcounts(u, [0, cumsum(prob)]);

%--------------------------------------------------------------------------
% STEP 3: Generate data
%--------------------------------------------------------------------------
data = randn(1, K * N) * sigma + mu(z);
hist(data, K*20)
end