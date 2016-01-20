% Input:
%     data        a row vector. Mixture of Gaussians.
%     alpha       a scalar. The concentration parameter of the DP.
%     sigma       a scalar. The standard deviation of the Gaussians.
%     actN        a scalar. The maximum number of activated atoms.
%     maxIter     a scalar. The maximum number of gibbs iterations.
% Output:
%     Z           a row vector of the same length with data. Spcify which
%                 cluster the corresponding observation is from.
%     mixing      a row vector of length actN. The mixing measure of the
%                 atoms. 
%     mu          a row vector. The length is length(uniqe(Z)).
function [Z, mixing, mu] = dp_post(data, sigma, alpha, actN, maxIter)
if nargin < 4
    actN = 100;
end
if nargin < 5
    maxIter = 1000;
end
%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------
Z = ones(1, length(data)); % init all in cluster 1
mu = randn(1, actN) * sigma;

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------
for iter = 1:maxIter
    % sample activated centers
    centers = accumarray(Z', data, [], @mean);
    ix = find(centers ~= 0);
    mu(ix) = centers(ix);
    
    % sample mixing measure
    a = histcounts(Z, 1:actN+1);
    b = cumsum(a, 'reverse');
    b = b(2:end);
    b = [b, 0];
    
    a = a + 1;
    b = b + alpha;
    
    V = betarnd(a, b);
    mixing = V;
    V = cumprod(1 - V);
    mixing(2:end) = mixing(2:end) .* V(1:end-1);
    
    % sample Z
    likelihood = normpdf(repmat(data', 1, actN), ...
                         repmat(mu, length(data), 1),...
                         sigma);
    post = likelihood .* repmat(mixing, length(data), 1);
    post = post ./ repmat(sum(post, 2), 1, actN); % normalize
    post = cumsum(post, 2);
    u = rand(length(data), 1);
    for i = 1:length(Z)
        [~, ~, Z(i)] = histcounts(u(i), [0, post(i,:)]);
    end
end

end