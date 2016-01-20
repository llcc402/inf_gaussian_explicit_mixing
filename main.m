clear
clc

K = 5;
N = 500;
sigma = 1;
alpha = .1;

[data, centers, prob] = data_generate(K, N, sigma);
tic;
[Z, mixing, mu] = dp_post(data, sigma, alpha);
toc

figure(2)
[~, ix] = sort(mu);
plot(centers, prob, 'o', mu(ix), mixing(ix), '*')
xlim([-10, 10])
hold on
line([centers; centers], [zeros(1, length(centers)); prob], 'color', 'blue')
line([mu; mu], [zeros(1, length(mu)); mixing], 'color', 'red')
legend('theretical', 'sampled')
title('The theoretical and the sampled mixing measures')
hold off