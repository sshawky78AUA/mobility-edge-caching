function R = gen_requests(P, seed)
% GEN_REQUESTS  Poisson requests per user with Zipf file choice.
%   R is a K-by-3 matrix [user, time, file], sorted by time. Time is the
%   second at which the request arrives, in 0..Tsim-1.
set_seed(10000 + seed);
p   = zipf_pmf(P.N, P.gamma);
cdf = cumsum(p);
R = zeros(0, 3);
for u = 1:P.U
  n = poissrnd_simple(P.lambda * P.Tsim);
  t = floor(rand(n, 1) * P.Tsim);
  f = arrayfun(@(r) find(cdf >= r, 1, 'first'), rand(n, 1));
  R = [R; [u*ones(n,1), t, f(:)]]; %#ok<AGROW>
end
R = sortrows(R, 2);
end

function n = poissrnd_simple(mu)
% Poisson sample by inversion (avoids the Statistics toolbox).
L = exp(-mu); k = 0; p = 1;
if mu > 50
  n = max(0, round(mu + sqrt(mu)*randn)); return;
end
while p > L, k = k + 1; p = p * rand; end
n = k - 1;
end
