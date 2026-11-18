function p = zipf_pmf(N, gamma)
% ZIPF_PMF  Zipf popularity over N files with exponent gamma (file 1 most popular).
p = (1:N).^(-gamma);
p = p / sum(p);
end
