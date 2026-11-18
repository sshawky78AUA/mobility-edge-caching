function Xc = place_mobility_aware(W, piB, p, P)
% PLACE_MOBILITY_AWARE  Greedy SBS placement that uses the reach matrix W.
% Objective: sum_b piB(b) sum_f p(f) [1 - prod_b2 (1 - x(b2,f) W(b,b2))],
% the probability that a request finds its file at some SBS it stays with for
% P.Ds seconds before the deadline, with visits treated as independent. The
% greedy step adds the (SBS, file) pair with the largest marginal gain until
% every cache holds P.C files (cf. the helper placement of Poularakis and
% Tassiulas, ISIT 2013, here without coding).
B = P.nSBS; N = P.N;
p  = p(:).';
M  = ones(B, N);                 % miss probability per origin SBS and file
Xc = false(B, N);
load = zeros(B, 1);
for k = 1:B*P.C
  G = (W.' * (repmat(piB, 1, N) .* M)) .* repmat(p, B, 1);  % gain of (b2, f)
  G(Xc) = -Inf;
  G(load >= P.C, :) = -Inf;
  [~, idx] = max(G(:));
  [b2, f] = ind2sub([B, N], idx);
  Xc(b2, f) = true;  load(b2) = load(b2) + 1;
  M(:, f) = M(:, f) .* (1 - W(:, b2));
end
end
