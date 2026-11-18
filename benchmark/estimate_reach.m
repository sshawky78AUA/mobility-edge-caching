function [W, piB] = estimate_reach(S, P)
% ESTIMATE_REACH  Empirical reach matrix of a trajectory set.
%   W(b,b2) estimates the probability that a user associated with SBS b at a
%   request instant stays with SBS b2 for at least P.Ds consecutive seconds
%   within the next P.Td seconds. piB(b) is the share of time spent at b.
%   This is the only mobility information the mobility-aware policy uses.
B = P.nSBS;
W = zeros(B, B); cnt = zeros(B, 1);
for u = 1:size(S,1)
  for t = 1:P.ma_step:P.Tsim
    win = S(u, t:t+P.Td);
    [vals, lens] = run_lengths(win);
    b0 = win(1);
    ok = unique(vals(lens >= P.Ds));
    W(b0, ok) = W(b0, ok) + 1;
    cnt(b0) = cnt(b0) + 1;
  end
end
cnt(cnt == 0) = 1;
W   = W ./ repmat(cnt, 1, B);
piB = histc(reshape(S(:, 1:P.Tsim), [], 1), 1:B);
piB = piB(:) / sum(piB);
end
