function res = run_one(mi, v, seed, P)
% RUN_ONE  One (model, speed, seed) point of the main sweep.
model = P.models{mi};
p = zipf_pmf(P.N, P.gamma);
% Mobility-aware placement is trained on an independent run of the same model
[Xt, Yt] = gen_mobility(model, v, P, 1000 + seed);
[W, piB] = estimate_reach(sbs_assoc(Xt, Yt, P), P);
Xma  = place_mobility_aware(W, piB, p, P);
Xmpc = repmat((1:P.N) <= P.C, P.nSBS, 1);
% User caches: each user stores c distinct files drawn in proportion to popularity
set_seed(20000 + seed);
Yd2d = false(P.U, P.N);
for u = 1:P.U
  q = p; for k = 1:P.c
    f = find(cumsum(q)/sum(q) >= rand, 1, 'first'); Yd2d(u, f) = true; q(f) = 0;
  end
end
[X, Y] = gen_mobility(model, v, P, seed);
S = sbs_assoc(X, Y, P);
R = gen_requests(P, seed);
out = eval_requests(X, Y, S, R, Xmpc, Xma, Yd2d, P);
st  = mobility_stats(X, Y, S, P, seed);
res = struct('model', model, 'v', v, 'seed', seed, 'snap', out.snap, ...
             'eff', out.eff, 'd2d', out.d2d, 'K', out.K, 'stats', st, 'W', W);
end
