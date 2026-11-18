function run_mixed(speeds, seeds)
% RUN_MIXED  Mixed population: P.mix_static of the users are static and the
% rest follow a random walk at speed v. Compares, on the same test runs,
%   MPC, the mobility-aware placement trained on the mixed population, and
%   the mobility-aware placement trained on an all-moving random walk
% (a mismatch that ignores the static users), plus LRU and D2D offloading.
% Results: results/mixed.mat and results/mixed.csv (one row per speed and seed).
P = bench_config();
if nargin < 1, speeds = [1 2 5 10 20]; end
if nargin < 2, seeds = 1:5; end
p = zipf_pmf(P.N, P.gamma);
Xmpc = repmat((1:P.N) <= P.C, P.nSBS, 1);
rows = [];
for v = speeds
  for s = seeds
    [Xt, Yt] = gen_mobility('mix_randwalk', v, P, 1000 + s);
    [Wm, pm] = estimate_reach(sbs_assoc(Xt, Yt, P), P);
    Xma_mix = place_mobility_aware(Wm, pm, p, P);
    [Xt, Yt] = gen_mobility('randwalk', v, P, 1000 + s);
    [Wr, pr] = estimate_reach(sbs_assoc(Xt, Yt, P), P);
    Xma_mov = place_mobility_aware(Wr, pr, p, P);
    set_seed(20000 + s);
    Yd2d = false(P.U, P.N);
    for u = 1:P.U
      q = p; for k = 1:P.c
        f = find(cumsum(q)/sum(q) >= rand, 1, 'first'); Yd2d(u, f) = true; q(f) = 0;
      end
    end
    [X, Y] = gen_mobility('mix_randwalk', v, P, s);
    S = sbs_assoc(X, Y, P);
    R = gen_requests(P, s);
    a = eval_requests(X, Y, S, R, Xmpc, Xma_mix, Yd2d, P);
    b = eval_requests(X, Y, S, R, Xmpc, Xma_mov, Yd2d, P);
    rows = [rows; v, s, a.eff(1), a.eff(2), b.eff(2), a.eff(3), a.d2d(:).']; %#ok<AGROW>
    fprintf('mixed v=%g s=%d: MPC %.3f MA(mixed) %.3f MA(moving) %.3f\n', v, s, a.eff(1), a.eff(2), b.eff(2));
  end
end
save('-v7', 'results/mixed.mat', 'rows', 'P');
fid = fopen('results/mixed.csv', 'w');
fprintf(fid, 'speed,seed,eff_mpc,eff_ma_mixed,eff_ma_moving,eff_lru');
fprintf(fid, ',d2d_Dd%d', P.Dd_list); fprintf(fid, '\n');
for k = 1:size(rows, 1), fprintf(fid, '%g,%d', rows(k,1:2)); fprintf(fid, ',%.6g', rows(k,3:end)); fprintf(fid, '\n'); end
fclose(fid);
end
