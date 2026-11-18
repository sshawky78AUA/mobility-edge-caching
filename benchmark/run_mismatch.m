function run_mismatch(v, seeds)
% RUN_MISMATCH  Train the mobility-aware placement on model i, test it on model j.
% Uses the reach matrices W stored by run_sweep (trained on independent runs)
% and reports the effective SBS hit ratio for every (train, test) pair.
% run_mismatch(v) uses every seed; run_mismatch(v, seeds) a subset, saved to
% a part file that merge_mismatch combines.
if nargin < 1, v = 5; end
P = bench_config();
part = nargin >= 2;
if part, P.seeds = seeds; end
nm = numel(P.models);
p = zipf_pmf(P.N, P.gamma);
H = zeros(nm, nm, numel(P.seeds));
Hmpc = zeros(nm, numel(P.seeds));
noD2D = false(P.U, P.N);
Xmpc = repmat((1:P.N) <= P.C, P.nSBS, 1);
for si = 1:numel(P.seeds)
  s = P.seeds(si);
  % placements from every training model
  Xma = cell(nm, 1);
  for i = 1:nm
    r = load(sprintf('results/main_%s_v%g_s%d.mat', P.models{i}, v, s));
    tr = r.res;
    [Xt, Yt] = gen_mobility(P.models{i}, v, P, 1000 + s); %#ok<ASGLU>
    piB = histc(reshape(sbs_assoc(Xt, Yt, P), [], 1), 1:P.nSBS); piB = piB(:)/sum(piB);
    Xma{i} = place_mobility_aware(tr.W, piB, p, P);
  end
  R = gen_requests(P, s);
  for j = 1:nm
    [X, Y] = gen_mobility(P.models{j}, v, P, s);
    S = sbs_assoc(X, Y, P);
    for i = 1:nm
      out = eval_requests(X, Y, S, R, Xmpc, Xma{i}, noD2D, P);
      H(i, j, si) = out.eff(2);
      Hmpc(j, si) = out.eff(1);
    end
    fprintf('mismatch seed %d test %s done\n', s, P.models{j});
  end
end
if part
  save('-v7', sprintf('results/mismatch_v%g_part%d.mat', v, seeds(1)), 'H', 'Hmpc', 'P');
else
  save('-v7', sprintf('results/mismatch_v%g.mat', v), 'H', 'Hmpc', 'P');
end
end
