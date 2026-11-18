function run_sweep(part, nparts)
% RUN_SWEEP  Main sweep: every model x speed x seed. Results go to results/.
% Run with run_sweep(1,1), or split across processes with run_sweep(k, n).
P = bench_config();
jobs = [];
for mi = 1:numel(P.models)
  for v = P.speeds
    for s = P.seeds
      jobs = [jobs; mi, v, s]; %#ok<AGROW>
    end
  end
end
for j = part:nparts:size(jobs,1)
  mi = jobs(j,1); v = jobs(j,2); s = jobs(j,3);
  fn = sprintf('results/main_%s_v%g_s%d.mat', P.models{mi}, v, s);
  if exist(fn, 'file'), continue; end
  res = run_one(mi, v, s, P); %#ok<NASGU>
  save('-v7', fn, 'res');
  fprintf('%s done\n', fn);
end
end
