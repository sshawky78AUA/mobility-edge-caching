function merge_mismatch(v)
% MERGE_MISMATCH  Combine the part files written by run_mismatch(v, seeds).
if nargin < 1, v = 5; end
f = dir(sprintf('results/mismatch_v%g_part*.mat', v));
H = []; Hmpc = []; seeds = [];
for k = 1:numel(f)
  r = load(fullfile('results', f(k).name));
  H = cat(3, H, r.H); Hmpc = [Hmpc, r.Hmpc]; seeds = [seeds, r.P.seeds]; %#ok<AGROW>
end
[~, o] = sort(seeds); H = H(:,:,o); Hmpc = Hmpc(:,o);
P = bench_config();
save('-v7', sprintf('results/mismatch_v%g.mat', v), 'H', 'Hmpc', 'P');
end
