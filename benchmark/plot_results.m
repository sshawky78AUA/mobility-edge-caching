function plot_results()
% PLOT_RESULTS  MATLAB version of the Section VI figures (same layout as
% plot_results.py, which rendered the figures in the paper). Run after
% collect_results.m and d2d_predict.m. Writes PDFs to ../figures/.
% Points are means over seeds; error bars are 95% confidence intervals.
S = load('results/summary.mat'); A = S.A; P = A.P;
v = P.speeds; nm = numel(P.models);
col = [138 138 133; 42 120 214; 235 104 52; 27 175 122; 237 161 0; 232 123 164; 98 80 214] / 255;
mk  = {'none','o','s','^','d','v','p'};
ls  = {'--','-','-','-','-','-','-'};
names = P.names; names{7} = 'L\''evy walk';
outdir = fullfile('..', 'paper', 'figures');
if ~isfolder(outdir), outdir = fullfile('..', 'figures'); end
mpc = mean(mean(A.eff(:,:,1)));

% Fig. 2: SBS caching
f = figure('Units', 'inches', 'Position', [1 1 7.16 3.55]);
panels = {A.eff(:,:,2), A.snap(:,:,2), A.eff(:,:,3)};
cis    = {A.eff_ci(:,:,2), A.t95 * A.snap_sd(:,:,2) / sqrt(numel(P.seeds)), A.eff_ci(:,:,3)};
ttl = {'(a) Mobility-aware placement, effective hit ratio', ...
       '(b) Mobility-aware placement, snapshot hit ratio', ...
       '(c) LRU at the serving SBS, effective hit ratio'};
yl = {[0.32 0.54], [0.32 0.54], [0.04 0.26]};
for k = 1:3
  ax = subplot(2, 2, k); hold(ax, 'on');
  for mi = 1:nm
    errorbar(v, panels{k}(mi,:), cis{k}(mi,:), 'LineStyle', ls{mi}, 'Marker', mk{mi}, ...
             'Color', col(mi,:), 'MarkerFaceColor', col(mi,:), 'MarkerSize', 4, 'LineWidth', 1.2);
  end
  if k < 3, yline(mpc, ':', 'Most popular', 'Color', [0.17 0.17 0.16]); end
  ylim(yl{k}); style_axis(ax, v); title(ttl{k}, 'FontWeight', 'normal'); ylabel('Hit ratio');
end
% (d) gain against reach
ax = subplot(2, 2, 4); hold(ax, 'on');
reach = A.stats(:,:,9); gain = 100 * (A.eff(:,:,2) ./ A.eff(:,:,1) - 1);
for mi = 1:nm
  plot(reach(mi,:), gain(mi,:), 'LineStyle', 'none', 'Marker', mk{max(mi,2)}, ...
       'Color', col(mi,:), 'MarkerFaceColor', col(mi,:), 'MarkerSize', 4);
end
fit = setdiff(2:nm, find(strcmp(P.models, 'rwp')));          % moving models except RWP
x = reach(fit,:) - 1; y = gain(fit,:);
kfit = sum(x(:) .* y(:)) / sum(x(:).^2);
xs = linspace(1, 3.4, 50); plot(xs, kfit * (xs - 1), ':', 'Color', [0.17 0.17 0.16]);
text(1.05, 23.5, sprintf('Dotted: %.1f%% per additional SBS', kfit), 'FontSize', 7);
xlim([0.95 3.4]); ylim([-1 27]); grid(ax, 'on'); box(ax, 'off');
xlabel('Reach: SBSs usable within the deadline (mean)'); ylabel('Gain (%)');
title('(d) Gain over most popular against reach', 'FontWeight', 'normal');
lg = legend(names, 'Orientation', 'horizontal', 'Box', 'off', 'Interpreter', 'tex');
lg.Position(1:2) = [0.08 0.0];
exportgraphics(f, fullfile(outdir, 'fig_sbs_hit.pdf'), 'ContentType', 'vector'); close(f);

% Fig. 3: D2D offloading
f = figure('Units', 'inches', 'Position', [1 1 7.16 3.55]);
show = [find(P.Dd_list == 2), find(P.Dd_list == 10), find(P.Dd_list == 30)];
yl = {[0.14 0.50], [0.11 0.27], [0.07 0.19]};
for k = 1:3
  ax = subplot(2, 2, k); hold(ax, 'on');
  for mi = 1:nm
    errorbar(v, A.d2d(mi,:,show(k)), A.d2d_ci(mi,:,show(k)), 'LineStyle', ls{mi}, 'Marker', mk{mi}, ...
             'Color', col(mi,:), 'MarkerFaceColor', col(mi,:), 'MarkerSize', 4, 'LineWidth', 1.2);
  end
  ylim(yl{k}); style_axis(ax, v); ylabel('D2D offloading ratio');
  title(sprintf('(%s) Transfer time D_d = %d s', char('a'+k-1), P.Dd_list(show(k))), 'FontWeight', 'normal');
end
% (d) every configuration against g(Dd), with the prediction of Eq. (5)
ax = subplot(2, 2, 4); hold(ax, 'on');
T = d2d_predict(); R = T.rows;                    % model, speed, seed, Dd, g, pred, sim
[keys, ~, id] = unique(R(:, [1 2 4]), 'rows');
gm = accumarray(id, R(:,5), [], @mean); sm = accumarray(id, R(:,7), [], @mean);
for mi = 1:nm
  sel = keys(:,1) == mi;
  plot(gm(sel), sm(sel), 'LineStyle', 'none', 'Marker', mk{max(mi,2)}, 'Color', col(mi,:), ...
       'MarkerFaceColor', col(mi,:), 'MarkerSize', 3.5);
end
p = zipf_pmf(P.N, P.gamma); p = p(:);
gs = logspace(-3.3, -0.7, 80); q = inclusion_prob_local(p, P.c);
plot(gs, arrayfun(@(g) sum(p .* (q + (1 - q) .* (1 - (1 - g*q).^(P.U - 1)))), gs), 'k-', 'LineWidth', 1);
set(ax, 'XScale', 'log'); xlim([5e-4 0.2]); ylim([0.07 0.52]); grid(ax, 'on'); box(ax, 'off');
xlabel('Usable-contact probability per pair, g(D_d)'); ylabel('D2D offloading ratio');
title('(d) All runs against g(D_d); solid line: prediction', 'FontWeight', 'normal');
lg = legend(names, 'Orientation', 'horizontal', 'Box', 'off');
lg.Position(1:2) = [0.08 0.0];
exportgraphics(f, fullfile(outdir, 'fig_d2d.pdf'), 'ContentType', 'vector'); close(f);
end

function q = inclusion_prob_local(p, c)
set_seed(1); reps = 40000; cnt = zeros(size(p));
for r = 1:reps
  w = p;
  for k = 1:c
    f = find(cumsum(w)/sum(w) >= rand, 1, 'first'); cnt(f) = cnt(f) + 1; w(f) = 0;
  end
end
q = cnt / reps;
end

function style_axis(ax, v)
set(ax, 'XScale', 'log', 'XTick', v, 'XTickLabel', arrayfun(@(x) sprintf('%g', x), v, 'UniformOutput', false), ...
    'XMinorTick', 'off', 'FontSize', 8);
xlim(ax, [0.4 25]); grid(ax, 'on'); box(ax, 'off');
xlabel(ax, 'Travel speed v (m/s)');
end
