function plot_results()
% PLOT_RESULTS  MATLAB version of the Section VI figures (same panels as
% plot_results.py, which rendered the figures in the paper). Run after
% collect_results.m and d2d_predict.m. Writes fig_sbs_hit.pdf, fig_d2d.pdf and
% fig_predict.pdf to ../paper/figures/.
% Curves are means over seeds; shaded bands are 95% confidence intervals.
P = bench_config();
[R, hdr] = read_raw('results/raw.csv');
H = read_hist('results/contact_hist.csv');
c = @(name) R.(name);
models = P.models; nm = numel(models); v = P.speeds; ns = numel(P.seeds);
t95 = tinv_local(ns - 1);
col = [138 138 133; 42 120 214; 235 104 52; 27 175 122; 237 161 0; 232 123 164; 98 80 214] / 255;
mk  = {'o','o','s','^','d','v','p'};
ddc = [158 197 240; 94 159 230; 42 120 214; 26 85 163; 13 50 102] / 255;
names = {'Static','Random walk','RWP','Gauss-Markov','Markov cell','Manhattan','L\''evy walk'};
outdir = fullfile('..', 'paper', 'figures');
if ~isfolder(outdir), outdir = fullfile('..', 'figures'); end

% per (model, speed) means and CI half-widths
mean_of = @(x) cellfun(@mean, x); ci_of = @(x) t95 * cellfun(@std, x) / sqrt(ns);
cells = @(x) by_group(x, R, nm, v);
G  = cells(100 * (c('eff_ma') ./ c('eff_mpc') - 1));
Gs = cells(100 * (c('snap_ma') ./ c('snap_mpc') - 1));
reach = mean_of(cells(c('mean_reach'))); pstay = mean_of(cells(c('p_stay_first')));
lru = mean_of(cells(c('eff_lru'))); snap_lru = mean(c('snap_lru'));

% ---------------- Fig. 2: SBS caching
f = figure('Units', 'inches', 'Position', [1 1 7.16 2.55]);
ax = subplot(1, 3, 1); hold(ax, 'on');
for mi = 1:nm
  band(ax, v, mean_of(G(mi,:)), ci_of(G(mi,:)), col(mi,:), mi == 1);
end
hline(ax, 0, [0.4 25]); ylim([-2 28]); style_speed(ax);
ylabel('Gain over most popular (%)');
title('(a) Gain over MPC against speed', 'FontWeight', 'normal');

ax = subplot(1, 3, 2); hold(ax, 'on');
patch([-26 0 0 -26], [-2 -2 28 28], [0.98 0.93 0.91], 'EdgeColor', 'none');
for mi = 2:nm
  plot(mean_of(Gs(mi,:)), mean_of(G(mi,:)), mk{mi}, 'Color', col(mi,:), ...
       'MarkerFaceColor', col(mi,:), 'MarkerSize', 4);
end
plot([-30 30], [-30 30], ':', 'Color', [0.42 0.42 0.4]); plot([0 0], [-2 28], 'k-'); hline(ax, 0, [-26 3]);
text(-25.5, 4, {'snapshot predicts a loss,', 'delivery shows a gain'}, 'FontSize', 7, ...
     'Color', [0.69 0.29 0.12]);
xlim([-26 3]); ylim([-2 28]); grid(ax, 'on'); box(ax, 'off');
xlabel('Gain under the snapshot hit ratio (%)'); ylabel('Gain under the effective hit ratio (%)');
title('(b) Snapshot vs effective gain', 'FontWeight', 'normal');

ax = subplot(1, 3, 3); hold(ax, 'on');
for mi = 1:nm
  plot(pstay(mi,:), lru(mi,:), mk{mi}, 'Color', col(mi,:), 'MarkerFaceColor', col(mi,:), 'MarkerSize', 4);
end
plot([0 1], snap_lru * [0 1], 'k-');
text(0.55, 0.075, 'line: LRU snapshot hit ratio \times P(stay \geq D_s)', 'FontSize', 7);
xlim([0.2 1.02]); ylim([0.04 0.26]); grid(ax, 'on'); box(ax, 'off');
xlabel('P(user stays \geq D_s in the first cell)'); ylabel('LRU effective hit ratio');
title('(c) LRU against handover risk', 'FontWeight', 'normal');

legend_models(f, col, mk, names, 1:nm);
save_pdf(f, fullfile(outdir, 'fig_sbs_hit.pdf'));

% ---------------- Fig. 3: D2D offloading
f = figure('Units', 'inches', 'Position', [1 1 7.16 2.55]);
ax = subplot(1, 3, 1); hold(ax, 'on');
rw = find(strcmp(models, 'randwalk'));
for k = 1:numel(P.Dd_list)
  key = sprintf('d2d_Dd%d', P.Dd_list(k)); s0 = mean(R.(key)(R.mi == 1));
  Q = cells(100 * (c(key) / s0 - 1));
  band(ax, v, mean_of(Q(rw,:)), ci_of(Q(rw,:)), ddc(k,:), false);
end
hline(ax, 0, [0.4 25]); ylim([-60 220]); style_speed(ax); ylabel('Offloading vs static (%)');
legend(arrayfun(@(d) sprintf('D_d = %d s', d), P.Dd_list, 'UniformOutput', false), ...
       'Location', 'northwest', 'Box', 'off');
title('(a) Offloading vs static, random walk', 'FontWeight', 'normal');

ax = subplot(1, 3, 2); hold(ax, 'on');
L = 1:size(H.h, 2); vs = [1 5 20];
for k = 1:3
  cnt = sum(H.h(H.mi == rw & H.speed == vs(k), :), 1);
  plot(L, fliplr(cumsum(fliplr(cnt))) / sum(cnt), 'Color', ddc(2*k - 1,:), 'LineWidth', 1.4);
end
for d = P.Dd_list
  plot([d d], [3e-3 1], ':', 'Color', [0.42 0.42 0.4], 'HandleVisibility', 'off');
  text(d, 1.25, sprintf('%d', d), 'FontSize', 7, 'HorizontalAlignment', 'center');
end
set(ax, 'XScale', 'log', 'YScale', 'log'); xlim([1 600]); ylim([3e-3 1]); grid(ax, 'on'); box(ax, 'off');
legend(arrayfun(@(s) sprintf('v = %g m/s', s), vs, 'UniformOutput', false), 'Location', 'southwest', 'Box', 'off');
xlabel('Contact duration l (s)'); ylabel('Share of contacts lasting \geq l');
title('(b) Contact durations, random walk', 'FontWeight', 'normal');

ax = subplot(1, 3, 3); hold(ax, 'on');
B = best_speeds(R, P);                 % [mi, Dd, vstar, interior]
for mi = 2:6
  b = B(B(:,1) == mi, :);
  plot(b(:,2), b(:,3), '-', 'Marker', mk{mi}, 'Color', col(mi,:), 'MarkerFaceColor', col(mi,:), ...
       'MarkerSize', 4, 'LineWidth', 1.2);
end
use = B(:,1) >= 2 & B(:,1) <= 6 & B(:,4) == 1;
kap = exp(mean(log(B(use,3) .* B(use,2))));
xs = linspace(1.8, 34, 50); plot(xs, kap ./ xs, 'k:');
text(2.1, 0.62, sprintf('dotted: v^\\star D_d = %.0f m', kap), 'FontSize', 7);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', P.Dd_list, 'YTick', [0.5 1 2 5 10 20]);
xlim([1.7 36]); ylim([0.4 26]); grid(ax, 'on'); box(ax, 'off');
xlabel('Transfer time D_d (s)'); ylabel('Best speed v^\star (m/s)');
title('(c) Best speed per transfer time', 'FontWeight', 'normal');

legend_models(f, col, mk, names, 2:6);
save_pdf(f, fullfile(outdir, 'fig_d2d.pdf'));

% ---------------- Fig. 4: what predicts the outcome
f = figure('Units', 'inches', 'Position', [1 1 7.16 2.45]);
ax = subplot(1, 2, 1); hold(ax, 'on');
gm = mean_of(G);
for mi = 1:nm
  plot(reach(mi,:), gm(mi,:), mk{mi}, 'Color', col(mi,:), 'MarkerFaceColor', col(mi,:), 'MarkerSize', 4);
end
fit = setdiff(2:nm, find(strcmp(models, 'rwp')));          % moving models except RWP
x = reach(fit,:) - 1; y = gm(fit,:);
kfit = sum(x(:) .* y(:)) / sum(x(:).^2); rc = corrcoef(x(:), y(:));
xs = linspace(1, 3.4, 50); plot(xs, kfit * (xs - 1), 'k:');
text(1.03, 25, sprintf('%.1f%% per additional SBS (r = %.2f)', kfit, rc(1,2)), 'FontSize', 7);
xlim([0.95 3.4]); ylim([-2 28]); grid(ax, 'on'); box(ax, 'off');
xlabel('Reach: SBSs usable within the deadline (mean)'); ylabel('Gain over most popular (%)');
title('(a) SBS caching: reach predicts the gain', 'FontWeight', 'normal');
legend(names, 'Location', 'northwest', 'Box', 'off');
ax = subplot(1, 2, 2); hold(ax, 'on');
T = d2d_predict(); X = T.rows;                     % model, speed, seed, Dd, g, pred, sim
[keys, ~, id] = unique(X(:, [1 2 4]), 'rows');
pm = accumarray(id, X(:,6), [], @mean); sm = accumarray(id, X(:,7), [], @mean);
for k = 1:numel(P.Dd_list)
  sel = keys(:,3) == P.Dd_list(k);
  plot(pm(sel), sm(sel), 'o', 'Color', ddc(k,:), 'MarkerFaceColor', ddc(k,:), 'MarkerSize', 3);
end
plot([0 0.6], [0 0.6], 'k-');
e = abs(pm - sm);
text(0.075, 0.46, {sprintf('%d configurations', numel(e)), ...
     sprintf('mean |error| %.3f, max %.3f', mean(e), max(e))}, 'FontSize', 7);
xlim([0.07 0.52]); ylim([0.07 0.52]); grid(ax, 'on'); box(ax, 'off');
xlabel('Predicted offloading ratio, Eq. (5)'); ylabel('Simulated offloading ratio');
title('(b) D2D caching: contact statistics predict offloading', 'FontWeight', 'normal');
save_pdf(f, fullfile(outdir, 'fig_predict.pdf'));
end

% ------------------------------------------------------------------ helpers
function G = by_group(x, R, nm, v)
% Cell array (model x speed) of the per-seed values of x.
G = cell(nm, numel(v));
for mi = 1:nm
  for k = 1:numel(v), G{mi,k} = x(R.mi == mi & R.speed == v(k)); end
end
end

function B = best_speeds(R, P)
% Speed that maximizes offloading, refined by a parabola in log v through the
% best grid point and its neighbours. Columns: model, Dd, v*, interior flag.
B = [];
x = log(P.speeds(:));
for mi = 2:numel(P.models)
  for d = P.Dd_list
    y = arrayfun(@(s) mean(R.(sprintf('d2d_Dd%d', d))(R.mi == mi & R.speed == s)), P.speeds(:));
    [~, i] = max(y); xb = x(i); inner = i > 1 && i < numel(y);
    if inner
      pp = polyfit(x(i-1:i+1), y(i-1:i+1), 2);
      if pp(1) < 0, xb = min(max(-pp(2) / (2 * pp(1)), x(i-1)), x(i+1)); end
    end
    B(end+1, :) = [mi, d, exp(xb), inner]; %#ok<AGROW>
  end
end
end

function hline(ax, y, xr)
plot(ax, xr, [y y], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

function band(ax, x, m, e, colr, dashed)
fill(ax, [x fliplr(x)], [m - e, fliplr(m + e)], colr, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
     'HandleVisibility', 'off');
if dashed, st = '--'; else, st = '-'; end
plot(ax, x, m, st, 'Color', colr, 'LineWidth', 1.4);
end

function style_speed(ax)
t = [0.5 1 2 5 10 20];
set(ax, 'XScale', 'log', 'XTick', t, 'XTickLabel', arrayfun(@(x) sprintf('%g', x), t, ...
    'UniformOutput', false), 'XMinorTick', 'off', 'FontSize', 8);
xlim(ax, [0.45 22]); grid(ax, 'on'); box(ax, 'off'); xlabel(ax, 'Travel speed v (m/s)');
end

function legend_models(f, col, mk, names, idx)
ax = axes(f, 'Position', [0 0 1 0.04], 'Visible', 'off'); hold(ax, 'on');
h = zeros(numel(idx), 1);
for k = 1:numel(idx)
  h(k) = plot(ax, NaN, NaN, '-', 'Marker', mk{idx(k)}, 'Color', col(idx(k),:), ...
              'MarkerFaceColor', col(idx(k),:), 'MarkerSize', 4);
end
legend(h, names(idx), 'Orientation', 'horizontal', 'Box', 'off', 'Location', 'south');
end

function save_pdf(f, fn)
if exist('exportgraphics', 'file'), exportgraphics(f, fn, 'ContentType', 'vector');
else, set(f, 'PaperPositionMode', 'auto'); print(f, fn, '-dpdf'); end
close(f);
end

function t = tinv_local(df)
tab = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262];
if df <= numel(tab), t = tab(df); else, t = 1.96; end
end

function [R, hdr] = read_raw(fn)
fid = fopen(fn); hdr = strsplit(fgetl(fid), ','); fclose(fid);
fmt = ['%s' repmat('%f', 1, numel(hdr) - 1)];
fid = fopen(fn); C = textscan(fid, fmt, 'Delimiter', ',', 'HeaderLines', 1, 'TreatAsEmpty', 'NaN'); fclose(fid);
P = bench_config();
R = struct('mi', cellfun(@(m) find(strcmp(P.models, m)), C{1}));
for k = 2:numel(hdr), R.(hdr{k}) = C{k}; end
end

function H = read_hist(fn)
P = bench_config();
fid = fopen(fn); first = fgetl(fid); fclose(fid);
nc = numel(strsplit(first, ','));
fid = fopen(fn); C = textscan(fid, ['%s' repmat('%f', 1, nc - 1)], 'Delimiter', ','); fclose(fid);
H.mi = cellfun(@(m) find(strcmp(P.models, m)), C{1});
H.speed = C{2}; H.h = cell2mat(C(5:end));
end
