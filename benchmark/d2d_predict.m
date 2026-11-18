function T = d2d_predict()
% D2D_PREDICT  Contact-based prediction of the D2D offloading ratio (Sec. VI-D).
% For every run of the main sweep and every transfer time Dd, computes the
% usable-contact probability g(Dd) of Eq. (4) from the run's contact-length
% histogram and the predicted offloading ratio of Eq. (5), and compares it
% with the simulated ratio. Returns a table-like struct and writes
% results/d2d_prediction_matlab.csv.
P = bench_config();
p = zipf_pmf(P.N, P.gamma); p = p(:);
q = inclusion_prob(p, P.c, 40000);
rows = [];
for mi = 1:numel(P.models)
  for v = P.speeds
    for s = P.seeds
      r = load(sprintf('results/main_%s_v%g_s%d.mat', P.models{mi}, v, s)); r = r.res;
      h = r.stats.contact_hist(:); l = (1:numel(h)).';
      for k = 1:numel(P.Dd_list)
        Dd = P.Dd_list(k);
        if strcmp(P.models{mi}, 'static')
          g = r.stats.contact_frac;            % pairs in range stay in range
        else
          w = (l >= Dd) .* ((l - Dd + 1) + (P.Td - Dd));
          g = sum(h .* w) / r.stats.n_pairs / P.Tsim;
        end
        pred = sum(p .* (q + (1 - q) .* (1 - (1 - g*q).^(P.U - 1))));
        rows = [rows; mi, v, s, Dd, g, pred, r.d2d(k)]; %#ok<AGROW>
      end
    end
  end
end
T.rows = rows;
err = rows(:,6) - rows(:,7);
fprintf('runs x Dd: %d, mean |error| %.4f, max |error| %.4f\n', size(rows,1), mean(abs(err)), max(abs(err)));
fid = fopen('results/d2d_prediction_matlab.csv', 'w');
fprintf(fid, 'model,speed,seed,Dd,g,pred,sim\n');
for k = 1:size(rows, 1)
  fprintf(fid, '%s,%g,%d,%d,%.6g,%.6g,%.6g\n', P.models{rows(k,1)}, rows(k,2:end));
end
fclose(fid);
end

function q = inclusion_prob(p, c, reps)
% Probability that a user holds each file under sequential sampling of c
% distinct files with probability proportional to popularity (Monte Carlo).
set_seed(1);
cnt = zeros(size(p));
for r = 1:reps
  w = p;
  for k = 1:c
    f = find(cumsum(w)/sum(w) >= rand, 1, 'first'); cnt(f) = cnt(f) + 1; w(f) = 0;
  end
end
q = cnt / reps;
end
