function A = collect_results()
% COLLECT_RESULTS  Average the main sweep over seeds and write CSV summaries.
% Output arrays are models x speeds; *_sd hold the standard deviation over seeds.
P = bench_config();
nm = numel(P.models); nv = numel(P.speeds); ns = numel(P.seeds); nD = numel(P.Dd_list);
snap = zeros(nm, nv, 3, ns); eff = snap; d2d = zeros(nm, nv, nD, ns);
st = zeros(nm, nv, 9, ns); lr = zeros(nm, nv, nD, ns); raw = []; hist = [];
for mi = 1:nm
  for vi = 1:nv
    for si = 1:ns
      r = load(sprintf('results/main_%s_v%g_s%d.mat', P.models{mi}, P.speeds(vi), P.seeds(si)));
      r = r.res;
      snap(mi,vi,:,si) = r.snap; eff(mi,vi,:,si) = r.eff; d2d(mi,vi,:,si) = r.d2d;
      s = r.stats;
      st(mi,vi,:,si) = [s.mean_speed, s.mean_sojourn, s.handovers_per_hour, ...
                        s.mean_contact, s.mean_intercontact, s.contact_frac, ...
                        s.n_contacts_per_pair_hour, s.p_stay_first, s.mean_reach];
      % rate of complete contacts lasting at least Dd, per pair-hour
      for i = 1:nD
        lr(mi,vi,i,si) = sum(s.contact_hist(P.Dd_list(i):end)) / s.n_pairs * 3600 / P.Tsim;
      end
      hist(end+1,:) = [mi, P.speeds(vi), P.seeds(si), s.n_pairs, s.contact_hist]; %#ok<AGROW>
      raw(end+1,:) = [mi, P.speeds(vi), P.seeds(si), r.snap, r.eff, r.d2d, ...
                      squeeze(st(mi,vi,:,si)).', squeeze(lr(mi,vi,:,si)).']; %#ok<AGROW>
    end
  end
end
A.P = P;
A.snap = mean(snap, 4); A.snap_sd = std(snap, 0, 4);
A.eff  = mean(eff, 4);  A.eff_sd  = std(eff, 0, 4);
A.d2d  = mean(d2d, 4);  A.d2d_sd  = std(d2d, 0, 4);
A.stats = mean(st, 4);  % NaN-aware mean is not needed: NaN only for the static model
A.stat_names = {'mean_speed','mean_sojourn','handovers_per_hour','mean_contact', ...
                'mean_intercontact','contact_frac','contacts_per_pair_hour', ...
                'p_stay_first','mean_reach'};
A.long_rate = mean(lr, 4);
% 95% confidence half-widths (Student t over seeds)
tq = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228];
A.t95 = tq(ns-1);
A.eff_ci = A.t95 * A.eff_sd / sqrt(ns); A.d2d_ci = A.t95 * A.d2d_sd / sqrt(ns);
% per-seed long table for external plotting
fid = fopen('results/raw.csv', 'w');
hdr = 'model,speed,seed,snap_mpc,snap_ma,snap_lru,eff_mpc,eff_ma,eff_lru';
for i = 1:nD, hdr = [hdr, sprintf(',d2d_Dd%d', P.Dd_list(i))]; end %#ok<AGROW>
for i = 1:numel(A.stat_names), hdr = [hdr, ',', A.stat_names{i}]; end %#ok<AGROW>
for i = 1:nD, hdr = [hdr, sprintf(',longrate_Dd%d', P.Dd_list(i))]; end %#ok<AGROW>
fprintf(fid, '%s\n', hdr);
for k = 1:size(raw, 1)
  fprintf(fid, '%s,%g,%d', P.models{raw(k,1)}, raw(k,2), raw(k,3));
  fprintf(fid, ',%.6g', raw(k,4:end)); fprintf(fid, '\n');
end
fclose(fid);
save('-v7', 'results/summary.mat', 'A');

% CSV: one row per (model, speed)
fid = fopen('results/summary.csv', 'w');
fprintf(fid, 'model,speed,snap_mpc,snap_ma,snap_lru,eff_mpc,eff_ma,eff_lru');
for i = 1:nD, fprintf(fid, ',d2d_Dd%d', P.Dd_list(i)); end
fprintf(fid, ',eff_ma_sd,d2d_sd_first,d2d_sd_last');
for i = 1:numel(A.stat_names), fprintf(fid, ',%s', A.stat_names{i}); end
fprintf(fid, '\n');
for mi = 1:nm
  for vi = 1:nv
    fprintf(fid, '%s,%g', P.models{mi}, P.speeds(vi));
    fprintf(fid, ',%.4f', squeeze(A.snap(mi,vi,:)));
    fprintf(fid, ',%.4f', squeeze(A.eff(mi,vi,:)));
    fprintf(fid, ',%.4f', squeeze(A.d2d(mi,vi,:)));
    fprintf(fid, ',%.4f,%.4f,%.4f', A.eff_sd(mi,vi,2), A.d2d_sd(mi,vi,1), A.d2d_sd(mi,vi,end));
    fprintf(fid, ',%.4f', squeeze(A.stats(mi,vi,:)));
    fprintf(fid, '\n');
  end
end
fclose(fid);
% contact-length histograms, one row per run: model, speed, seed, pairs, h_1..h_600
fid = fopen('results/contact_hist.csv', 'w');
for k = 1:size(hist, 1)
  fprintf(fid, '%s,%g,%d,%d', P.models{hist(k,1)}, hist(k,2), hist(k,3), hist(k,4));
  fprintf(fid, ',%d', hist(k,5:end)); fprintf(fid, '\n');
end
fclose(fid);
end
