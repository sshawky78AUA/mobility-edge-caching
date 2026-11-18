function run_rwp_warmup(v)
% RUN_RWP_WARMUP  Speed decay of classical RWP (speed ~ U(0, 2v], uniform start,
% no warm-up) versus the benchmark RWP (speed ~ U[0.5v, 1.5v], warm-up removed).
% Reports, per 10-minute bin, the mean user speed, the handover rate and the
% effective SBS hit ratio of most-popular caching, and the D2D offloading ratio
% for short (2 s) and long (30 s) transfers.
if nargin < 1, v = 5; end
P = bench_config();
P.Tsim = 4*3600; P.U = 200;
nb = P.Tsim/600;
p = zipf_pmf(P.N, P.gamma);
out = struct();
for variant = 1:2
  spd = zeros(numel(P.seeds), nb); ho = spd; hit = spd; d2s = spd; d2l = spd;
  for si = 1:numel(P.seeds)
    set_seed(300 + P.seeds(si));
    T = P.Tsim + P.Td + 1;
    X = zeros(P.U, T); Y = X;
    for u = 1:P.U
      if variant == 1, t = 0; else, t = -P.Twarm; end
      x = rand*P.L; y = rand*P.L; tw = t; xw = x; yw = y;
      while t < T
        xd = rand*P.L; yd = rand*P.L;
        if variant == 1, s = 2*v*rand + 1e-9; else, s = v*(0.5 + rand); end
        t = t + max(hypot(xd-x, yd-y)/s, 1e-6); x = xd; y = yd;
        tw(end+1) = t; xw(end+1) = x; yw(end+1) = y; %#ok<AGROW>
        pz = P.rwp_pause*rand;
        if pz > 1e-6, t = t + pz; tw(end+1) = t; xw(end+1) = x; yw(end+1) = y; end %#ok<AGROW>
      end
      X(u,:) = interp1(tw, xw, 0:T-1); Y(u,:) = interp1(tw, yw, 0:T-1);
    end
    S = sbs_assoc(X, Y, P);
    sp = hypot(diff(X,1,2), diff(Y,1,2));
    hc = double(S(:,2:end) ~= S(:,1:end-1));
    % user caches: c distinct files per user, drawn in proportion to popularity
    set_seed(20000 + P.seeds(si));
    Yd2d = false(P.U, P.N);
    for u = 1:P.U
      q = p; for k = 1:P.c
        f = find(cumsum(q)/sum(q) >= rand, 1, 'first'); Yd2d(u, f) = true; q(f) = 0;
      end
    end
    set_seed(500 + P.seeds(si));
    % requests: one per user every 30 s, file drawn from Zipf
    cdf = cumsum(p);
    for b = 1:nb
      cols = (b-1)*600 + (1:600);
      spd(si,b) = mean(reshape(sp(:,cols), [], 1));
      ho(si,b)  = sum(reshape(hc(:,cols), [], 1)) / P.U * 6;   % per user-hour
      hits = 0; n = 0; ds = 0; dl = 0;
      for u = 1:P.U
        for t = cols(1):30:cols(end)
          f = find(cdf >= rand, 1, 'first');
          [vals, lens] = run_lengths(S(u, t:t+P.Td)); %#ok<ASGLU>
          hits = hits + (f <= P.C) * any(lens >= P.Ds); n = n + 1;
          if Yd2d(u, f)
            ds = ds + 1; dl = dl + 1;
          else
            J = find(Yd2d(:, f)); J(J == u) = [];
            if ~isempty(J)
              win = t : t + P.Td;
              dx = X(J, win) - repmat(X(u, win), numel(J), 1);
              dy = Y(J, win) - repmat(Y(u, win), numel(J), 1);
              cnt = double(dx.^2 + dy.^2 <= P.rD2D^2);
              ds = ds + any(any(conv2(cnt, ones(1, 2), 'valid') >= 1.5));
              dl = dl + any(any(conv2(cnt, ones(1, 30), 'valid') >= 29.5));
            end
          end
        end
      end
      hit(si,b) = hits / n; d2s(si,b) = ds / n; d2l(si,b) = dl / n;
    end
  end
  out(variant).speed = spd; out(variant).ho = ho; out(variant).hit = hit;
  out(variant).d2d_short = d2s; out(variant).d2d_long = d2l;
end
save('-v7', sprintf('results/rwp_warmup_v%g.mat', v), 'out', 'P');
end
