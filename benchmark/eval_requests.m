function out = eval_requests(X, Y, S, R, Xmpc, Xma, Yd2d, P)
% EVAL_REQUESTS  Serve every request under each policy and return hit metrics.
%   Snapshot hit   : file cached at the SBS serving the user at request time.
%   Effective hit  : file delivered from an SBS cache, i.e. the user stays with
%                    some SBS caching the file for P.Ds consecutive seconds
%                    within P.Td seconds of the request (a handover interrupts
%                    the transfer).
%   D2D offloading : file in the user's own cache, or a user caching it stays
%                    within P.rD2D for Dd consecutive seconds within P.Td.
%   LRU            : reactive cache at the serving SBS only, updated on every
%                    request it receives.
K  = size(R, 1);
nD = numel(P.Dd_list);
snap = zeros(K, 3); eff = zeros(K, 3); d2d = zeros(K, nD);
lruF = zeros(P.nSBS, P.C); lruT = zeros(P.nSBS, P.C);
for k = 1:K
  u = R(k,1); t = R(k,2); f = R(k,3);
  win = t+1 : t+1+P.Td;
  s = S(u, win);
  [vals, lens] = run_lengths(s);
  long = lens >= P.Ds;
  % MPC (same content everywhere)
  snap(k,1) = Xmpc(s(1), f);
  eff(k,1)  = any(long & Xmpc(vals, f).');
  % mobility-aware
  snap(k,2) = Xma(s(1), f);
  eff(k,2)  = any(long & Xma(vals, f).');
  % LRU at the serving SBS
  b = s(1);
  j = find(lruF(b,:) == f, 1);
  if ~isempty(j)
    snap(k,3) = 1; eff(k,3) = long(1); lruT(b,j) = k;
  else
    [~, j] = min(lruT(b,:)); lruF(b,j) = f; lruT(b,j) = k;
  end
  % D2D
  if Yd2d(u, f)
    d2d(k,:) = 1;
  else
    J = find(Yd2d(:, f)); J(J == u) = [];
    if ~isempty(J)
      dx = X(J, win) - repmat(X(u, win), numel(J), 1);
      dy = Y(J, win) - repmat(Y(u, win), numel(J), 1);
      cnt = double(dx.^2 + dy.^2 <= P.rD2D^2);
      for i = 1:nD
        Dd = P.Dd_list(i);
        d2d(k,i) = any(any(conv2(cnt, ones(1, Dd), 'valid') >= Dd - 0.5));
      end
    end
  end
end
out.snap = mean(snap, 1);   % [MPC, MA, LRU]
out.eff  = mean(eff, 1);
out.d2d  = mean(d2d, 1);    % one entry per Dd
out.K    = K;
end
