function [X, Y] = gen_mobility(model, v, P, seed)
% GEN_MOBILITY  User trajectories sampled every second.
%   [X,Y] = GEN_MOBILITY(model, v, P, seed) returns U-by-T position matrices,
%   T = P.Tsim + P.Td + 1. Every model is run for P.Twarm seconds first and
%   that warm-up is discarded, so all models start close to their stationary
%   regime (Yoon et al., INFOCOM 2003). v is the travel speed while moving.

set_seed(seed);
T  = P.Tsim + P.Td + 1;
t0 = -P.Twarm;
tq = 0:T-1;
L  = P.L;
X  = zeros(P.U, T);  Y = zeros(P.U, T);

switch model
  case 'gaussmarkov'
    [X, Y] = gauss_markov(v, P, T);
    return;
end

for u = 1:P.U
  switch model
    case 'static'
      x = rand*L; y = rand*L;
      tw = [t0, T]; xw = [x, x]; yw = [y, y];

    case 'randwalk'
      [tw, xw, yw] = wp_randwalk(v, P, t0, T);

    case 'rwp'
      [tw, xw, yw] = wp_rwp(v, P, t0, T);

    case 'markovcell'
      [tw, xw, yw] = wp_markovcell(v, P, t0, T);

    case 'manhattan'
      [tw, xw, yw] = wp_manhattan(v, P, t0, T);

    case 'levy'
      [tw, xw, yw] = wp_levy(v, P, t0, T);

    otherwise
      error('unknown model %s', model);
  end
  X(u,:) = interp1(tw, xw, tq, 'linear');
  Y(u,:) = interp1(tw, yw, tq, 'linear');
end
end

% ---------------------------------------------------------------------------
function [tw, xw, yw] = wp_randwalk(v, P, t0, T)
L = P.L; tau = P.rw_epoch;
t = t0; x = rand*L; y = rand*L;
tw = t; xw = x; yw = y;
while t < T
  th = 2*pi*rand;
  [tk, xk, yk] = flight(x, y, v*tau*cos(th), v*tau*sin(th), tau, L);
  tw = [tw, t + tk]; xw = [xw, xk]; yw = [yw, yk]; %#ok<AGROW>
  t = t + tau; x = xk(end); y = yk(end);
end
end

function [tw, xw, yw] = wp_rwp(v, P, t0, T)
L = P.L;
t = t0; x = rand*L; y = rand*L;
tw = t; xw = x; yw = y;
while t < T
  xd = rand*L; yd = rand*L;
  s  = v*(0.5 + rand);                   % speed ~ U[0.5v, 1.5v], v_min > 0
  t  = t + max(hypot(xd-x, yd-y)/s, 1e-6);
  x = xd; y = yd;
  tw(end+1) = t; xw(end+1) = x; yw(end+1) = y; %#ok<AGROW>
  p = P.rwp_pause*rand;
  if p > 1e-6
    t = t + p; tw(end+1) = t; xw(end+1) = x; yw(end+1) = y; %#ok<AGROW>
  end
end
end

function [tw, xw, yw] = wp_markovcell(v, P, t0, T)
L = P.L; n = P.nSide; w = L/n;
ci = randi(n); cj = randi(n);
x = (ci-1+rand)*w; y = (cj-1+rand)*w; t = t0;
tw = t; xw = x; yw = y;
while t < T
  if rand > P.mc_stay
    nb = [ci-1 cj; ci+1 cj; ci cj-1; ci cj+1];
    ok = nb(:,1) >= 1 & nb(:,1) <= n & nb(:,2) >= 1 & nb(:,2) <= n;
    nb = nb(ok,:); k = randi(size(nb,1));
    ci = nb(k,1); cj = nb(k,2);
  end
  xd = (ci-1+rand)*w; yd = (cj-1+rand)*w;
  t = t + max(hypot(xd-x, yd-y)/v, 1e-6);
  x = xd; y = yd;
  tw(end+1) = t; xw(end+1) = x; yw(end+1) = y; %#ok<AGROW>
end
end

function [tw, xw, yw] = wp_manhattan(v, P, t0, T)
L = P.L; b = P.mh_block; m = round(L/b);   % intersections 0..m
gi = randi(m+1)-1; gj = randi(m+1)-1;
d  = randi(4);                              % 1:+x 2:-x 3:+y 4:-y
t = t0; tw = t; xw = gi*b; yw = gj*b;
dirs = [1 0; -1 0; 0 1; 0 -1];
left  = [3 4 2 1]; right = [4 3 1 2]; back = [2 1 4 3];
while t < T
  r = rand;
  if r < P.mh_straight, nd = d;
  elseif r < P.mh_straight + (1-P.mh_straight)/2, nd = left(d);
  else, nd = right(d);
  end
  % stay on the grid: if the chosen street leaves the area, pick a valid one
  cand = [nd, d, left(d), right(d), back(d)];
  for k = cand
    ni = gi + dirs(k,1); nj = gj + dirs(k,2);
    if ni >= 0 && ni <= m && nj >= 0 && nj <= m, nd = k; break; end
  end
  gi = gi + dirs(nd,1); gj = gj + dirs(nd,2); d = nd;
  t = t + b/v;
  tw(end+1) = t; xw(end+1) = gi*b; yw(end+1) = gj*b; %#ok<AGROW>
end
end

function [tw, xw, yw] = wp_levy(v, P, t0, T)
L = P.L;
t = t0; x = rand*L; y = rand*L;
tw = t; xw = x; yw = y;
while t < T
  l  = trunc_pareto(P.levy_alpha, P.levy_lmin, P.levy_lmax);
  th = 2*pi*rand;
  [tk, xk, yk] = flight(x, y, l*cos(th), l*sin(th), l/v, L);
  tw = [tw, t + tk]; xw = [xw, xk]; yw = [yw, yk]; %#ok<AGROW>
  t = t + l/v; x = xk(end); y = yk(end);
  p = trunc_pareto(P.levy_beta, P.levy_pmin, P.levy_pmax);
  t = t + p; tw(end+1) = t; xw(end+1) = x; yw(end+1) = y; %#ok<AGROW>
end
end

function [tk, xk, yk] = flight(x, y, dx, dy, dur, L)
% Straight flight of duration dur from (x,y) by (dx,dy), reflected at the
% walls. Waypoints every second (at least two) so that the linear
% interpolation follows the reflected path instead of cutting corners.
n = max(2, ceil(dur));
f = (1:n) / n;
[xk, yk] = reflect(x + f*dx, y + f*dy, L);
tk = f * dur;
end

function [X, Y] = gauss_markov(v, P, T)
% Gauss-Markov model as described by Camp et al. (2002): speed and direction
% are AR(1) processes; near an edge the mean direction points to the center,
% and users reflect off the walls.
L = P.L; U = P.U; a = P.gm_alpha; Tt = P.Twarm + T;
x = rand(U,1)*L; y = rand(U,1)*L;
s = v*ones(U,1); d = 2*pi*rand(U,1);
sig_s = 0.3*v; sig_d = pi/4;
X = zeros(U, T); Y = zeros(U, T);
for k = 1:Tt
  dm = d;                                    % mean direction defaults to current
  edge = x < P.gm_edge | x > L-P.gm_edge | y < P.gm_edge | y > L-P.gm_edge;
  dm(edge) = atan2(L/2 - y(edge), L/2 - x(edge));
  s = a*s + (1-a)*v + sqrt(1-a^2)*sig_s*randn(U,1);
  s = max(s, 0);
  d = a*d + (1-a)*dm + sqrt(1-a^2)*sig_d*randn(U,1);
  % reflect at the walls and mirror the heading. (Clamping the position to
  % the wall, used in the first version, pinned fast users against the edge
  % for several steps and crowded them into the border strip.)
  x = x + s.*cos(d); y = y + s.*sin(d);
  bx = x < 0 | x > L; by = y < 0 | y > L;
  d(bx) = pi - d(bx); d(by) = -d(by);
  [x, y] = reflect(x, y, L);
  if k > P.Twarm
    X(:, k-P.Twarm) = x; Y(:, k-P.Twarm) = y;
  end
end
end

% ---------------------------------------------------------------------------
function [x, y] = reflect(x, y, L)
x = mod(x, 2*L); x(x > L) = 2*L - x(x > L);
y = mod(y, 2*L); y(y > L) = 2*L - y(y > L);
end

function z = trunc_pareto(alpha, zmin, zmax)
% inverse-CDF sample from a Pareto(alpha) law truncated to [zmin, zmax]
u = rand;
a = zmin^(-alpha); b = zmax^(-alpha);
z = (a - u*(a - b))^(-1/alpha);
end
