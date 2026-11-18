function spatial_check()
% SPATIAL_CHECK  Share of user-time near the edges and at the center of the
% area, for the Gauss-Markov, RWP and random walk models (seed 1). Used in
% Section 6 to explain the Gauss-Markov and RWP contact results.
P = bench_config();
for v = [2 20]
  for m = {'gaussmarkov', 'rwp', 'randwalk', 'levy', 'markovcell', 'manhattan'}
    [X, Y] = gen_mobility(m{1}, v, P, 1);
    x = X(:, 1:P.Tsim); y = Y(:, 1:P.Tsim);
    edge = mean(x(:) < 50 | x(:) > 950 | y(:) < 50 | y(:) > 950);
    cen  = mean(abs(x(:) - 500) < 150 & abs(y(:) - 500) < 150);
    fprintf('%-12s v=%2d  border strip 50 m: %.3f (area 0.19)  central 300 m square: %.3f (area 0.09)\n', ...
            m{1}, v, edge, cen);
  end
end
end
