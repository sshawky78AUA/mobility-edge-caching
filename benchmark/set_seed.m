function set_seed(s)
% SET_SEED  Seed the generators in both MATLAB and GNU Octave.
if exist('OCTAVE_VERSION', 'builtin')
  rand('state', s); randn('state', s); %#ok<RAND>
else
  rng(s, 'twister');
end
end
