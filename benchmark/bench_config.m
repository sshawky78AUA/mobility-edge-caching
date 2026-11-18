function P = bench_config()
% BENCH_CONFIG  Parameters of the unified mobility benchmark (Paper B, Sec. 6).
% All lengths in meters, times in seconds (one simulation slot = 1 s).

P.L        = 1000;          % side of the square area
P.U        = 200;           % number of users
P.Tsim     = 3600;          % simulated period in which requests arrive
P.Td       = 60;            % delivery deadline (window after a request)
P.Twarm    = 2000;          % warm-up discarded for every model

% Small cells: 5 x 5 grid of SBSs, nearest-SBS association
P.nSide    = 5;
P.nSBS     = P.nSide^2;

% Content and caches
P.N        = 1000;          % library size
P.gamma    = 0.8;           % Zipf exponent
P.C        = 50;            % SBS cache size (files)
P.c        = 10;            % user (D2D) cache size (files)

% Requests
P.lambda   = 1/120;         % requests per user per second

% Delivery requirements
P.Ds       = 10;            % seconds of uninterrupted SBS association to deliver a file
P.Dd_list  = [2 5 10 20 30];     % seconds of uninterrupted D2D contact to deliver a file
P.rD2D     = 50;            % D2D range

% Sweeps
P.speeds   = [0.5 1 2 5 10 20];   % travel speed while moving (m/s)
P.seeds    = 1:10;
P.models   = {'static','randwalk','rwp','gaussmarkov','markovcell','manhattan','levy'};
P.names    = {'Static','Random walk','RWP','Gauss-Markov','Markov cell','Manhattan','Levy walk'};

% Model-specific parameters
P.rw_epoch    = 30;          % random walk: seconds per straight flight
P.rwp_pause   = 60;          % RWP: pause ~ U[0, rwp_pause]
P.gm_alpha    = 0.85;        % Gauss-Markov memory
P.gm_edge     = 50;          % Gauss-Markov: steer to center within this distance of an edge
P.mc_stay     = 0.3;         % Markov cell: probability of staying in the same cell
P.mh_block    = 100;         % Manhattan: street spacing
P.mh_straight = 0.5;         % Manhattan: probability of going straight
P.levy_alpha  = 1.5;         % Levy: flight-length tail exponent
P.levy_beta   = 1.5;         % Levy: pause-time tail exponent
P.levy_lmin   = 10;  P.levy_lmax = 1000;
P.levy_pmin   = 10;  P.levy_pmax = 1800;

% Mobility-aware placement: subsampling step used to estimate the reach matrix
P.ma_step  = 10;
end
