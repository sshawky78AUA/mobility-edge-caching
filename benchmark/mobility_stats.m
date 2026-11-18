function st = mobility_stats(X, Y, S, P, seed)
% MOBILITY_STATS  Quantities that link a mobility model to the caching gain.
T = P.Tsim;
sp = hypot(diff(X(:,1:T), 1, 2), diff(Y(:,1:T), 1, 2));
st.mean_speed = mean(sp(:));                       % includes pauses
% SBS sojourn time (complete runs only) and handovers per user-hour
soj = []; nho = 0;
for u = 1:P.U
  [~, lens] = run_lengths(S(u, 1:T));
  nho = nho + numel(lens) - 1;
  if numel(lens) > 2, soj = [soj, lens(2:end-1)]; end %#ok<AGROW>
end
if isempty(soj), st.mean_sojourn = NaN; else, st.mean_sojourn = mean(soj); end
st.handovers_per_hour = nho / P.U * 3600 / T;
% D2D contacts on a random set of pairs
set_seed(50000 + seed);
npair = 3000;
a = randi(P.U, npair, 1); b = randi(P.U, npair, 1);
keep = a ~= b; a = a(keep); b = b(keep);
cd = []; icd = []; frac = 0;
for i = 1:numel(a)
  c = (X(a(i),1:T) - X(b(i),1:T)).^2 + (Y(a(i),1:T) - Y(b(i),1:T)).^2 <= P.rD2D^2;
  frac = frac + mean(c);
  [vals, lens] = run_lengths(double(c));
  if numel(lens) > 2
    inner = 2:numel(lens)-1;
    cd  = [cd,  lens(inner(vals(inner) == 1))]; %#ok<AGROW>
    icd = [icd, lens(inner(vals(inner) == 0))]; %#ok<AGROW>
  end
end
st.contact_frac = frac / numel(a);
st.mean_contact = mean_or_nan(cd);
st.mean_intercontact = mean_or_nan(icd);
st.n_contacts_per_pair_hour = numel(cd) / numel(a) * 3600 / T;
% contact-length histogram (complete contacts, lengths capped at 600 s), so
% the rate of contacts long enough for any transfer time can be computed
st.contact_hist = accumarray(min(cd(:), 600), 1, [600 1]).';
st.n_pairs = numel(a);
% SBS-side reach: for random (user, time) pairs, the probability that the
% user stays at least Ds in its current cell, and the number of distinct
% SBSs it stays with for at least Ds within the deadline window
set_seed(60000 + seed);
nw = 4000; uu = randi(P.U, nw, 1); tt = randi(T - P.Td - 1, nw, 1);
first = zeros(nw, 1); nreach = zeros(nw, 1);
for i = 1:nw
  [vals, lens] = run_lengths(S(uu(i), tt(i)+1 : tt(i)+1+P.Td));
  first(i) = lens(1) >= P.Ds;
  nreach(i) = numel(unique(vals(lens >= P.Ds)));
end
st.p_stay_first = mean(first);
st.mean_reach = mean(nreach);
end

function m = mean_or_nan(x)
if isempty(x), m = NaN; else, m = mean(x); end
end
