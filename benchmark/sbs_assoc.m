function S = sbs_assoc(X, Y, P)
% SBS_ASSOC  Index of the nearest SBS for every user and second.
% SBSs sit at the centers of an nSide x nSide grid, so the nearest SBS is the
% grid cell that contains the user.
w  = P.L / P.nSide;
ci = min(floor(X / w) + 1, P.nSide);
cj = min(floor(Y / w) + 1, P.nSide);
S  = (cj - 1) * P.nSide + ci;
end
