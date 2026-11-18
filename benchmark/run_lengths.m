function [vals, lens] = run_lengths(s)
% RUN_LENGTHS  Values and lengths of the maximal constant runs of a row vector.
s = s(:).';
starts = [1, find(s(2:end) ~= s(1:end-1)) + 1];
ends   = [starts(2:end) - 1, numel(s)];
vals   = s(starts);
lens   = ends - starts + 1;
end
