function p_value = mc_perm_test(x, y, perm)
% MC_PERM_TEST permutation test
%
% p-value is calculated using Monte Carlo simulations when perm is
% different from the total number of unique permutations.
%
arguments
    x (:,1) {mustBeNumeric,mustBeReal}
    y (:,1) {mustBeNumeric,mustBeReal,mustBeEqualSize(x,y)}
    perm (1,1) double = 100000 % number of permutations
end

[~, ~, ~, stats_initial] =  ttest(y, x);
initial = stats_initial.tstat;

xy = [x; y];
max_perm = factorial(length(xy));

try % try generating all possible permutations
    xy_perm = perms(xy);

    if size(xy_perm, 1) < perm
        perm = size(xy_perm, 1);
    end

    xy_perm = datasample(xy_perm,perm,1,'Replace',false);
catch

    xy_perm = num2cell(repmat(xy',[perm, 1]),2);
    xy_perm = cell2mat(cellfun(@(vec) datasample(vec, length(vec),'Replace',false), xy_perm, ...
        'UniformOutput', false));
end

locs_xy = [ones(size(x)); 2.*ones(size(x))];

perm_results = nan(1,perm);

for p = 1:perm

    [~, ~, ~, stats] = ttest(xy_perm(p,locs_xy==2)', xy_perm(p,locs_xy==1)');
    perm_results(p) = stats.tstat;

end

num = sum(abs(perm_results) >= abs(initial));

if perm < max_perm
    % in case of random permutations (i.e., montecarlo sample, and NOT full
    % permutation), the minimum p-value should not be 0, but 1/N
    num = num + 1;
    perm = perm + 1;
end

p_value = num/perm; % Monte-Carlo p-value

end

% Custom validation function
function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid,msg))
end
end