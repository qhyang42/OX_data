function p = OX_rsa_empirical_p(observed, null)
% Two-sided inclusive percentile test, add-one correction in each tail.
assert(isscalar(observed) && isfinite(observed) && all(isfinite(null(:))) && ~isempty(null));
p = min(1,2*min(1+sum(null(:)<=observed),1+sum(null(:)>=observed))/(numel(null)+1));
end
