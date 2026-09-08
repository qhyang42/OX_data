function q = OX_rsa_bh(p)
% Benjamini-Hochberg adjustment; input must contain the complete test family.
assert(all(isfinite(p(:))) && all(p(:)>=0 & p(:)<=1));
[v,order] = sort(p(:)); n = numel(v);
adjusted = min(1,flipud(cummin(flipud(v.*n./(1:n)'))));
q = zeros(size(p)); q(order) = adjusted;
end
