function [V, f] = modes(K, M, n)
%MODES Lowest natural frequencies and mass-normalised mode shapes.
%   [V, f] = fmbd.modes(K, M, n) with f in Hz. K may be singular (free
%   systems); rigid-body modes then have f close to 0.

n = min(n, size(K, 1)); % number of modes

if size(K, 1) <= 200 % small systems: use dense eig
    [V, L] = eig(full(K), full(M), 'vector');
else % large systems: use sparse eigs
    r = full(diag(K)./diag(M)); % shift below zero, scaled to the softest dofs
    r = r(r > 0 & isfinite(r)); % (keeps K - sigma*M well conditioned for sparse eigs)
    sigma = -1e-4;
    if ~isempty(r), sigma = -1e-4*min(r); end
    [V, L] = eigs(K, M, n, sigma);  L = diag(L);
end

[lam, i] = sort(real(L));
V   = real(V(:, i(1:n))); % (real) eigenvectors corresponding to the lowest eigenvalues
lam = lam(1:n);           % eigenvalues upto n modes
V = V./sqrt(dot(V, M*V)); % mass-normalise the modes
f = sqrt(max(lam, 0))/(2*pi);
end
