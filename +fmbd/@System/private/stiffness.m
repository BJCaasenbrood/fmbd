function K = stiffness(s, q, u, lambda, nu, cols)
% Tangent stiffness -d(F - G'*lambda)/dq by central differences in the
% increments of Eq. (1.2), for the columns cols (default: all). Eq. (9.1).
    if nargin < 5 || isempty(nu), nu = zeros(s.nv, 1); end
    if nargin < 6, cols = 1:s.nv; end
    n = s.nv;
    K = zeros(n);
    e = 1e-5; % ~ eps^(1/3): balances truncation and round-off
    for j = cols(:)'
        dq = zeros(n, 1);
        dq(j) = e;
        rp = residual(retract(s, q, dq));
        rm = residual(retract(s, q, -dq));
        K(:, j) = -(rp - rm)/(2*e);
    end

    function r = residual(q)
        [~, G] = constraints(s, q);
        r = forces(s, q, nu, u, false) - G'*lambda;
    end
end
