function [f, V] = modes(s, x, n, opts)
%MODES Natural frequencies and mode shapes about a state at rest.
%   [f, V] = s.modes(x, n) returns n frequencies (Hz) and mode shapes.
%   Unstable directions have negative f. Option Input as in linearize.
    arguments
        s
        x struct = struct([])
        n (1,1) double = inf
        opts.Input double = []
    end
    lin = s.linearize(x, Input = opts.Input);
    K = (lin.K + lin.K')/2;

    % undamped modes K*P = M*P*lam, lowest first
    [P, lam] = eig(K, lin.M, 'vector');
    [lam, order] = sort(real(lam));
    lam(abs(lam) < 1e-10*max(abs(lam))) = 0; % rigid-body and mechanism modes
    n = min(n, numel(lam));

    f = sign(lam(1:n)).*sqrt(abs(lam(1:n)))/(2*pi);
    V = lin.T*real(P(:, order(1:n))); % mode shapes in the system velocities
end
