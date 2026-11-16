function sys = statespace(s, x, opts)
%STATESPACE Modal state-space model.
%   sys = s.statespace(x) linearises about x and returns A, B, C, D in
%   modal coordinates [eta; deta], with signal names and frequencies f.
%   Options:
%     Modes             number of modes (default: all)
%     Damping           modal damping ratio, added to the physical dampers
%     StaticCorrection  keep the static gain when modes are cut (default true)
%     Input             inputs at the operating point (default: x.u)
%     Format            "struct" (default) or "ss" (Control System Toolbox)
    arguments
        s
        x struct = struct([])
        opts.Modes (1,1) double {mustBePositive} = inf
        opts.Damping double = 0
        opts.StaticCorrection (1,1) logical = true
        opts.Input double = []
        opts.Format (1,1) string {mustBeMember(opts.Format, ["struct", "ss"])} = "struct"
    end
    lin = s.linearize(x, Input = opts.Input);
    K = (lin.K + lin.K')/2;
    M = lin.M;
    B = lin.B;
    nz = size(K, 1);

    % undamped modes K*V = M*V*Omega^2, lowest first
    [V, lam] = eig(K, M, 'vector');
    [lam, order] = sort(real(lam));
    V = real(V(:, order));
    lam = zeroSmall(lam);

    % keep the first modes; add the static shapes of the inputs when modes are cut
    nModes = min(opts.Modes, nz);
    truncated = nModes < nz;
    if opts.StaticCorrection && truncated && ~isempty(B)
        [V, lam] = staticCorrection(K, M, B, V(:, 1:nModes), lam, nModes);
    else
        V = V(:, 1:nModes);
        lam = lam(1:nModes);
    end
    V = V./sqrt(dot(V, M*V)); % V'*M*V = I
    nr = size(V, 2);

    % modal equations, Eq. (9.3): eta'' + Z*eta' + Om2*eta = Bm*u
    zeta = opts.Damping(:);
    if isscalar(zeta), zeta = zeta*ones(nr, 1); end
    omega = sqrt(max(lam, 0));
    Z = diag(2*zeta.*omega) + V'*lin.C*V;
    Om2 = V'*lin.K*V;
    Bm = V'*B;

    % state space with the state [eta; eta']. Acceleration outputs use the
    % modal equations for eta'', which gives the feedthrough D.
    A = [zeros(nr), eye(nr); -Om2, -Z];
    Bs = [zeros(nr, size(B, 2)); Bm];
    Ca = lin.Ca*V;
    Cs = [lin.Cz*V - Ca*Om2, lin.Cv*V - Ca*Z];
    Ds = Ca*Bm;

    f = sign(lam).*sqrt(abs(lam))/(2*pi); % negative for unstable modes
    sys = struct('A', A, 'B', Bs, 'C', Cs, 'D', Ds, 'InputName', {lin.InputName}, ...
                 'OutputName', {lin.OutputName}, 'f', f);
    if opts.Format == "ss"
        sys = ss(A, Bs, Cs, Ds, 'InputName', lin.InputName, 'OutputName', lin.OutputName);
    end
end

function [V, lam] = staticCorrection(K, M, B, V, lamAll, nModes)
% Add the static shapes K\B to the kept modes V, so that the static gain is
% preserved, and diagonalise again in the enlarged basis.
    % shift K when it is singular (rigid-body or mechanism modes)
    sigma = 0;
    if any(lamAll(1:nModes) <= 0), sigma = -1e-6*max(abs(lamAll)); end
    X = (K - sigma*M)\B;

    % remove what V already contains, and M-orthonormalise the rest
    X = X - V*(V'*(M*X));
    G = X'*M*X;
    [P, mu] = eig((G + G')/2, 'vector');
    keep = mu > 1e-12*max([mu; eps]);
    W = [V, X*P(:, keep)./sqrt(mu(keep))'];

    Kw = W'*K*W;
    Mw = W'*M*W;
    [P, lam] = eig((Kw + Kw')/2, (Mw + Mw')/2, 'vector');
    [lam, order] = sort(real(lam));
    V = W*real(P(:, order));
    lam = zeroSmall(lam);
end

function lam = zeroSmall(lam)
% Round-off eigenvalues of rigid-body and mechanism modes: exactly zero.
    lam(abs(lam) < 1e-10*max(abs(lam))) = 0;
end
