function [F, Kt, Ct, M] = forces(s, q, nu, u, tangents)
% Generalised forces F (Eqs. 2.2, 4.5, 6.1) and the mass matrix M. With
% tangents (default true) also Kt = -dF/dq and Ct = -dF/dnu, Eq. (7.4).
    if nargin < 5, tangents = true; end
    n = s.nv;
    F = zeros(n, 1);
    Kt = zeros(n);
    Ct = zeros(n);
    M = zeros(n);
    g = s.Gravity(:);

    % inertia, gravity and elasticity of each unit
    for i = 1:numel(s.units)
        B = s.units(i);
        c = B.cols;
        R = q.R(:,:,i);
        if B.elastic
            eta = q.e(B.ecols);
            if tangents
                [M(c, c), F(c), Ct(c, c)] = elasticTerms(B, R, eta, nu(c), g);
                Kt(c(7:end), c(7:end)) = B.Kc;
            else
                [M(c, c), F(c)] = elasticTerms(B, R, eta, nu(c), g);
            end
        else
            omega = nu(c(4:6));
            M(c, c) = blkdiag(B.m*eye(3), B.J);
            F(c(1:3)) = B.m*g;
            F(c(4:6)) = -cross(omega, B.J*omega); % Eq. (2.2)
            Ct(c(4:6), c(4:6)) = hat(omega)*B.J - hat(B.J*omega);
        end
    end

    % springs and dampers on the joint coordinates, Eq. (6.1)
    for j = 1:numel(s.jointFrames)
        J = s.jointFrames(j);
        active = J.free & (J.k' ~= 0 | J.c' ~= 0);
        if ~any(active), continue, end
        [coords, Js] = jointKin(s, j, q);
        coords = coords(active);
        Js = Js(active, :);
        k = J.k(active);
        c = J.c(active);
        F = F - Js'*(k.*coords + c.*(Js*nu));
        Kt = Kt + Js'*(k.*Js);
        Ct = Ct + Js'*(c.*Js);
    end

    % elastic transmissions: spring and damper on the stretch
    isElastic = arrayfun(@(T) ~isempty(T.k), s.transData);
    for T = s.transData(isElastic)
        [delta, Jd] = stretch(s, T, q);
        F = F - Jd'*(T.k*delta + T.c*(Jd*nu));
        Kt = Kt + T.k*(Jd'*Jd);
        Ct = Ct + T.c*(Jd'*Jd);
    end

    % inputs
    for i = 1:numel(s.inSignals)
        sig = s.inSignals(i);
        isMoment = strcmp(sig.kind, 'ref') && sig.dir > 3;
        if isMoment % moment in world axes on a reference point
            [~, Ra, ~, Jw] = anchorKin(s, s.refFrames(sig.idx).anchor, q);
            mu = zeros(3, 1);
            mu(sig.dir - 3) = u(i);
            F = F + Jw'*(Ra'*mu);
        else % force or torque on a joint coordinate, point or node set
            [~, Ju] = coord(s, sig.kind, sig.idx, q);
            F = F + Ju(sig.dir, :)'*u(i);
        end
    end
end
