function [coords, Jc, g, G] = jointKin(s, j, q)
% Joint coordinates [d; phi] (Eq. 5.2) with their Jacobian (Eq. 5.3), and
% the rotation constraints (Table 5.1) with their Jacobian (Eq. 5.4).
    J = s.jointFrames(j);
    [rA, RA, JvA, JwA] = anchorKin(s, J.A, q);
    [rB, RB, JvB, JwB] = anchorKin(s, J.B, q);

    % offset D = pB - pA between the joint points, and its Jacobian
    D = (rB + RB*J.b) - (rA + RA*J.a);
    JD = (JvB - RB*hat(J.b)*JwB) - (JvA - RA*hat(J.a)*JwA);

    % offset in the joint frame of A, and rotation vector of B relative to A
    Rrel = J.QA'*RA'*RB*J.QB;
    d = J.QA'*RA'*D;
    phi = logso3(Rrel);
    Jd = J.QA'*(hat(RA'*D)*JwA + RA'*JD);
    Jphi = dexpInv(phi)*(J.QB'*JwB - Rrel'*J.QA'*JwA);
    coords = [d; phi];
    Jc = [Jd; Jphi];

    % rotation constraints: axis i of A normal to axis j of B
    np = size(J.pairs, 1);
    g = zeros(np, 1);
    G = zeros(np, s.nv);
    for p = 1:np
        ua = J.QA(:, J.pairs(p, 1));
        wb = J.QB(:, J.pairs(p, 2));
        u = RA*ua;
        w = RB*wb;
        g(p) = u'*w;
        G(p, :) = -cross(RB'*u, wb)'*JwB - cross(RA'*w, ua)'*JwA;
    end
end
