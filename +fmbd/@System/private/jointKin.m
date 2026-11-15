function [c, Jc, g, Gr] = jointKin(s, j, q)
% Joint coordinates [d; phi] (Eq. 5.2) with their Jacobian (Eq. 5.3), and
% the rotation constraints (Table 5.1) with their Jacobian (Eq. 5.4).
    J = s.jointFrames(j);
    [rA, RA, JvA, JwA] = anchorKin(s, J.A, q);
    [rB, RB, JvB, JwB] = anchorKin(s, J.B, q);
    D = rB + RB*J.b - rA - RA*J.a;
    JD = (JvB - RB*hat(J.b)*JwB) - (JvA - RA*hat(J.a)*JwA);
    Rrel = J.QA'*RA'*RB*J.QB; % rotation of B relative to A
    phi = logso3(Rrel);
    c = [J.QA'*RA'*D; phi];
    Jc = [J.QA'*(hat(RA'*D)*JwA + RA'*JD); dexpInv(phi)*(J.QB'*JwB - Rrel'*J.QA'*JwA)];
    pairs = J.pairs;
    g = zeros(size(pairs, 1), 1);  Gr = zeros(size(pairs, 1), s.nv);
    for p = 1:size(pairs, 1)
        ua = J.QA(:, pairs(p,1));  wb = J.QB(:, pairs(p,2));
        u = RA*ua;  w = RB*wb;
        g(p) = u'*w;
        Gr(p, :) = -cross(RB'*u, wb)'*JwB - cross(RA'*w, ua)'*JwA;
    end
end
