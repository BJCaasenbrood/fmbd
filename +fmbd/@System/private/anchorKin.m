function [ra, Ra, Jv, Jw] = anchorKin(s, k, q)
% Position ra and rotation Ra of anchor k (0 = ground) with the
% Jacobians of its velocity (world) and angular velocity (anchor frame),
% Eqs. (3.3) and (5.1).
    Jv = zeros(3, s.nv);
    Jw = zeros(3, s.nv);
    if k == 0 % ground
        ra = zeros(3, 1);
        Ra = eye(3);
        return
    end

    A = s.anchors(k);
    B = s.units(A.body);
    c = B.cols;
    r = q.r(:, A.body);
    R = q.R(:,:,A.body);
    Jv(:, c(1:3)) = eye(3);

    if A.frame % the unit frame itself
        ra = r;
        Ra = R;
        Jw(:, c(4:6)) = eye(3);
        return
    end

    % coupling: moves with the deformation of the unit
    eta = q.e(B.ecols);
    PsiT = A.Psi(1:3, :);
    PsiR = A.Psi(4:6, :);
    ut = A.x0 + PsiT*eta; % position in the unit frame
    psi = PsiR*eta; % rotation vector relative to the unit frame
    E = expso3(psi);
    ra = r + R*ut;
    Ra = R*E;
    Jv(:, c(4:6)) = -R*hat(ut);
    Jv(:, c(7:end)) = R*PsiT;
    Jw(:, c(4:6)) = E';
    Jw(:, c(7:end)) = dexp(psi)*PsiR;
end

function T = dexp(p)
%DEXP Right Jacobian T(p) of SO(3), Eq. (1.3): omega = T*dp/dt.
    th = norm(p);
    P = hat(p);
    if th < 1e-6 % ensure accuracy for small angles
        T = eye(3) - P/2 + (P*P)/6;
    else
        T = eye(3) - (1 - cos(th))/th^2*P + (th - sin(th))/th^3*(P*P);
    end
end
