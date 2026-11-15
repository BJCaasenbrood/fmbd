function Ti = dexpInv(p)
%DEXPINV Inverse right Jacobian of SO(3), Eq. (1.4).
th = norm(p);  P = hat(p);
if th < 1e-6
    Ti = eye(3) + P/2 + (P*P)/12;
else
    Ti = eye(3) + P/2 + (1/th^2 - (1 + cos(th))/(2*th*sin(th)))*(P*P);
end
end
