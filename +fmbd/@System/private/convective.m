function Gdnu = convective(s, q, nu)
% dG/dt*nu, the velocity part of the constraints differentiated twice:
% G*nuDot + dG/dt*nu = 0. Central differences along nu.
    if norm(nu) == 0
        g = constraints(s, q);
        Gdnu = zeros(size(g));
        return
    end
    e = 1e-6/norm(nu);
    [~, Gp] = constraints(s, retract(s, q, e*nu));
    [~, Gm] = constraints(s, retract(s, q, -e*nu));
    Gdnu = (Gp - Gm)*nu/(2*e);
end
