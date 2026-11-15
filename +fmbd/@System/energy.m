function E = energy(s, x)
%ENERGY Energies along a simulation.
%   E = s.energy(x) returns kinetic, gravity, strain, springs and total.
    build(s); % ensure that the system is built
    nt = size(x.nu, 2);  
    z  = zeros(1, nt);
    E  = struct('kinetic', z, 'gravity', z, 'strain', z, 'springs', z, 'total', z);
    g  = s.Gravity(:);

    for k = 1:nt
        q = struct('r', x.r(:,:,k), 'R', x.R(:,:,:,k), 'e', x.e(:,k));
        nu = x.nu(:,k); % velocities
        
        % compute kinetic energy
        E.kinetic(k) = 0.5*nu'*massMatrix(s, q)*nu;

        % compute gravity and strain energies (potential energy) for each body
        for i = 1:numel(s.units)
            B = s.units(i);  S = zeros(3, 1);
            if B.elastic
                eta = q.e(B.ecols);  S = B.s*[1; eta];
                E.strain(k) = E.strain(k) + 0.5*eta'*B.Kc*eta;
            end
            E.gravity(k) = E.gravity(k) - g'*(B.m*q.r(:,i) + q.R(:,:,i)*S);
        end

        % compute spring energies for each joint and translational spring
        for j = 1:numel(s.jointFrames)
            J = s.jointFrames(j);
            if any(J.k), c = jointKin(s, j, q);  E.springs(k) = E.springs(k) + 0.5*J.k'*c.^2; end
        end

        % compute spring energies for each translational spring
        for T = s.transData(arrayfun(@(T) ~isempty(T.k), s.transData))
            E.springs(k) = E.springs(k) + 0.5*T.k*stretch(s, T, q)^2;
        end
    end

    % compute total energy
    E.total = E.kinetic + E.gravity + E.strain + E.springs;
end
