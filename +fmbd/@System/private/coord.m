function [v, Jv] = coord(s, kind, idx, q)
% Coordinates and Jacobian of a joint, reference point or node set.
    switch kind
        case 'joint' % the free joint coordinates
            [coords, Jc] = jointKin(s, idx, q);
            free = s.jointFrames(idx).free;
            v = coords(free);
            Jv = Jc(free, :);

        case 'ref' % displacement and rotation vector from the reference configuration
            ref = s.refFrames(idx);
            [ra, Ra, Ja, Jw] = anchorKin(s, ref.anchor, q);
            phi = logso3(Ra);
            v = [ra + Ra*ref.s - ref.x0; phi];
            Jv = [Ja - Ra*hat(ref.s)*Jw; dexpInv(phi)*Jw];

        otherwise % node set: mean displacement of its nodes
            set = s.nodeSets(idx);
            B = s.units(set.body);
            c = B.cols;
            R = q.R(:,:,set.body);
            u = set.u0 + set.Phi*q.e(B.ecols);
            v = q.r(:, set.body) + R*u - set.x0;
            Jv = zeros(3, s.nv);
            Jv(:, c(1:3)) = eye(3);
            Jv(:, c(4:6)) = -R*hat(u);
            Jv(:, c(7:end)) = R*set.Phi;
    end
end
