function [g, G] = constraints(s, q, reduce)
% Constraints g(q) = 0 and their Jacobian G (Section 5). With reduce
% (default true) only the independent rows s.keepRows are returned.
    if nargin < 3, reduce = true; end
    g = zeros(0, 1);
    G = zeros(0, s.nv);

    % joints: locked translations and rotation constraints
    for j = 1:numel(s.jointFrames)
        [coords, Jc, gRot, GRot] = jointKin(s, j, q);
        locked = ~s.jointFrames(j).free(1:3);
        g = [g; coords(locked); gRot]; %#ok<AGROW>
        G = [G; Jc(locked, :); GRot]; %#ok<AGROW>
    end

    % rigid transmissions: stretch = 0
    isRigid = arrayfun(@(T) isempty(T.k), s.transData);
    for T = s.transData(isRigid)
        [delta, Jd] = stretch(s, T, q);
        g = [g; delta]; %#ok<AGROW>
        G = [G; Jd]; %#ok<AGROW>
    end

    % supported joint coordinates
    for L = s.locks
        [coords, Jc] = coord(s, 'joint', L.joint, q);
        g = [g; coords(L.dir) - L.value]; %#ok<AGROW>
        G = [G; Jc(L.dir, :)]; %#ok<AGROW>
    end

    if reduce
        g = g(s.keepRows);
        G = G(s.keepRows, :);
    end
end
