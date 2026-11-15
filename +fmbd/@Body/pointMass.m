function b = pointMass(m, c, J, opts)
%POINTMASS Rigid body from mass properties.
%   b = fmbd.Body.pointMass(m, c, J) with J about c. It has no surface;
%   attach it with a joint. Options: Name, RefPoint (default c).
    arguments
        m (1,1) double
        c (1,3) double
        J (3,3) double = zeros(3)
        opts.Name (1,:) char = 'mass'
        opts.RefPoint double = []
    end

    % create a rigid body with no mesh (point mass)
    if isempty(opts.RefPoint) 
        opts.RefPoint = c;
    end

    % create a rigid body with no mesh (point mass)
    b = fmbd.Body(struct('X', zeros(0, 3), 'T', zeros(0, 4)), struct('rho', 0), ...
                 Name = opts.Name, Rigid = true, RefPoint = opts.RefPoint);
    b.MassProps = struct('m', m, 'c', c, 'J', (J + J')/2);
end
