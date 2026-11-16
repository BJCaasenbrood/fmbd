function name = joint(s, A, B, type, opts)
%JOINT Joint between two reference points.
%   name = s.joint(A, B, type) connects B to A ('ground' for A). The free
%   directions are joint coordinates, zero where the bodies are defined.
%   type: 'fixed' (default), 'revolute', 'prismatic', 'cylindrical',
%   'spherical', 'planar' or 'free'.
%   Options:
%     Point      joint location (default: reference point of B)
%     Axis       local z axis (default [0 0 1]), or Frame = [ex ey ez]
%     Free       mask [ux uy uz rx ry rz] in the joint frame, overrides type
%     Stiffness  spring on the joint coordinates, relaxed at zero
%     Damping    damper on the joint coordinates
%     Position   initial joint coordinates
%     Velocity   initial joint velocities
%     Name       default 'joint<k>'
    arguments
        s
        A
        B
        type (1,:) char {mustBeMember(type, {'fixed', 'revolute', 'prismatic', ...
            'cylindrical', 'spherical', 'planar', 'free'})} = 'fixed'
        opts.Point double = []
        opts.Axis (1,3) double = [0 0 1]
        opts.Frame double = []
        opts.Free = []
        opts.Stiffness double = 0
        opts.Damping double = 0
        opts.Position double = []
        opts.Velocity double = []
        opts.Name (1,:) char = ''
    end
    iA = refIndex(s, A);
    iB = refIndex(s, B);
    assert(iB > 0 && iA ~= iB, 'fmbd:joint', 'B must be a reference point other than A.');

    % free directions [ux uy uz rx ry rz] in the joint frame
    mask = struct('fixed', [0 0 0 0 0 0], 'revolute', [0 0 0 0 0 1], ...
        'prismatic', [0 0 1 0 0 0], 'cylindrical', [0 0 1 0 0 1], ...
        'spherical', [0 0 0 1 1 1], 'planar', [1 1 0 0 0 1], 'free', [1 1 1 1 1 1]);
    free = logical(mask.(type));
    if ~isempty(opts.Free), free = logical(opts.Free(:)'); end
    nFree = nnz(free);

    % joint frame [ex ey ez] and location
    R = opts.Frame;
    if isempty(R), R = frameFromAxis(opts.Axis); end
    p = opts.Point;
    if isempty(p), p = s.Refs(iB).x; end

    for option = {'Position', 'Velocity'}
        value = opts.(option{1});
        assert(isempty(value) || numel(value) == nFree, 'fmbd:joint', ...
            '%s: %d value(s) expected.', option{1}, nFree);
    end

    name = opts.Name;
    if isempty(name), name = sprintf('joint%d', numel(s.Joints) + 1); end
    assert(~ismember(name, names(s)), 'fmbd:joint', 'The name %s is already used.', name);

    s.Joints(end+1) = struct('Name', name, 'A', iA, 'B', iB, 'p', p(:)', 'R', R, 'free', free, ...
        'k', opts.Stiffness(:).*ones(nFree, 1), 'damp', opts.Damping(:).*ones(nFree, 1), ...
        'position', opts.Position(:), 'velocity', opts.Velocity(:));
    s.built = false;
end

function R = frameFromAxis(axis)
% Right-handed frame [ex ey ez] with ez along axis; ex is the world x axis
% (or y, when the axis is close to x) made normal to ez.
    ez = axis/norm(axis);
    helper = [1 0 0];
    if abs(ez(1)) > 0.9, helper = [0 1 0]; end
    ex = helper - dot(helper, ez)*ez;
    ex = ex/norm(ex);
    R = [ex; cross(ez, ex); ez]';
end

function k = refIndex(s, id)
% Reference point index of id; 0 for ground.
    k = 0;
    if (ischar(id) || isstring(id)) && strcmpi(id, 'ground'), return, end
    [kind, k] = lookup(s, id, true);
    assert(strcmp(kind, 'ref'), 'fmbd:joint', ...
        'Joint ends must be rigid bodies, couplings (s.couple) or ''ground''.');
end
