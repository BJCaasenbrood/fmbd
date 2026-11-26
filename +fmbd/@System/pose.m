function pose(s, varargin)
%POSE Set joint positions, e.g. before equilibrium, modes or simulate.
%   s.pose('elbow', pi/3) sets the Position of joint 'elbow' (one value per
%   joint coordinate); s.pose('shoulder', 0.5, 'elbow', pi/3) sets several.
%   [] clears a Position. The model is rebuilt in that pose on next use.
    assert(mod(numel(varargin), 2) == 0, 'fmbd:pose', 'Give joint names and positions in pairs.');
    for i = 1:2:numel(varargin)
        name = varargin{i};
        value = varargin{i+1};
        [kind, k] = lookup(s, name, true);
        assert(strcmp(kind, 'joint'), 'fmbd:pose', '%s is not a joint.', string(name));
        nFree = nnz(s.Joints(k).free);
        assert(isempty(value) || numel(value) == nFree, 'fmbd:pose', ...
            '%s: %d value(s) expected.', string(name), nFree);
        s.Joints(k).position = value(:);
    end
    s.built = false;
end
