function [kind, k] = lookup(s, id, whole)
% Resolve id (a name or a Body) to kind 'joint', 'ref' or 'body' with its
% index k. With whole, id is used as a whole (no node selection): joints and
% reference points then take precedence, and a rigid body resolves to its
% reference point.
    if isa(id, 'fmbd.Body')
        k = find(cellfun(@(b) b == id, s.Bodies));
        assert(isscalar(k), 'fmbd:System', '%s is not part of the system.', id.Name);
        kind = 'body';
        if whole && id.Rigid
            isCoupling = ~cellfun(@isempty, {s.Refs.nodes});
            k = find([s.Refs.body] == k & ~isCoupling);
            kind = 'ref';
        end
        return
    end

    if whole
        kind = 'joint';
        k = find(strcmp({s.Joints.Name}, id));
        if isscalar(k), return, end
        kind = 'ref';
        k = find(strcmp({s.Refs.Name}, id));
        if isscalar(k), return, end
    end

    kind = 'body';
    bodyNames = cellfun(@(b) b.Name, s.Bodies, 'UniformOutput', false);
    k = find(strcmp(bodyNames, id));
    assert(isscalar(k), 'fmbd:System', 'Unknown name %s.', string(id));
    if whole && s.Bodies{k}.Rigid
        [kind, k] = lookup(s, s.Bodies{k}, true);
    end
end
