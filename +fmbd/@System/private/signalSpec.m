function sig = signalSpec(s, name, id, sel, dirs, type)
% Signals of an input or output: one per direction, called name, or
% name_<dir> when there are several.
    used = [s.InputName, s.OutputName];
    assert(~any(strcmp(used, name) | startsWith(used, [name '_'])), 'fmbd:io', ...
        'Signal name %s is already used.', name);

    % what the signal acts on, and its default directions
    [kind, k] = lookup(s, id, isempty(sel));
    nodes = [];
    switch kind
        case 'body' % node set: all nodes unless selected
            if isempty(sel)
                nodes = (1:s.Bodies{k}.nNode)';
            else
                nodes = s.Bodies{k}.nodes(sel);
            end
            assert(~isempty(nodes), 'fmbd:io', '%s: no nodes selected.', name);
            if isempty(dirs), dirs = 1:3; end
        case 'joint'
            if isempty(dirs), dirs = 1:nnz(s.Joints(k).free); end
        otherwise % reference point
            if isempty(dirs), dirs = 1:6; end
    end

    sig = struct('Name', {}, 'kind', {}, 'idx', {}, 'dir', {}, 'nodes', {}, 'type', {});
    for d = dirs(:)'
        signalName = name;
        if numel(dirs) > 1, signalName = sprintf('%s_%d', name, d); end
        sig(end+1) = struct('Name', signalName, 'kind', kind, 'idx', k, 'dir', d, ...
                            'nodes', nodes, 'type', type); %#ok<AGROW>
    end
end
