function build(s)
% Compile the model (bodies, ties, joints, signals) into the data that the
% solvers use. Runs on first use, and again after the model has changed.
    if s.built, return, end
    s.kinCache = []; % joint kinematics of the old model (see jointKin)

    q = buildUnits(s);
    buildReferencePoints(s);
    buildJoints(s);

    s.transData = struct('terms', {}, 'coef', {}, 'k', {}, 'c', {});
    for T = s.Transmissions
        s.transData(end+1) = struct('terms', T.terms, 'coef', T.coef, 'k', T.k, 'c', T.damp);
    end

    s.nodeSets = struct('body', {}, 'u0', {}, 'Phi', {}, 'x0', {});
    s.inSignals = compileSignals(s, s.Inputs);
    s.outSignals = compileSignals(s, s.Outputs);

    % initial configuration: joints at their Position, elastic
    % transmissions unstretched, elastic bodies undeformed
    s.keepRows = [];
    q = initialConfiguration(s, q);

    % drop redundant constraints: keep the rows of G that a pivoted QR of
    % G' finds independent
    [~, G] = constraints(s, q, false);
    if ~isempty(G)
        [~, R, order] = qr(G', 'vector');
        nd = min(size(R));
        pivots = abs(R(sub2ind(size(R), 1:nd, 1:nd)));
        rankG = nnz(pivots > 1e-10*max(pivots));
        s.keepRows = sort(order(1:rankG));
    end

    s.q0 = q;
    s.nu0 = initialVelocities(s, q);
    s.built = true;
end

function q = buildUnits(s)
% Units: each rigid body, and each group of elastic bodies tied together.
% A unit gets its columns in nu and an anchor at its frame. Returns the
% reference configuration q.
    nb = numel(s.Bodies);
    group = tiedGroups(s);
    groups = unique(group, 'stable');
    nUnits = numel(groups);

    s.bodyUnit = zeros(1, nb);
    s.bodyOffset = zeros(1, nb);
    s.anchors = struct('body', {}, 'frame', {}, 'x0', {}, 'Psi', {});
    s.units = [];
    nCols = 0; % velocities assigned so far
    nModal = 0; % modal coordinates assigned so far
    for u = 1:nUnits
        members = find(group == groups(u));
        s.bodyUnit(members) = u;
        if isscalar(members) && s.Bodies{members}.Rigid
            B = rigidUnit(s.Bodies{members});
        else
            nodeOffset = cumsum([0, cellfun(@(b) b.nNode, s.Bodies(members))]);
            s.bodyOffset(members) = nodeOffset(1:end-1);
            B = elasticGroup(s, members);
        end

        % columns of the unit in nu = [v; omega; etaDot], and in q.e
        B.cols = nCols + (1:6 + B.ne);
        B.ecols = nModal + (1:B.ne);
        nCols = nCols + 6 + B.ne;
        nModal = nModal + B.ne;

        s.units = [s.units, B];
        s.anchors(u) = struct('body', u, 'frame', true, 'x0', zeros(3, 1), 'Psi', zeros(6, B.ne));
    end
    s.nv = nCols;

    q = struct('r', reshape([s.units.c0], 3, nUnits), ...
               'R', repmat(eye(3), 1, 1, nUnits), ...
               'e', zeros(nModal, 1));
end

function group = tiedGroups(s)
% Bodies that are tied together form a group; group(i) is the lowest body
% index in the group of body i.
    nb = numel(s.Bodies);
    parent = 1:nb;
    for T = s.Ties
        a = groupRoot(parent, T.master);
        b = groupRoot(parent, T.slave);
        parent(max(a, b)) = min(a, b);
    end
    group = zeros(1, nb);
    for i = 1:nb
        group(i) = groupRoot(parent, i);
    end
end

function i = groupRoot(parent, i)
    while parent(i) ~= i
        i = parent(i);
    end
end

function B = rigidUnit(b)
% Unit of a rigid body: frame at the centre of mass, no modes.
    [m, c, J] = b.massProperties();
    B = struct('elastic', false, 'm', m, 'J', J, 'c0', c(:), 'Xl', b.X - c(:)', ...
               'ne', 0, 'Phi', zeros(3*b.nNode, 0), 'Q', [], 's', [], 'CR', [], ...
               'Kc', [], 'Dc', [], 'Mee', [], 'Psi', zeros(6, 0, 0), 'cpl', [], ...
               'grounded', false, 'Q9', [], 'CRx', [], 'Qaz', [], 'Qju', [], 'Qjz', []);
end

function buildReferencePoints(s)
% Attach every reference point to an anchor. A rigid body uses the anchor
% of its unit. So does the first coupling of an elastic unit, which carries
% the unit frame. Every other coupling gets an anchor of its own.
    s.refFrames = struct('anchor', {}, 's', {}, 'x0', {});
    for k = 1:numel(s.Refs)
        u = s.bodyUnit(s.Refs(k).body);
        B = s.units(u);
        x = s.Refs(k).x(:);
        isCoupling = ~isempty(s.Refs(k).nodes);
        if ~isCoupling
            s.refFrames(k) = struct('anchor', u, 's', x - B.c0, 'x0', x);
        elseif ~B.grounded && B.cpl(1) == k
            s.refFrames(k) = struct('anchor', u, 's', zeros(3, 1), 'x0', x);
        else
            s.anchors(end+1) = struct('body', u, 'frame', false, 'x0', x - B.c0, ...
                                      'Psi', B.Psi(:,:,B.cpl == k));
            s.refFrames(k) = struct('anchor', numel(s.anchors), 's', zeros(3, 1), 'x0', x);
        end
    end
end

function buildJoints(s)
% Joints between anchors, with zero coordinates in the reference
% configuration. Supports become joints to the ground (reference points)
% or locks (joint coordinates).
    s.jointFrames = struct('A', {}, 'B', {}, 'a', {}, 'b', {}, 'QA', {}, 'QB', {}, ...
                           'free', {}, 'k', {}, 'c', {}, 'pairs', {});
    for J = s.Joints
        addJoint(s, J.A, J.B, J.p(:), J.R, J.free, J.k, J.damp);
    end

    s.locks = struct('joint', {}, 'dir', {}, 'value', {});
    for f = s.Fixes
        dirs = f.dirs;
        switch f.kind
            case 'nodes' % handled in the reduction of the elastic unit
            case 'ref'
                if isempty(dirs), dirs = 1:6; end
                free = ~ismember(1:6, dirs);
                addJoint(s, 0, f.idx, s.Refs(f.idx).x(:), eye(3), free, zeros(nnz(free), 1), 0);
            otherwise % joint coordinates, held at their Position
                J = s.Joints(f.idx);
                if isempty(dirs), dirs = 1:nnz(J.free); end
                for d = dirs(:)'
                    value = 0;
                    if ~isempty(J.position), value = J.position(d); end
                    s.locks(end+1) = struct('joint', f.idx, 'dir', d, 'value', value);
                end
        end
    end

    % an elastic unit with node supports has its frame fixed to the ground
    for u = find([s.units.grounded])
        s.jointFrames(end+1) = struct('A', 0, 'B', u, 'a', s.units(u).c0, 'b', zeros(3, 1), ...
                                      'QA', eye(3), 'QB', eye(3), 'free', false(1, 6), ...
                                      'k', zeros(6, 1), 'c', zeros(6, 1), ...
                                      'pairs', rotationPairs(false(1, 6)));
    end
end

function addJoint(s, iA, iB, p, R, free, k, c)
% Joint between reference points iA (0 = ground) and iB at point p with
% frame R (world axes, in the reference configuration, where the anchors
% are aligned with the world).
% joint point in the anchor of A and of B
    anchorA = 0;
    a = p;
    if iA > 0
        refA = s.refFrames(iA);
        anchorA = refA.anchor;
        a = p - refA.x0 + refA.s;
    end
    refB = s.refFrames(iB);
    anchorB = refB.anchor;
    b = p - refB.x0 + refB.s;
    assert(anchorA ~= anchorB, 'fmbd:joint', 'A joint connects a reference point to itself.');

    % stiffness and damping on the free coordinates
    nFree = nnz(free);
    stiffness = zeros(6, 1);
    damping = zeros(6, 1);
    stiffness(free) = k(:).*ones(nFree, 1);
    damping(free) = c(:).*ones(nFree, 1);

    s.jointFrames(end+1) = struct('A', anchorA, 'B', anchorB, 'a', a, 'b', b, 'QA', R, 'QB', R, ...
                                  'free', logical(free), 'k', stiffness, 'c', damping, ...
                                  'pairs', rotationPairs(free));
end

function pairs = rotationPairs(free)
% Axis pairs (i, j) of the rotation constraints e_i(A) . e_j(B) = 0 of a
% joint, Table 5.1.
    freeRot = free(4:6);
    switch nnz(freeRot)
        case 0 % the frames stay aligned
            pairs = [2 3; 3 1; 1 2];
        case 1 % free about k: axis k of B stays on axis k of A
            k = find(freeRot);
            i = setdiff(1:3, k);
            pairs = [i(:), [k; k]];
        case 2 % k locked: axis i of A stays normal to axis j of B, (i, j, k) cyclic
            k = find(~freeRot);
            pairs = circshift([1 2 3], -k);
            pairs = pairs(1:2);
        otherwise
            pairs = zeros(0, 2);
    end
end

function sig = compileSignals(s, list)
% Signals as (kind, idx, dir, type). A signal on the nodes of a body
% becomes a node set: the mean displacement of those nodes.
    sig = struct('kind', {list.kind}, 'idx', {list.idx}, 'dir', {list.dir}, 'type', {list.type});
    for c = find(strcmp({list.kind}, 'body'))
        i = list(c).idx; % body
        u = s.bodyUnit(i);
        B = s.units(u);
        nodes = s.bodyOffset(i) + list(c).nodes(:); % in the numbering of the unit

        % mean mode shape of the nodes
        Phi = zeros(3, B.ne);
        if B.ne > 0
            dofs = 3*nodes' - [2; 1; 0];
            Phi = reshape(mean(reshape(B.Phi(dofs(:), :), 3, numel(nodes), B.ne), 2), 3, B.ne);
        end

        s.nodeSets(end+1) = struct('body', u, 'u0', mean(B.Xl(nodes, :), 1)', 'Phi', Phi, ...
                                   'x0', mean(s.Bodies{i}.X(list(c).nodes, :), 1)');
        sig(c).kind = 'node';
        sig(c).idx = numel(s.nodeSets);
    end
end

function nu = initialVelocities(s, q)
% Velocities that satisfy the constraints and the given joint velocities;
% elastic bodies start without deformation rates.
% equations A*nu = rhs: constraints G*nu = 0, then the joint velocities
    [~, G] = constraints(s, q);
    A = G;
    rhs = zeros(size(G, 1), 1);
    for j = 1:numel(s.Joints)
        velocity = s.Joints(j).velocity;
        if isempty(velocity), continue, end
        [~, Jc] = coord(s, 'joint', j, q);
        A = [A; Jc]; %#ok<AGROW>
        rhs = [rhs; velocity(:)]; %#ok<AGROW>
    end
    if all(rhs == 0)
        nu = zeros(s.nv, 1);
        return
    end

    % no deformation rates: etaDot = 0
    nModal = numel([s.units.ecols]);
    E = zeros(nModal, s.nv);
    for i = 1:numel(s.units)
        B = s.units(i);
        E(B.ecols, B.cols(7:end)) = eye(B.ne);
    end
    A = [A; E];
    rhs = [rhs; zeros(nModal, 1)];

    nu = lsqminnorm(A, rhs);
    if norm(A*nu - rhs) > 1e-8*max(1, norm(rhs))
        warning('fmbd:System', 'Joint velocities are not consistent with the constraints; using the closest.');
    end
end
