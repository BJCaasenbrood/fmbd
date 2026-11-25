function B = elasticGroup(s, members)
% Collect the couplings, ties, supports and loads of an elastic group, in
% the node numbering of the group, and reduce it.
    bodies = s.Bodies(members);
    nn = sum(cellfun(@(b) b.nNode, bodies));
    offset = zeros(1, numel(s.Bodies)); % node offset of each member in the group
    offset(members) = s.bodyOffset(members);

    % couplings
    isCoupling = ~cellfun(@isempty, {s.Refs.nodes});
    cpl = find(ismember([s.Refs.body], members) & isCoupling);
    couplings = struct('Name', {}, 'nodes', {}, 'x', {});
    for k = cpl
        ref = s.Refs(k);
        couplings(end+1) = struct('Name', ref.Name, 'nodes', offset(ref.body) + ref.nodes(:), ...
                                  'x', ref.x); %#ok<AGROW>
    end

    % ties
    ties = s.Ties(ismember([s.Ties.master], members));
    for t = 1:numel(ties)
        ties(t).snodes = offset(ties(t).slave) + ties(t).snodes;
        ties(t).mnodes = offset(ties(t).master) + ties(t).mnodes;
    end

    % supported dofs
    supported = zeros(0, 1);
    isNodeFix = strcmp({s.Fixes.kind}, 'nodes') & ismember([s.Fixes.idx], members);
    for f = s.Fixes(isNodeFix)
        dofs = 3*(offset(f.idx) + f.nodes(:)') - 3 + f.dirs(:);
        supported = [supported; dofs(:)]; %#ok<AGROW>
    end

    % unit loads of the node-set inputs, for their static shapes
    loads = zeros(3*nn, 0);
    isNodeInput = strcmp({s.Inputs.kind}, 'body') & ismember([s.Inputs.idx], members);
    for in = s.Inputs(isNodeInput)
        f = zeros(3*nn, 1);
        f(3*(offset(in.idx) + in.nodes(:)) - 3 + in.dir) = 1/numel(in.nodes);
        loads = [loads, f]; %#ok<AGROW>
    end

    name = strjoin(cellfun(@(b) b.Name, bodies, 'UniformOutput', false), ' + ');
    nModes = sum(cellfun(@(b) b.Modes, bodies));
    % reuse a recent reduction with the same inputs: e.g. one System per
    % pose that differ only in a joint Position. The key holds copies of
    % the body data, so a Body edited in place is reduced again.
    % shortcut: keeps the last 4 reductions for the MATLAB session (clear
    % with "clear functions"); make it an option if memory becomes an issue
    persistent cache
    if isempty(cache), cache = struct('key', {}, 'unit', {}); end
    bodyData = cellfun(@(b) {b.X, b.T, b.Mat, b.Rigid, b.Damping}, bodies, 'UniformOutput', false);
    key = {bodyData, couplings, ties, unique(supported), loads, nModes, name};
    hit = find(arrayfun(@(c) isequal(c.key, key), cache), 1);
    if isempty(hit)
        B = reduceGroup(bodies, couplings, ties, unique(supported), loads, nModes, name);
        cache = [struct('key', {key}, 'unit', B), cache(1:min(end, 3))];
    else
        B = cache(hit).unit;
    end
    B.cpl = cpl;
end



function B = reduceGroup(bodies, cpl, ties, sup, loads, nModes, name)
% Reduce an elastic group (bodies tied together) to Craig-Bampton modes in a
% floating frame at its first coupling, and compute the inertia invariants.
% See Sections 3.1 and 3.3 of paper.md.
    X = cell2mat(cellfun(@(b) b.X, bodies(:), 'UniformOutput', false));
    nn = size(X, 1);
    nc = numel(cpl);
    grounded = ~isempty(sup);
    assert(grounded || nc > 0, 'fmbd:System', ['%s: an elastic body needs a coupling (s.couple), ' ...
           'which carries its frame, or node supports (s.fix).'], name);

    [K, M, D] = groupMatrices(bodies, ties, nn);

    % Step 1: dependent nodes, u = T*p. The first coupling is clamped and
    % carries the frame; the other couplings are interfaces. A grounded group is
    % clamped at its supports instead, so all its couplings are interfaces.
    iface = 1:nc;
    if ~grounded, iface = 2:nc; end
    [T, nFree] = dependentNodes(X, cpl, iface, ties, sup, name);
    Kr = T'*K*T;
    Mr = T'*M*T;
    Dr = T'*D*T;
    Kr = (Kr + Kr')/2;
    Mr = (Mr + Mr')/2;

    % Step 2: Craig-Bampton modes p = V*eta, Eq. (3.1), plus the static shapes
    % of the loads on node sets
    nq = size(T, 2);
    interior = (1:nFree)';
    boundary = (nFree + 1:nq)';
    [~, ~, V] = fmbd.craigBampton(Kr, Mr, boundary, nModes);
    if ~isempty(loads)
        V = [V, loadShapes(Kr, Mr, V, T'*loads, interior)];
    end
    ne = size(V, 2);
    Kc = V'*Kr*V;
    Mc = V'*Mr*V;
    Dc = V'*Dr*V;
    Phi = full(T*V); % mode shapes at the nodes

    % rows of V at each interface coupling: its displacement and rotation
    Psi = zeros(6, ne, nc);
    for i = 1:numel(iface)
        Psi(:,:,iface(i)) = V(nFree + 6*(i - 1) + (1:6), :);
    end

    % frame origin: the first coupling, or the centroid of a grounded group
    if nc > 0
        x0 = cpl(1).x(:);
    else
        x0 = mean(X, 1)';
    end
    Xl = X - x0';
    [Q, s, CR, Ms] = invariants(Xl, Phi, M);

    B = struct('elastic', true, 'm', full(sum(Ms(:))), 'J', [], 'c0', x0, 'Xl', Xl, 'ne', ne, ...
               'Phi', Phi, 'Q', Q, 's', s, 'CR', CR, 'Kc', (Kc + Kc')/2, 'Dc', (Dc + Dc')/2, ...
               'Mee', (Mc + Mc')/2, 'Psi', Psi, 'cpl', [], 'grounded', grounded);
    B = addReshapes(B);
end

function [K, M, D] = groupMatrices(bodies, ties, nn)
% Stiffness, mass and Rayleigh damping of the bodies side by side, with the
% springs of the elastic interfaces, K += G'*Ds*G.
    nb = numel(bodies);
    K = cell(1, nb);
    M = cell(1, nb);
    D = cell(1, nb);
    for i = 1:nb
        [K{i}, M{i}] = bodies{i}.matrices();
        D{i} = bodies{i}.Damping(1)*M{i} + bodies{i}.Damping(2)*K{i};
    end
    K = blkdiag(K{:});
    M = blkdiag(M{:});
    D = blkdiag(D{:});

    for t = ties(arrayfun(@(t) ~isempty(t.k), ties))
        ns = numel(t.snodes);
        kn = t.k(1); % normal stiffness per unit area
        kt = t.k(end); % tangential
        nnT = permute(t.n, [2 3 1]).*permute(t.n, [3 2 1]); % n*n' per slave node
        Db = reshape(t.A, 1, 1, []).*(kt*eye(3) + (kn - kt)*nnT);
        [p, q, i] = ndgrid(1:3, 1:3, 1:ns);
        Ds = sparse(3*(i(:) - 1) + p(:), 3*(i(:) - 1) + q(:), Db(:));
        G = gapOperator(t, nn);
        K = K + G'*Ds*G;
    end
end

function [T, nFree] = dependentNodes(X, cpl, iface, ties, sup, name)
% u = T*p: coupled nodes follow their coupling, tie slaves follow the master
% surface and supported dofs are zero. p holds the nFree free dofs, then six
% coordinates per interface coupling.
    nn = size(X, 1);
    nc = numel(cpl);

    % role of each node: coupled (owner > 0), tie slave, supported or free
    owner = zeros(nn, 1);
    for c = 1:nc
        if any(owner(cpl(c).nodes))
            error('fmbd:System', '%s: couplings %s and %s share nodes.', name, ...
                  cpl(max(owner(cpl(c).nodes))).Name, cpl(c).Name);
        end
        owner(cpl(c).nodes) = c;
    end
    rigidTies = ties(arrayfun(@(t) isempty(t.k), ties));
    slave = false(nn, 1);
    for t = rigidTies
        slave(t.snodes) = true;
    end
    assert(~any(slave & owner > 0), 'fmbd:System', '%s: a tied node is also in a coupling.', name);
    for t = rigidTies
        assert(~any(slave(t.mnodes(:))), 'fmbd:System', '%s: chained ties (a master node is tied).', name);
    end
    fixed = false(3*nn, 1);
    fixed(sup) = true;
    dependent = owner > 0 | slave;
    assert(~any(fixed(reshape(3*find(dependent)' - [2; 1; 0], [], 1))), 'fmbd:System', ...
           '%s: supported nodes cannot be coupled or tied.', name);

    % free dofs: identity
    freeDof = find(~fixed & ~reshape(repmat(dependent', 3, 1), [], 1));
    nFree = numel(freeDof);
    nq = nFree + 6*numel(iface);
    rows = freeDof;
    cols = (1:nFree)';
    vals = ones(nFree, 1);

    % coupled nodes: u = u_c + theta_c x (x - x_c)
    for i = 1:numel(iface)
        c = iface(i);
        nodes = cpl(c).nodes(:);
        np = numel(nodes);
        d = X(nodes, :) - cpl(c).x;
        first = nFree + 6*(i - 1);
        r = 3*nodes' - [2; 1; 0];
        H = zeros(3, 3, np); % -hat(d) per node
        H(1,2,:) = d(:,3);
        H(1,3,:) = -d(:,2);
        H(2,1,:) = -d(:,3);
        H(2,3,:) = d(:,1);
        H(3,1,:) = d(:,2);
        H(3,2,:) = -d(:,1);
        rows = [rows; r(:); reshape(repmat(r, 3, 1), [], 1)]; %#ok<AGROW>
        cols = [cols; repmat(first + (1:3)', np, 1); reshape(repmat(first + [4 5 6], 3, np), [], 1)]; %#ok<AGROW>
        vals = [vals; ones(3*np, 1); H(:)]; %#ok<AGROW>
    end
    T = sparse(rows, cols, vals, 3*nn, nq);

    % tie slaves follow the master surface
    if ~isempty(rigidTies)
        P = sparse(3*nn, 3*nn);
        for t = rigidTies
            P = P - gapOperator(t, nn, true);
        end
        T = T + P*T;
    end
end

function Xs = loadShapes(Kr, Mr, V, f, interior)
% Residual static shapes of the loads f: the static response with the
% interfaces held, made M-orthogonal to V and M-orthonormal.
    nq = size(Kr, 1);
    Xs = zeros(nq, size(f, 2));
    Xs(interior, :) = decomposition(Kr(interior, interior), 'chol')\full(f(interior, :));
    Xs = Xs - V*((V'*Mr*V)\(V'*Mr*Xs));
    G = full(Xs'*Mr*Xs);
    [P, lam] = eig((G + G')/2, 'vector');
    keep = lam > 1e-10*max([lam; eps]);
    Xs = Xs*P(:, keep)./sqrt(lam(keep))';
end

function [Q, s, CR, Ms] = invariants(Xl, Phi, M)
% Inertia invariants of the shapes U = [X0, Phi_1, ..., Phi_ne], Eq. (3.4):
% Q(:,:,a,b) = U_a'*Ms*U_b, s(:,a) = U_a'*Ms*1, and CR(:,a,j) = cr(Q_a,j).
    nn = size(Xl, 1);
    ne = size(Phi, 2);
    Ms = M(1:3:end, 1:3:end); % scalar mass matrix
    U = reshape(cat(3, Xl, permute(reshape(Phi, 3, nn, ne), [2 1 3])), nn, []);
    MU = Ms*U;
    Q = permute(reshape(full(U'*MU), 3, 1 + ne, 3, 1 + ne), [1 3 2 4]);
    s = reshape(full(sum(MU, 1)), 3, 1 + ne);
    CR = zeros(3, 1 + ne, ne);
    for a = 1:1 + ne
        for j = 1:ne
            CR(:, a, j) = cr(Q(:,:,a,j+1));
        end
    end
end

function G = gapOperator(t, nn, masterOnly)
% Sparse gap operator of a tie, G*u = u_slave - sum_j N_j u_master,j, with
% one row per slave dof. With masterOnly, only the master part, in rows at
% the slave dofs of the full numbering (3*nn x 3*nn).
    if nargin < 3, masterOnly = false; end
    xyz = reshape(1:3, 1, 1, 3); % the three dofs of a node, along dimension 3
    ns = numel(t.snodes);
    if masterOnly
        rows = 3*(t.snodes(:) - 1) + xyz + 0*t.mnodes;
        cols = 3*t.mnodes - 3 + xyz;
        vals = -t.N + 0*xyz;
        G = sparse(rows(:), cols(:), vals(:), 3*nn, 3*nn);
    else
        nodes = [t.snodes(:), t.mnodes];
        rows = 3*(0:ns - 1)' + xyz + 0*nodes;
        cols = 3*nodes - 3 + xyz;
        vals = [ones(ns, 1), -t.N] + 0*xyz;
        G = sparse(rows(:), cols(:), vals(:), 3*ns, 3*nn);
    end
end

function B = addReshapes(B)
% Reshaped inertia invariants for elasticTerms (computed once).
    ne = B.ne;
    B.Q9  = reshape(B.Q, 9, []); % Q_ab, (a,b) pairs
    B.CRx = reshape(permute(B.CR, [1 3 2]), 3*ne, 1 + ne);
    B.Qaz = reshape(B.Q(:,:,:,2:end), 9, []); % Q_{a,j+1}
    B.Qju = reshape(B.Q(:,:,2:end,:), 9*ne, 1 + ne); % Q_{j+1,b}
    B.Qjz = reshape(B.Q(:,:,2:end,2:end), 9*ne, ne); % Q_{j+1,k+1}
end

function v = cr(Q)
% sum_pq Ms_pq y_p x z_q from Q = sum_pq Ms_pq y_p z_q'.
    v = [Q(2,3) - Q(3,2); Q(3,1) - Q(1,3); Q(1,2) - Q(2,1)];
end
