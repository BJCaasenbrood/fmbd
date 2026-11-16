function ns = tie(s, master, slave, opts)
%TIE Bond the surface of slave to the surface of master.
%   n = s.tie(master, slave) ties the slave surface nodes within Tol of the
%   master surface and returns their number. Use the finer mesh as slave.
%   Tied elastic bodies are reduced together. A tie to a rigid body makes
%   the tied nodes a coupling with a fixed joint.
%   Options: Tol, Stiffness = [kn kt] (elastic interface, per unit area),
%   Select (selector on the slave nodes).
    arguments
        s
        master
        slave
        opts.Tol double = []
        opts.Stiffness double = []
        opts.Select = []
    end
    [~, iMaster] = lookup(s, master, false);
    [~, iSlave] = lookup(s, slave, false);
    bm = s.Bodies{iMaster};
    bs = s.Bodies{iSlave};
    assert(~isempty(bm.F) && ~isempty(bs.F), 'fmbd:tie', 'A point mass has no surface; use a joint.');
    assert(~(bm.Rigid && bs.Rigid), 'fmbd:tie', 'Two rigid bodies cannot be tied; use a fixed joint.');
    if bs.Rigid % the rigid body is the master
        [bm, bs, iMaster, iSlave] = deal(bs, bm, iSlave, iMaster);
    end
    assert(~bm.Rigid || isempty(opts.Stiffness), 'fmbd:tie', ['An elastic interface to a ' ...
        'rigid body is not supported; use a coupling and a joint with Stiffness.']);

    [slaveNodes, masterNodes, N, normals, areas] = tieNodes(bm, bs, opts.Tol, opts.Select);
    ns = numel(slaveNodes);
    if ns == 0
        warning('fmbd:tie', 'No slave nodes within the tolerance.');
        return
    end

    if bm.Rigid % couple the tied nodes and fix the coupling to the rigid body
        c = s.couple(bs.Name, slaveNodes, Name = sprintf('%s_tie%d', bs.Name, numel(s.Refs) + 1));
        s.joint(bm.Name, c, 'fixed', Point = s.Refs(end).x, Name = sprintf('%s_tie', c));
    else
        s.Ties(end+1) = struct('master', iMaster, 'slave', iSlave, 'snodes', slaveNodes, ...
                               'mnodes', masterNodes, 'N', N, 'n', normals, 'A', areas, ...
                               'k', opts.Stiffness);
        s.built = false;
    end
end

function [slaveNodes, masterNodes, N, normals, areas] = tieNodes(bm, bs, tol, select)
% Surface nodes of the slave body bs within tol of the surface of the
% master body bm. Per slave node: the nodes of its master face, the shape
% function weights N at its projection, the face normal and its nodal area.
    Fm = bm.F;
    A = bm.X(Fm(:,1),:); % corners of the master faces
    B = bm.X(Fm(:,2),:);
    C = bm.X(Fm(:,3),:);

    % default tolerance: a tenth of the smaller median edge length
    if isempty(tol)
        edge = @(b, F) median(vecnorm(b.X(F(:,[2 3 1]),:) - b.X(F(:,1:3),:), 2, 2));
        tol = 0.1*min(edge(bm, Fm), edge(bs, bs.F));
    end

    % candidates: slave surface nodes inside the bounding box of the master
    slaveNodes = unique(bs.F(:));
    if ~isempty(select)
        slaveNodes = slaveNodes(select(bs.X(slaveNodes,:)));
    end
    lo = min(min(A, B), C) - tol; % bounding box of each master face
    hi = max(max(A, B), C) + tol;
    inBox = all(bs.X(slaveNodes,:) >= min(lo) & bs.X(slaveNodes,:) <= max(hi), 2);
    slaveNodes = slaveNodes(inBox);

    % closest master face of each candidate
    ns = numel(slaveNodes);
    distance = inf(ns, 1);
    face = zeros(ns, 1);
    bary = zeros(ns, 3); % barycentric coordinates of the projection
    for i = 1:ns
        p = bs.X(slaveNodes(i),:);
        near = find(all(p >= lo & p <= hi, 2));
        if isempty(near), continue, end
        [d, l] = closestOnTri(p, A(near,:), B(near,:), C(near,:));
        [distance(i), k] = min(d);
        face(i) = near(k);
        bary(i,:) = l(k,:);
    end

    % keep the nodes within the tolerance
    tied = distance <= tol;
    slaveNodes = slaveNodes(tied);
    face = face(tied);
    bary = bary(tied,:);
    masterNodes = Fm(face, :);
    N = triShape(bary, size(Fm, 2));
    normals = cross(B(face,:) - A(face,:), C(face,:) - A(face,:), 2);
    normals = normals./vecnorm(normals, 2, 2);

    % nodal areas of the interface, from the slave faces that touch a tied node
    % and face the master surface
    Fs = bs.F;
    faceNormal = cross(bs.X(Fs(:,2),:) - bs.X(Fs(:,1),:), bs.X(Fs(:,3),:) - bs.X(Fs(:,1),:), 2);
    nodeNormal = zeros(bs.nNode, 3);
    nodeNormal(slaveNodes,:) = normals;
    masterNormal = reshape(sum(reshape(nodeNormal(Fs,:), [size(Fs) 3]), 2), [], 3); % per slave face
    opposed = dot(faceNormal, masterNormal, 2) < ...
              -0.5*vecnorm(faceNormal, 2, 2).*vecnorm(masterNormal, 2, 2);
    Fs = Fs(opposed,:);
    faceArea = vecnorm(faceNormal(opposed,:), 2, 2)/2;
    weights = [1 1 1]/3; % share of the face area per node (Tet4 faces)
    if size(Fs, 2) == 6, weights = [3 3 3 16 16 16]/57; end % Tet10 faces
    nodeArea = accumarray(Fs(:), reshape(faceArea.*weights, [], 1), [bs.nNode 1]);
    areas = nodeArea(slaveNodes);
end

function [d, l] = closestOnTri(p, A, B, C)
% Distance from p to triangles (A, B, C) and the closest-point coordinates.

    e1 = B - A;
    e2 = C - A;
    v = p - A;
    a = dot(e1, e1, 2);
    b = dot(e1, e2, 2);
    c = dot(e2, e2, 2);
    r1 = dot(e1, v, 2);
    r2 = dot(e2, v, 2);
    s = (c.*r1 - b.*r2)./(a.*c - b.^2);
    t = (a.*r2 - b.*r1)./(a.*c - b.^2);
    l = [1 - s - t, s, t];

    out = any(l < 0, 2);
    if any(out) % closest point on an edge
        clamp = @(x) min(max(x, 0), 1);
        e3 = C(out,:) - B(out,:);
        u1 = clamp(r1(out)./a(out));
        u2 = clamp(r2(out)./c(out));
        u3 = clamp(dot(p - B(out,:), e3, 2)./dot(e3, e3, 2));
        z = zeros(nnz(out), 1);
        cand = cat(3, [1-u1, u1, z], [1-u2, z, u2], [z, 1-u3, u3]);
        q = @(k) cand(:,1,k).*A(out,:) + cand(:,2,k).*B(out,:) + cand(:,3,k).*C(out,:);
        dk = [vecnorm(p - q(1), 2, 2), vecnorm(p - q(2), 2, 2), vecnorm(p - q(3), 2, 2)];
        [~, k] = min(dk, [], 2);
        m = numel(k);
        l(out,:) = cand((1:m)' + (k - 1)*3*m + [0, m, 2*m]);
    end
    d = vecnorm(p - (l(:,1).*A + l(:,2).*B + l(:,3).*C), 2, 2);
end

function N = triShape(l, nf)
%TRISHAPE Tri3/Tri6 shape functions at barycentric points l (m x 3).
%   Tri6 edge nodes sit on edges (1,2) (2,3) (3,1).

    if nf == 3
        N = l;
    else
        N = [l.*(2*l - 1), 4*l(:,1).*l(:,2), 4*l(:,2).*l(:,3), 4*l(:,3).*l(:,1)];
    end
end
