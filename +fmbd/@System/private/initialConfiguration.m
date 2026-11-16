function q = initialConfiguration(s, q)
% Move the frames until the joints fit, Eq. (8.1). The constraints and the
% given joint positions are solved exactly (hard); the other joint
% coordinates stay as close to zero as that allows (soft).
    % hard targets [joint, coordinate, value]; soft coordinates [joint, coordinate]
    target = zeros(0, 3);
    soft = zeros(0, 2);
    for j = 1:numel(s.Joints)
        position = s.Joints(j).position;
        nFree = nnz(s.Joints(j).free);
        if isempty(position)
            soft = [soft; j*ones(nFree, 1), (1:nFree)']; %#ok<AGROW>
        else
            target = [target; j*ones(nFree, 1), (1:nFree)', position(:)]; %#ok<AGROW>
        end
    end

    % elastic transmissions start unstretched, unless all their coordinates
    % are given joint positions
    isGiven = @(T) all(arrayfun(@(term) strcmp(term.kind, 'joint') && ...
                       ~isempty(s.Joints(term.idx).position), T.terms));
    unstretched = s.transData(arrayfun(@(T) ~isempty(T.k) && ~isGiven(T), s.transData));

    % only the frames move; the modal coordinates stay zero
    frames = arrayfun(@(B) B.cols(1:6), s.units, 'UniformOutput', false);
    frames = [frames{:}];

    % approach large positions in increments of at most 0.3
    nIncrements = max(1, ceil(max([abs(target(:,3)); 0])/0.3));
    for inc = 1:nIncrements
        fraction = inc/nIncrements;
        for it = 1:50
            % hard equations: constraints, joint positions, unstretched transmissions
            [g, G] = constraints(s, q, false);
            for i = 1:size(target, 1)
                [coords, Jc] = coord(s, 'joint', target(i,1), q);
                g = [g; coords(target(i,2)) - fraction*target(i,3)]; %#ok<AGROW>
                G = [G; Jc(target(i,2), :)]; %#ok<AGROW>
            end
            for T = unstretched
                [delta, Jd] = stretch(s, T, q);
                g = [g; delta]; %#ok<AGROW>
                G = [G; Jd]; %#ok<AGROW>
            end

            % soft equations: the other joint coordinates
            gs = zeros(size(soft, 1), 1);
            Gs = zeros(size(soft, 1), s.nv);
            for i = 1:size(soft, 1)
                [coords, Jc] = coord(s, 'joint', soft(i,1), q);
                gs(i) = coords(soft(i,2));
                Gs(i,:) = Jc(soft(i,2), :);
            end

            % Newton step: minimum norm for the hard equations, then the
            % soft ones in the null space of the hard ones
            Gh = G(:, frames);
            Js = Gs(:, frames);
            step = -lsqminnorm(Gh, g);
            if ~isempty(gs)
                Nh = null(full(Gh));
                if isempty(Gh), Nh = eye(numel(frames)); end
                step = step - Nh*lsqminnorm(Js*Nh, gs + Js*step);
            end
            if max(abs([g; 0])) < 1e-12 && norm(step) < 1e-12, break, end

            dq = zeros(s.nv, 1);
            dq(frames) = step;
            q = retract(s, q, dq);
        end
    end
    assert(isempty(g) || max(abs(g)) < 1e-9, 'fmbd:System', ['The initial configuration ' ...
        'cannot satisfy the joints, transmissions and joint positions (residual %.1e).'], max(abs(g)));
end
