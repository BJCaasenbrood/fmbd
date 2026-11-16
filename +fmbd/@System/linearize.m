function lin = linearize(s, x, opts)
%LINEARIZE Linear model about a state at rest.
%   lin = s.linearize(x) returns M, C, K, B, Cz, Cv, Ca and T of
%     M z'' + C z' + K z = B u,   y = Cz z + Cv z' + Ca z''
%   in independent coordinates z, with dq = T*z. K includes gravity,
%   constraint forces and preload. Default x: the initial configuration.
%   Option Input (default: x.u).
    arguments
        s
        x struct = struct([])
        opts.Input double = []
    end
    build(s);
    n = s.nv;
    nIn = numel(s.inSignals);
    nOut = numel(s.outSignals);

    % operating point: state and inputs
    [q, nu] = startState(s, x);
    if any(nu)
        warning('fmbd:linearize', 'The state is moving; the model is linearised at its velocities.');
    end
    u0 = opts.Input;
    if isempty(u0) && ~isempty(x) && isfield(x, 'u') && ~isempty(x.u)
        u0 = x.u(:, end);
    end
    u0 = inputVector(s, u0);

    % independent coordinates dq = T*z: null space of the constraints, Eq. (9.2)
    [~, G] = constraints(s, q);
    nc = size(G, 1);
    T = null(G);
    if isempty(G), T = eye(n); end

    % multipliers at the operating point, from Eq. (4.6). The dynamic solve
    % also works away from equilibrium, where the accelerations are not zero.
    [F, ~, ~, M] = forces(s, q, nu, u0);
    sol = lsqminnorm([M, G'; G, zeros(nc)], [F; -convective(s, q, nu)]);
    nuDot = sol(1:n);
    lambda = sol(n+1:end);
    if ~any(nu) && norm(M*nuDot) > 1e-6*max(norm(F), 1)
        warning('fmbd:linearize', ['The state is not an equilibrium (e.g. a preloaded spring ' ...
            'without a holding input); use x = s.equilibrium(Input = u0).']);
    end

    % stiffness, damping and input matrices, Eq. (9.1). The differences for
    % C and B use a step of 1 and are exact: the forces are quadratic in
    % nu and linear in u.
    K = stiffness(s, q, u0, lambda, nu);
    C = zeros(n);
    I = eye(n);
    for j = 1:n
        C(:, j) = -(forces(s, q, nu + I(:,j), u0) - forces(s, q, nu - I(:,j), u0))/2;
    end
    B = zeros(n, nIn);
    I = eye(nIn);
    for i = 1:nIn
        B(:, i) = forces(s, q, nu, u0 + I(:,i)) - forces(s, q, nu, u0);
    end

    % outputs y = Yq*dq + Yv*dq' + Ya*dq''
    Yq = zeros(nOut, n);
    Yv = zeros(nOut, n);
    Ya = zeros(nOut, n);
    for i = 1:nOut
        sig = s.outSignals(i);
        [~, Jc] = coord(s, sig.kind, sig.idx, q);
        switch sig.type
            case 'displacement', Yq(i,:) = Jc(sig.dir, :);
            case 'velocity',     Yv(i,:) = Jc(sig.dir, :);
            otherwise,           Ya(i,:) = Jc(sig.dir, :);
        end
    end

    % project onto the independent coordinates, Eq. (9.2)
    lin = struct('M', T'*M*T, 'C', T'*C*T, 'K', T'*K*T, 'B', T'*B, ...
                 'Cz', Yq*T, 'Cv', Yv*T, 'Ca', Ya*T, 'T', T, ...
                 'InputName', {s.InputName}, 'OutputName', {s.OutputName});
    lin.M = (lin.M + lin.M')/2;
end
