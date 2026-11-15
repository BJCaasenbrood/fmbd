function [t, y, x] = simulate(s, tspan, opts)
%SIMULATE Simulate the system in time.
%   [t, y, x] = s.simulate([t0 t1]) starts in the initial configuration and
%   returns the time t, the outputs y (struct) and the states x, with
%   fields r, R (frame positions and rotations per unit), e (modal
%   coordinates), nu (velocities) and u (inputs).
%   Options:
%     Step     time step (default 1e-3)
%     Input    controller u = f(t, y), using the outputs of the previous step
%     Initial  start state, e.g. from equilibrium
%     Rho      numerical damping, 0..1 (default 0.8; 1 = none)
%     Tol      Newton tolerance (default 1e-8)

   arguments
        s
        tspan (1,2) double
        opts.Step (1,1) double {mustBePositive} = 1e-3
        opts.Input = []
        opts.Rho (1,1) double {mustBeInRange(opts.Rho, 0, 1)} = 0.8
        opts.Tol (1,1) double {mustBePositive} = 1e-8
        opts.Initial struct = struct([])
    end
    build(s);
    nt = max(1, round(diff(tspan)/opts.Step));
    h = diff(tspan)/nt;
    t = tspan(1) + h*(0:nt)';
    nb = numel(s.units);
    n = s.nv;

    % generalised-alpha parameters, Eq. (7.1)
    rho = opts.Rho;
    alphaM = (2*rho - 1)/(rho + 1);
    alphaF = rho/(rho + 1);
    gam = 0.5 + alphaF - alphaM;
    bet = 0.25*(gam + 0.5)^2;
    ca = (1 - alphaF)/(1 - alphaM);

    % consistent initial accelerations and multipliers from Eq. (4.6)
    [q, nu] = startState(s, opts.Initial);
    u = zeros(numel(s.inSignals), 1);
    [F, ~, ~, M] = forces(s, q, nu, u);
    [~, G] = constraints(s, q);
    nc = size(G, 1);
    sol = [M, G'; G, zeros(nc)]\[F; -convective(s, q, nu)];
    nuDot = sol(1:n);
    y0 = outputs(s, q, nu, nuDot);
    if ~isempty(opts.Input) % the first input changes the accelerations
        u = opts.Input(t(1), y0);
        F = forces(s, q, nu, u);
        sol = [M, G'; G, zeros(nc)]\[F; -convective(s, q, nu)];
        nuDot = sol(1:n);
        y0 = outputs(s, q, nu, nuDot);
    end
    lambda = sol(n+1:end);
    a = nuDot; % auxiliary accelerations, Eq. (7.2)

    x = struct('r', zeros(3, nb, nt+1), 'R', zeros(3, 3, nb, nt+1), ...
               'e', zeros(numel(q.e), nt+1), 'nu', zeros(n, nt+1), 'u', zeros(numel(s.inSignals), nt+1));
    Y = zeros(numel(s.outSignals), nt+1);
    x.r(:,:,1) = q.r;
    x.R(:,:,:,1) = q.R;
    x.e(:,1) = q.e;
    x.nu(:,1) = nu;
    x.u(:,1) = u;
    Y(:,1) = cell2mat(struct2cell(y0));
    yk = y0;
    for k = 1:nt
        if ~isempty(opts.Input), u = opts.Input(t(k), yk); end

        % Newton iterations on the new accelerations and multipliers
        nuDotNew = nuDot;
        prev = inf;
        for it = 1:30
            [qNew, nuNew] = step(nuDotNew);
            if it == 1 % tangents once per step (modified Newton)
                [F, Kt, Ct, M] = forces(s, qNew, nuNew, u);
            else
                [F, ~, ~, M] = forces(s, qNew, nuNew, u, false);
            end
            [g, G] = constraints(s, qNew);
            res = [M*nuDotNew - F + G'*lambda; g/(h^2*bet*ca)]; % Eq. (4.6)
            S = [M + h*gam*ca*Ct + h^2*bet*ca*Kt, G'; G, zeros(nc)]; % Eq. (7.4)
            dx = -S\res;
            nuDotNew = nuDotNew + dx(1:n);
            lambda = lambda + dx(n+1:end);
            e = norm(dx(1:n))/(1 + norm(nuDotNew));
            if e <= opts.Tol || (e >= prev && e < 1e3*opts.Tol), break, end % converged, or at round-off
            prev = e;
        end
        if it == 30 && e > 1e3*opts.Tol
            warning('fmbd:simulate', 'Newton did not converge at t = %g s; try a smaller Step.', t(k+1));
        end

        % update, and project the velocities onto the constraints, Eq. (7.5)
        [q, nu, aNew] = step(nuDotNew);
        if nc > 0
            [~, G] = constraints(s, q);
            M = massMatrix(s, q);
            p = [M, G'; G, zeros(nc)]\[zeros(n, 1); -G*nu];
            nu = nu + p(1:n);
        end
        a = aNew;
        nuDot = nuDotNew;
        yk = outputs(s, q, nu, nuDot);
        x.r(:,:,k+1) = q.r;
        x.R(:,:,:,k+1) = q.R;
        x.e(:,k+1) = q.e;
        x.nu(:,k+1) = nu;
        x.u(:,k+1) = u;
        Y(:,k+1) = cell2mat(struct2cell(yk));
    end
    y = struct();
    for i = 1:numel(s.outSignals)
        y.(matlab.lang.makeValidName(s.Outputs(i).Name)) = Y(i,:)';
    end

    function [qNew, nuNew, aNew] = step(nuDotNew)
    % Configuration and velocities at t(k+1) for the accelerations
    % nuDotNew, Eqs. (7.2) and (7.3).
        aNew = ((1 - alphaF)*nuDotNew + alphaF*nuDot - alphaM*a)/(1 - alphaM);
        qNew = retract(s, q, h*(nu + h*(0.5 - bet)*a + h*bet*aNew));
        nuNew = nu + h*(1 - gam)*a + h*gam*aNew;
    end
end
