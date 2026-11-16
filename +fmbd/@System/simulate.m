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
    n = s.nv;
    nUnits = numel(s.units);
    nIn = numel(s.inSignals);
    maxIter = 30;

    % time grid
    nSteps = max(1, round(diff(tspan)/opts.Step));
    h = diff(tspan)/nSteps;
    t = tspan(1) + h*(0:nSteps)';

    % generalised-alpha parameters, Eq. (7.1)
    rho = opts.Rho;
    alphaM = (2*rho - 1)/(rho + 1);
    alphaF = rho/(rho + 1);
    gam = 0.5 + alphaF - alphaM;
    bet = 0.25*(gam + 0.5)^2;
    ca = (1 - alphaF)/(1 - alphaM);

    % start: accelerations and multipliers that satisfy Eq. (4.6)
    [q, nu] = startState(s, opts.Initial);
    u = zeros(nIn, 1);
    [nuDot, lambda] = startAccelerations();
    yk = outputs(s, q, nu, nuDot);
    if ~isempty(opts.Input) % the first input changes the accelerations
        u = opts.Input(t(1), yk);
        [nuDot, lambda] = startAccelerations();
        yk = outputs(s, q, nu, nuDot);
    end
    nc = numel(lambda);
    a = nuDot; % auxiliary accelerations, Eq. (7.2)

    % storage
    x = struct('r', zeros(3, nUnits, nSteps+1), 'R', zeros(3, 3, nUnits, nSteps+1), ...
               'e', zeros(numel(q.e), nSteps+1), 'nu', zeros(n, nSteps+1), ...
               'u', zeros(nIn, nSteps+1));
    Y = zeros(numel(s.outSignals), nSteps+1);
    record(1);

    for k = 1:nSteps
        if ~isempty(opts.Input)
            u = opts.Input(t(k), yk);
        end

        % Newton iterations on the new accelerations and multipliers
        nuDotNew = nuDot;
        previous = inf;
        for it = 1:maxIter
            [qNew, nuNew] = step(nuDotNew);
            if it == 1 % tangents once per step (modified Newton)
                [F, Kt, Ct, M] = forces(s, qNew, nuNew, u);
            else
                [F, ~, ~, M] = forces(s, qNew, nuNew, u, false);
            end
            [g, G] = constraints(s, qNew);

            residual = [M*nuDotNew - F + G'*lambda; g/(h^2*bet*ca)]; % Eq. (4.6)
            S = [M + h*gam*ca*Ct + h^2*bet*ca*Kt, G'; G, zeros(nc)]; % Eq. (7.4)
            dx = -S\residual;
            nuDotNew = nuDotNew + dx(1:n);
            lambda = lambda + dx(n+1:end);

            change = norm(dx(1:n))/(1 + norm(nuDotNew));
            converged = change <= opts.Tol;
            atRoundOff = change >= previous && change < 1e3*opts.Tol;
            if converged || atRoundOff, break, end
            previous = change;
        end
        if it == maxIter && change > 1e3*opts.Tol
            warning('fmbd:simulate', 'Newton did not converge at t = %g s; try a smaller Step.', t(k+1));
        end

        % accept the step
        [q, nu, a] = step(nuDotNew);
        nuDot = nuDotNew;

        % project the velocities onto the constraints, Eq. (7.5)
        if nc > 0
            [~, G] = constraints(s, q);
            M = massMatrix(s, q);
            sol = [M, G'; G, zeros(nc)]\[zeros(n, 1); -G*nu];
            nu = nu + sol(1:n);
        end

        yk = outputs(s, q, nu, nuDot);
        record(k + 1);
    end

    y = struct();
    for i = 1:numel(s.outSignals)
        y.(matlab.lang.makeValidName(s.Outputs(i).Name)) = Y(i,:)';
    end

    function [nuDot0, lambda0] = startAccelerations()
    % Accelerations and multipliers at the start, from Eq. (4.6) and the
    % constraints differentiated twice, G*nuDot + dG/dt*nu = 0.
        [F0, ~, ~, M0] = forces(s, q, nu, u);
        [~, G0] = constraints(s, q);
        m = size(G0, 1);
        sol0 = [M0, G0'; G0, zeros(m)]\[F0; -convective(s, q, nu)];
        nuDot0 = sol0(1:n);
        lambda0 = sol0(n+1:end);
    end

    function [qNew, nuNew, aNew] = step(nuDotNew)
    % Configuration and velocities at t(k+1) for the accelerations
    % nuDotNew, Eqs. (7.2) and (7.3).
        aNew = ((1 - alphaF)*nuDotNew + alphaF*nuDot - alphaM*a)/(1 - alphaM);
        qNew = retract(s, q, h*(nu + h*(0.5 - bet)*a + h*bet*aNew));
        nuNew = nu + h*(1 - gam)*a + h*gam*aNew;
    end

    function record(k)
    % Store the state and the outputs of step k.
        x.r(:,:,k) = q.r;
        x.R(:,:,:,k) = q.R;
        x.e(:,k) = q.e;
        x.nu(:,k) = nu;
        x.u(:,k) = u;
        Y(:,k) = cell2mat(struct2cell(yk));
    end
end

function y = outputs(s, q, nu, nuDot)
% Output signals as a struct, with fields in the order of OutputName.
% nuDot is needed for acceleration outputs.
    nOut = numel(s.outSignals);
    values = zeros(nOut, 1);
    for i = 1:nOut
        sig = s.outSignals(i);
        [coords, Jc] = coord(s, sig.kind, sig.idx, q);
        J = Jc(sig.dir, :);
        switch sig.type
            case 'displacement'
                values(i) = coords(sig.dir);
            case 'velocity'
                values(i) = J*nu;
            otherwise % acceleration: J*nuDot + dJ/dt*nu
                values(i) = J*nuDot;
                if any(nu) % dJ/dt by central differences along nu
                    e = 1e-6/norm(nu);
                    [~, Jp] = coord(s, sig.kind, sig.idx, retract(s, q, e*nu));
                    [~, Jm] = coord(s, sig.kind, sig.idx, retract(s, q, -e*nu));
                    values(i) = values(i) + (Jp(sig.dir, :) - Jm(sig.dir, :))*nu/(2*e);
                end
        end
    end
    y = cell2struct(num2cell(values), matlab.lang.makeValidName(s.OutputName(:)), 1);
    if nOut == 0, y = struct(); end
end
