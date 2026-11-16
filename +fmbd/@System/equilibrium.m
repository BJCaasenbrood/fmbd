function [x, info] = equilibrium(s, opts)
%EQUILIBRIUM Static equilibrium.
%   [x, info] = s.equilibrium() solves for the rest state under gravity,
%   springs and constant inputs. Free mechanisms keep their position. x
%   can be used for simulate, linearize, statespace and modes; info has
%   the fields converged, iterations, force and constraint.
%   Options: Input, Initial, Tol (default 1e-10), MaxIter (default 100).
    arguments
        s
        opts.Input double = []
        opts.Initial struct = struct([])
        opts.Tol (1,1) double {mustBePositive} = 1e-10
        opts.MaxIter (1,1) double = 100
    end
    build(s); % ensure that the system is built
    [q, ~] = startState(s, opts.Initial);
    n = s.nv;
    nu = zeros(n, 1);
    u = inputVector(s, opts.Input);
    L = max([1; abs(q.r(:))]); % model size, for the constraint tolerance
    F = forces(s, q, nu, u, false);  [g, G] = constraints(s, q);
    lambda = lsqminnorm(G', F);
    r = [F - G'*lambda; g]; % force balance, constraints
    err = @(r, F) [norm(r(1:n))/max(norm(F), 1), max(abs([r(n+1:end); 0]))/L];
    merit = @(e) e(1) + e(2);

    % solve for the equilibrium with a modified Newton method
    best = inf; stall = 0; conv = false;  Kt = []; age = 0;
    for it = 1:opts.MaxIter
        
        e = err(r, F); % compute the force and constraint residuals
        if all(e < opts.Tol), conv = true;  break, end % converged!
        m = max(e/opts.Tol); % progress measure

        if m < 0.5*best
            best = m;  
            stall = 0;
        elseif age > 1
            Kt = []; % slow with an old tangent: refresh it first
        else
            stall = stall + 1; % slow even with a fresh tangent
        end

        if stall >= 3 % plateau: round-off floor, or no equilibrium
            conv = all(e < 1e3*opts.Tol);
            break
        end

        % (re)-compute the tangent stiffness matrix if needed (improves convergence)
        if isempty(Kt) || age >= 10 % modified Newton: reuse the tangent
            Kt = frameStiffness(s, q, u, lambda);  age = 0;
        end

        % Newton step split by the constraints, Eq. (8.2): correct them, then
        % balance the forces in the directions they allow (null space of G)
        r1 = r(1:n);
        dc = -lsqminnorm(G, g);
        N  = null(G);  if isempty(G), N = eye(n); end
        z  = lsqminnorm(N'*Kt*N, N'*(r1 - Kt*dc)); % free mechanisms stay
        dq = dc + N*z;
        dLambda = lsqminnorm(G', r1 - Kt*dq);
        a = 1;  m0 = merit(e);

        for ls = 1:6 % backtracking line search
            qNew = retract(s, q, a*dq);  lambdaNew = lambda + a*dLambda;
            FNew = forces(s, qNew, nu, u, false);  [gNew, GNew] = constraints(s, qNew);
            rNew = [FNew - GNew'*lambdaNew; gNew];
            if merit(err(rNew, FNew)) < m0 || ls == 6, break, end
            a = a/2;
        end

        [q, lambda, F, G, g, r] = deal(qNew, lambdaNew, FNew, GNew, gNew, rNew);
        age = age + 1; % count iterations since the last tangent update
    end
    e = err(r, F); % compute the force and constraint residuals
    info = struct('converged', conv, 'iterations', it, 'force', e(1), 'constraint', e(2)*L);
    
    % no equilibrium found warning
    if ~conv
        warning('fmbd:equilibrium', ['No equilibrium found after %d iterations (force residual ' ...
            '%.1e, constraint residual %.1e); a free coordinate may carry a net force.'], ...
            it, e(1), e(2)*L);
    end

    % return the state at equilibrium
    x = state(s, q, nu, u);
end

function K = frameStiffness(s, q, u, lambda)
% Tangent stiffness for equilibrium: differences over the frame coordinates
% only, plus the modal stiffness Kc of the elastic units.
    frames = arrayfun(@(B) B.cols(1:6), s.units, 'UniformOutput', false);
    K = stiffness(s, q, u, lambda, [], [frames{:}]);
    for B = s.units(logical([s.units.elastic]))
        modal = B.cols(7:end);
        K(modal, modal) = K(modal, modal) + B.Kc;
    end
end
