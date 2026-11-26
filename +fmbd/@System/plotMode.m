function h = plotMode(s, x, v, opts)
%PLOTMODE Draw a mode shape, coloured by the displacement of the nodes.
%   s.plotMode(x0, v) with x0 and v = V(:,k) from modes.
%   Options: Amplitude (largest node displacement, default 0.05), Title.
    arguments
        s
        x struct
        v (:,1) double
        opts.Amplitude (1,1) double {mustBePositive} = 0.05
        opts.Title = []
    end
    [x, a, V0] = modeScale(s, x, v, opts.Amplitude);
    nb = numel(s.Bodies);
    xm = s.perturb(x, a*v);
    V = arrayfun(@(i) vertices(s, i, xm, 1, 1), 1:nb, 'UniformOutput', false);

    % colour: displacement of each node from x0 (point masses stay grey)
    C = cell(1, nb);
    largest = eps;
    for i = 1:nb
        if isempty(s.Bodies{i}.F)
            C{i} = [0.5 0.5 0.5];
        else
            C{i} = vecnorm(V{i} - V0{i}, 2, 2);
            largest = max([largest; C{i}]);
        end
    end
    h = drawBodies(s.Bodies, V, C);
    colormap(gca, turbo);
    clim([0, largest]);
    cb = colorbar;
    cb.Label.String = 'displacement';

    if ~isempty(opts.Title), title(opts.Title); end
    if nargout == 0, clear h, end
end
