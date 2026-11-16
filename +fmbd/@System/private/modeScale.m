function [x, a, V0] = modeScale(s, x, v, amplitude)
% State x (default: initial), the factor a that makes the largest node
% displacement of a*v equal to amplitude, and the node positions V0.
    build(s);
    nb = numel(s.Bodies);
    if isempty(x), x = state(s, s.q0, s.nu0); end
    V0 = arrayfun(@(i) vertices(s, i, x, size(x.nu, 2), 1), 1:nb, 'UniformOutput', false);

    % node motion per unit of v, from a small step along v
    e = 1e-6;
    x1 = s.perturb(x, e*v);
    motion = arrayfun(@(i) (vertices(s, i, x1, 1, 1) - V0{i})/e, 1:nb, 'UniformOutput', false);
    largest = max([cellfun(@(d) max(vecnorm(d, 2, 2)), motion), eps]);
    a = amplitude/largest;
end
