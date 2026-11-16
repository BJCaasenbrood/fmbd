function [Phi, X0] = modeShapes(s, body)
%MODESHAPES Nodal mode shapes of an elastic body.
%   [Phi, X0] = s.modeShapes(body): the nodes are at r + R*(X0 + Phi*eta).
    build(s);
    [kind, i] = lookup(s, body, false);
    B = s.units(s.bodyUnit(i));
    assert(strcmp(kind, 'body') && B.elastic, 'fmbd:modeShapes', 'Not an elastic body.');

    % rows of the body in its unit
    nodes = s.bodyOffset(i) + (1:s.Bodies{i}.nNode);
    dofs = 3*nodes - [2; 1; 0];
    Phi = B.Phi(dofs(:), :);
    X0 = B.Xl(nodes, :);
end
