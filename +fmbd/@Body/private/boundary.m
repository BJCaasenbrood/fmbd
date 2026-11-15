function F = boundary(T)
% Outward boundary faces: faces that belong to exactly one Tetrahedron. T is ne x 4 or ne x 10. F is nf x 3 or nf x 6.

if size(T, 2) == 4 % linear tetrahedra (Tet4)
    loc = [1 3 2; 1 2 4; 1 4 3; 2 3 4];
else % quadratic tetrahedra (Tet10)
    loc = [1 3 2 7 6 5; 1 2 4 5 10 8; 1 4 3 8 9 7; 2 3 4 6 9 10];
end

F = reshape(permute(reshape(T(:, loc'), [], size(loc, 2), 4), [1 3 2]), [], size(loc, 2));
[~, i, j] = unique(sort(F(:, 1:3), 2), 'rows');
F = F(i(accumarray(j, 1) == 1), :);
end
