function [delta, Jd] = stretch(s, T, q)
% Transmission stretch delta = sum c_i s_i and its Jacobian, Section 5.3.
    delta = 0;
    Jd = zeros(1, s.nv);
    for i = 1:numel(T.terms)
        term = T.terms(i);
        [v, Jv] = coord(s, term.kind, term.idx, q);
        delta = delta + T.coef(i)*v(term.dir);
        Jd = Jd + T.coef(i)*Jv(term.dir, :);
    end
end
