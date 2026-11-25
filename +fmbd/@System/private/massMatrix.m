function M = massMatrix(s, q)
% Mass matrix M(q), Eqs. (2.2) and (4.1).
    M = zeros(s.nv);
    for i = 1:numel(s.units)
        B = s.units(i);
        c = B.cols;
        if B.elastic
            M(c, c) = elasticTerms(B, q.R(:,:,i), q.e(B.ecols));
        else
            M(c(1:3), c(1:3)) = B.m*eye(3); % M is zero elsewhere: no blkdiag
            M(c(4:6), c(4:6)) = B.J;
        end
    end
end
