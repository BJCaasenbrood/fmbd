function M = massMatrix(s, q)
% Mass matrix M(q), Eqs. (2.2) and (4.1).
    M = zeros(s.nv);
    for i = 1:numel(s.units)
        B = s.units(i);
        c = B.cols;
        if B.elastic
            M(c, c) = elasticTerms(B, q.R(:,:,i), q.e(B.ecols));
        else
            M(c, c) = blkdiag(B.m*eye(3), B.J);
        end
    end
end
