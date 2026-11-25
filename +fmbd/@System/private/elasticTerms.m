function [M, F, Ct] = elasticTerms(B, R, eta, nu, g)
% Mass matrix (Eq. 4.1) and generalised forces (Eqs. 4.2-4.5) of an
% elastic unit, from the inertia invariants of Eq. (3.4). With a third
% output also the velocity tangent Ct = -dF/dnu (Eq. 7.4).
% kron(a, b) is written reshape(b*a.', [], 1) and trace(W) as
% W(1) + W(5) + W(9): this runs in the inner loop, where the call overhead
% of kron and trace dominates.
ne = B.ne;
xi = [1; eta];

% first moment S and inertia J about the frame origin, Eq. (3.6)
S = B.s*xi;
S1 = B.s(:, 2:end);
Quu = reshape(B.Q9*reshape(xi*xi.', [], 1), 3, 3);
J = (Quu(1) + Quu(5) + Quu(9))*eye(3) - Quu;
J = (J + J')/2;

% mass matrix, Eq. (4.1)
MomegaEta = reshape(B.CRx*xi, 3, ne);
M = [B.m*eye(3),    -R*hat(S), R*S1
     (-R*hat(S))',  J,         MomegaEta
     (R*S1)',       MomegaEta', B.Mee];
M = (M + M')/2;
if nargout < 2, return, end

omega = nu(4:6);
etaDot = nu(7:end);

% velocity forces on the frame, Eqs. (4.2) and (4.3), with
% sum Ms u x (omega x y) = omega*tr(W) - W'*omega for W = sum Ms u*y'
Quz = reshape(B.Qaz*reshape(xi*etaDot.', [], 1), 3, 3);
W1 = -Quu*hat(omega); % centripetal
W2 = Quz; % Coriolis
trW1 = W1(1) + W1(5) + W1(9);
trW2 = W2(1) + W2(5) + W2(9);
a = S1*etaDot;
hv = R*(cross(omega, cross(omega, S)) + 2*cross(omega, a));
hOmega = omega*trW1 - W1'*omega + 2*(omega*trW2 - W2'*omega);

% velocity forces on the modes, Eq. (4.4); column j holds vec(Q_ju), vec(Q_jz)
Qju = reshape(B.Qju*xi, 9, ne);
Qjz = reshape(B.Qjz*etaDot, 9, ne);
trQju = Qju(1,:) + Qju(5,:) + Qju(9,:);
crQjz = Qjz([8 3 4], :) - Qjz([6 7 2], :); % Eq. (3.5)
hEta = (reshape(omega*omega.', 1, [])*Qju)' - (omega'*omega)*trQju' - 2*(omega'*crQjz)';

% gravity, elasticity and damping, Eq. (4.5)
gb = R'*g;
F = [B.m*g - hv
     cross(S, gb) - hOmega
     S1'*gb - hEta - B.Kc*eta - B.Dc*etaDot];
if nargout < 3, return, end

% velocity tangent Ct = d(hv; hOmega; hEta + Dc*etaDot)/d(v; omega; etaDot).
% The forces are quadratic in nu, so these derivatives are exact.
I3 = eye(3);
Ct = zeros(6 + ne);
e = 7:6 + ne;
Ct(1:3, 4:6) = R*((omega'*S)*I3 + omega*S' - 2*S*omega' - 2*hat(a));
Ct(1:3, e) = 2*R*hat(omega)*S1;

% hOmega: omega*tr(W1) = -omega*(c'*omega) with c = cr(Quu), and
% W1'*omega = omega x (Quu'*omega); W2 = sum_j etaDot_j Z_j with
% Z_j = Q_ju' (since Q(:,:,b,a) = Q(:,:,a,b)')
c = [Quu(2,3) - Quu(3,2); Quu(3,1) - Quu(1,3); Quu(1,2) - Quu(2,1)];
QuuTw = Quu'*omega;
Qj = reshape(Qju, 3, 3, ne);
QjW = reshape(pagemtimes(Qj, omega), 3, ne); % Q_ju*omega = Z_j'*omega
QjTW = reshape(pagemtimes(Qj, 'transpose', omega, 'none'), 3, ne); % Q_ju'*omega
Ct(4:6, 4:6) = -(c'*omega)*I3 - omega*c' + hat(QuuTw) - hat(omega)*Quu' ...
               + 2*(trW2*I3 - W2');
Ct(4:6, e) = 2*(omega*trQju - QjW);

% hEta_j = omega'*Q_ju*omega - |omega|^2 tr(Q_ju) - 2 omega'*cr(Q_jz)
Qjzk = reshape(B.Qjz, 9, ne, ne); % (vec, j, k): Q_{j+1,k+1}
crK = Qjzk([8 3 4], :, :) - Qjzk([6 7 2], :, :);
Ct(e, 4:6) = (QjW + QjTW)' - 2*trQju'*omega' - 2*crQjz';
Ct(e, e) = -2*reshape(omega'*reshape(crK, 3, []), ne, ne) + B.Dc;
end
