function [M, F] = elasticTerms(B, R, eta, nu, g)
% Mass matrix (Eq. 4.1) and generalised forces (Eqs. 4.2-4.5) of an
% elastic unit, from the inertia invariants of Eq. (3.4).
ne = B.ne;
xi = [1; eta];

% first moment S and inertia J about the frame origin, Eq. (3.6)
S = B.s*xi;
S1 = B.s(:, 2:end);
Quu = reshape(B.Q9*kron(xi, xi), 3, 3);
J = trace(Quu)*eye(3) - Quu;
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
Quz = reshape(B.Qaz*kron(etaDot, xi), 3, 3);
W1 = -Quu*hat(omega); % centripetal
W2 = Quz; % Coriolis
hv = R*(cross(omega, cross(omega, S)) + 2*cross(omega, S1*etaDot));
hOmega = omega*trace(W1) - W1'*omega + 2*(omega*trace(W2) - W2'*omega);

% velocity forces on the modes, Eq. (4.4); column j holds vec(Q_ju), vec(Q_jz)
Qju = reshape(B.Qju*xi, 9, ne);
Qjz = reshape(B.Qjz*etaDot, 9, ne);
trQju = Qju(1,:) + Qju(5,:) + Qju(9,:);
crQjz = Qjz([8 3 4], :) - Qjz([6 7 2], :); % Eq. (3.5)
hEta = (kron(omega, omega)'*Qju)' - (omega'*omega)*trQju' - 2*(omega'*crQjz)';

% gravity, elasticity and damping, Eq. (4.5)
gb = R'*g;
F = [B.m*g - hv
     cross(S, gb) - hOmega
     S1'*gb - hEta - B.Kc*eta - B.Dc*etaDot];
end
