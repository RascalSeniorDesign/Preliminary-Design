function x = linearOrbitCoast(x0, t, t0, n)
% function x = linearOrbitCoast(x0, t, t0, n) propogates the coasting
% linearized orbit from initial condition x0 (at time t0) across the time vector t.
%
% x         - (6 x m) linearized state vector [x y z xdot ydot zdot]'
% x0        - initial conditions (6 x 1)
% t         - vector of m time points
% t0        - time step where x=x0 (note that t0 > t(1) is allowed!)
% n         - angular rate of the circular reference orbit
%
% linearOrbitCoast.m Created 02/12/2003      M. Swartwout
% linearOrbitCoast.m Updated 02/28/2003      M. Swartwout

% Define parameters
m = length(t);
% Propogate
x = zeros(6, m);
for i=1:m
    % Propogate from original location
    x(:,i) = linearOrbitMatrix(t(i), t0, n)*x0;
end