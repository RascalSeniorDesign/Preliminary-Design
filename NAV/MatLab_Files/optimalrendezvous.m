function deltaV = optimalrendezvous(x0, t0, xf, tf, dVtime, n)
% deltaV = optimalRendezvous(x0, t0, xf, tf, dVtime, n) computes the 
% impulsive thrusts for optimal transfer from x0 (at time t0) to xf  
% (at time tf) given thrusts at the time vector dVtime.
%
% The state vector is [x y z xdot ydot zdot]'
%
% Optimal equations based on CHEN and NI, improved for matrix work and an
% arbitrary final position.
%
% optimalRendezvous.m   Created 2/28/2003   M. Swartwout
% optimalRendezvous.m   Updated 3/03/2003   M. Swartwout
numThrusts = length(dVtime);
deltaV = zeros(3, numThrusts);
phi = zeros(6, 6, numThrusts);

% Define parameters
coastArc = linearOrbitCoast(x0, tf, t0, n);
trans = [zeros(3,3) eye(3)];

% Compute G matrix
G = zeros(6,6);
for i = 1:numThrusts
    % Compute the transition matrix from this thrust time until the end
    % time
    phi(:, :, i) = linearOrbitMatrix(tf, dVtime(i), n);
    % This equation comes from the partial derivatives of the cost
    % function.  Yes, it's ugly.
    G = G + phi(:, :, i)*[zeros(3, 6); trans*phi(:, :, i)'];
end
% Compute Lagrange multiplier
lambda = inv(G)*(coastArc-xf);
% Computer each deltaV
for i=1:numThrusts
    deltaV(:,i) = -trans*phi(:, :, i)'*lambda;
end