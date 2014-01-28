function phi = linearOrbitMatrix(t, t0, n)
% phi = linearOrbitMatrix(t, t0, n) returns the state transition matrix
% from time t0 to t, given orbital rate n under linearized gravity
% approximations.
%
% phi       - 6 x 6 transition matrix
% t         - vector of m time points
% t0        - time step where x=x0 (note that t0 > t(1) is allowed!)
% n         - angular rate of the circular reference orbit
%
% linearOrbitMatrix.m Created 02/28/2003      M. Swartwout

c = cos(n*(t-t0));
s = sin(n*(t-t0));

phi = [(4-3*c) 0 0  s/n 2*(1-c)/n 0;
    6*(s-n*(t-t0)) 1 0 -2*(1-c)/n (4*s-3*n*(t-t0))/n 0;
    0 0 c 0 0 s/n;
    3*n*s 0 0 c 2*s 0;
    -6*n*(1-c) 0 0 -2*s 4*c-3 0;
    0 0 -n*s 0 0 c];