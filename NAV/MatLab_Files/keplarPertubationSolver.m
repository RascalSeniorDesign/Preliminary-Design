function [r_,v_] = keplarPertubationSolver(r0_,v0_,tf)
%The keplarPertubationSolver function takes in position and velcity vectors 
%(r0_, v0_) and, provided an end time (tf), finds the position and
%velocity vectors at some other point on the spacecraft orbit or
%trajectory. This works for all types of orbits/trajectories and makes use
%of universal variables and canonical time in order to speed up
%operations. This particular function also accounts for pertubations due to
%the Earth's oblateness.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      r0_       Starting positon vector      3-vector            km
%      v0_       Starting Velocity vector     3-vector           km/s
%      tf        Flight Time                   Scalar             min
%      r_        Ending Position Vector       3-vector            km
%      v_        Ending Velocity Vector       3-vector           km/s
%==========================================================================
%Initial Release, keplarPertubationSolver.m, Tom Moline, 2/01/2014

%Begin Code

%==========================================================================
%                  Find Initial Keps and Eccentric Anamoly
%==========================================================================
[p0,a0,e0,i0,Omega0,omega0,theta0] = keplarElements(r0_,v0_);
[E0] = thetaToAnomaly(e0,theta0);
M0=E0-e0*sind(E0);
a0=a0/6378.1;
p0=p0/6378.1;
n0=sqrt(1/a0^3);

%==========================================================================
%                      Update for Pertubations
%==========================================================================
J2=0.0010826267;
Omega=Omega0-(3*n0*J2)/(2*p0^2)*cosd(i0)*tf;
omega=omega0+(3*n0*J2)/(4*p0^2)*(4-5*(sind(i0))^2)*tf;
[E] = keplarEquationE(M0,e0);
[theta] = anamolyToTheta(e0,E);

%==========================================================================
%                Convert Back to Positon and Velocity
%==========================================================================
p0=p0*6378.1;

[r_,v_] = keplarElementsToRV(p0,e0,i0,Omega,omega,theta);






