function [deltaVi,deltaVf,dr_,dv_]=CWSolverIdeal(drf_,dv0_,dr0_,t,Nt,rtgt0_)
%The CWSolverESCIdeal function 
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%   drf_         Final Relative Position     3x1 Vector           km
%   dv0_         Initial Realtive Velocity   3x1 Vector          km/s
%   t            Transfer Time                 Scalar             min
%   Nt           Number of Time Elements       Scalar          Unitless
%   rtgt0_       Target Spacecraft Orbit     3x1 Vector           km

%Initial Release, CWSolverSK.m, Tom Moline, 3/3/2014

%Begin Code

%% Initialize Variables
t=linspace(0.00001,t*60,Nt);
dr_=zeros(3,Nt);
dv_=zeros(3,Nt);

%% Find Initial and Final DeltaV's
[deltaVi,deltaVf] = CWPrussingDeltaVSolver(dr0_,dv0_,rtgt0_,drf_,t(Nt));
%% Find Relative Velocities and Positions Based on DeltaV's
for i=1:Nt
    [dr_(:,i),dv_(:,i)] = CWPrussingSolver(dr0_,dv0_,rtgt0_,...
                                           deltaVi,t(i));
end
end

