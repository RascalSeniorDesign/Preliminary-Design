function [dr_,dv_] = CWPrussingSolver(dr0_,dv0_,rtgt0_,deltaVi,t)
%The CWPrussingRendezvousSolver takes in the initial relative position and
%velocities between a target and chaser spacecraft, the final relative
%position between each of them, the inertial target spacecraft orbital
%position, and a time, and outputs the deltaV's required to perform a
%maneuver in the provided time to the final position.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      dr0_      Initial Relative Positoin    3x1 Vector          km
%      dv0_      Initial Relative Velocity    3x1 Vector         km/s
%      rtgt0_    Initial Inertial Target Pos  3x1 Vector          km
%      t         Rendezvous Transfer Time      Scalar             s
%      drf_      Final Relative Position      3x1 CVector         km
%==========================================================================
%Initial Release, CWPrussingRendezvousSolver, Tom Moline, 3/16/2014

%Begin Code

%% Initialize Variables
% Convert Inputs to Canonical Units
dr0_=dr0_./6378.137;
deltaVi=deltaVi./7.9053838;
dv0_=dv0_./7.9053838;
%% Create State Transition Matrix Elements
[M,N,S,T] = CWStateMatrix(rtgt0_,t);
%% Find Relative Position and Velocity (Canonical Units)
dr_=M*dr0_+N*(deltaVi+dv0_);
dv_=S*dr0_+T*(deltaVi+dv0_);
%% Convert Results to SI Units (km, and km/s)
dr_=dr_*6378.137;
dv_=dv_*7.9053838;