function [p,a,e,i,Omega,omega,theta] = keplarSolver(r0_,v0_)
%The keplarSolver function finds the keplarian elements associated with the
%position and velocity vectors (r0_, v0_) at a particular point in time.
%These parameters include the semiparameter, semi-major axis, eccentricity,
%inclination, right ascension of the ascending node, argument of perigee,
%and true anamoly of a geocentric orbit, as defined in the Earth-Fixed
%coordinate system.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      r0_       Starting positon vector      3-vector            km
%      v0_       Starting Velocity vector     3-vector           km/s
%      p         Semiparameter                 Scalar             km
%      a         Semi-Major Axis Length        Scalar             km
%      e         Eccentricity                  Scalar          Unitless
%      i         Inclination                   Scalar             deg
%      Omega     RAAN                          Scalar             deg
%      omega     Argument of Perigee           Scalar             deg
%      theta     True Anamoly                  Scalar             deg
%==========================================================================
%Initial Release, keplarSolver.m, Tom Moline, 1/31/2014

%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
r0_=r0_./6378.1; %Converts from km to canonical units, DU
v0_=v0_./7.9053838; %Convets from km/s to DU/TU
r0=sqrt(sum(abs(r0_).^2));
v0=sqrt(sum(abs(v0_).^2));
%==========================================================================
%             Define Specific Momentum and Nodal Line Vectors
%==========================================================================
h_=cross(r0_,v0_); %Specific angular momentum vector, DU^2/TU
h=sqrt(sum(abs(h_).^2)); %Specific angular momentum magnitude

Khat=[0 0 1];
n_=cross(Khat,h_); %Nodal vector line
n=sqrt(sum(abs(n_).^2));

%==========================================================================
%                       Find Keplarian Elements
%==========================================================================
e_=((v0^2-(1/r0))*r0_-dot(r0_,v0_)*v0_); %Eccentricity vector
e=sqrt(sum(abs(e_).^2)); %Eccentricity
zi=(v0^2/2)-(1/r0); %Energy
if e ~=1.0
    a=-1/(2*zi); %Semi Marjor Axis, DU
    p=a*(1-e^2); %Semiparameter, DU
else
    p=h^2;
    a=inf;
end

a=a*6378.1; %Convert from canonical units to km
p=p*6378.1; %Ditto

i=acosd(h_(3)/h); %Inclination, deg
Omega=acosd(n_(1)/n); %RAAN, deg
if n_(2)<0
    Omega=360-Omega; %Correction for quadrent change
end
omega=acosd(dot(n_,e_)/(n*e)); %Argument of perigee, deg
if e_(3)<0
    omega=360-omega;
end
theta=acosd(dot(e_,r0_)/(r0*e)); %True Anomoly
if dot(r0_,v0_)<0
    theta=360-theta;
end



