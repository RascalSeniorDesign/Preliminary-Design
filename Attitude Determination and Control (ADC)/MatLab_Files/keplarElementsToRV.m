function [r_,v_] = keplarElementsToRV(p,e,i,Omega,omega,theta)
%The keplarElementsToRV function takes in all six keplerian elements
%(semi-parameter, eccentricity, inclination, right ascension of the
%ascending node, argument of perigee, and true anamoly) and returns the
%geocentric positon and velocity vectors of a given spacecraft.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      r_       Starting positon vector      3-vector             km
%      v_       Starting Velocity vector     3-vector            km/s
%      p         Semiparameter                 Scalar             km
%      e         Eccentricity                  Scalar          Unitless
%      i         Inclination                   Scalar             deg
%      Omega     RAAN                          Scalar             deg
%      omega     Argument of Perigee           Scalar             deg
%      theta     True Anamoly                  Scalar             deg
%==========================================================================
%Initial Release, keplarElementsToRV.m, Tom Moline, 2/01/2014

%Begin Code

%==========================================================================
%                      Convert to Canonical Units
%==========================================================================
p=p/6378.1; %Converts semi-parameter from km to DU

%==========================================================================
%              Define Perifocal Distance and Velocity Vectors
%==========================================================================
rper_=[(p*cosd(theta))/(1+e*cosd(theta))...
    (p*sind(theta))/(1+e*cosd(theta)) 0];

vper_=[-sqrt(1/p)*sind(theta) sqrt(1/p)*(e+cosd(theta)) 0];

%==========================================================================
%                      Convert to Geocentric Frame
%==========================================================================
transform=[cosd(Omega)*cosd(omega)-sind(Omega)*sind(omega)*cosd(i),...%1,1
     -cosd(Omega)*sind(omega)-sind(Omega)*cosd(omega)*cosd(i),...%(1,2)
     sind(Omega)*sind(i);...%(1,3)
     sind(Omega)*cosd(omega)+cosd(Omega)*sind(omega)*cosd(i),...%(2,1)
     -sind(Omega)*sind(omega)+cosd(Omega)*cosd(omega)*cosd(i),...%(2,2)
     -cosd(Omega)*sind(i);...%(2,3)
     sind(omega)*sind(i), cosd(omega)*sind(i), cosd(i)];%(3,1)-(3,3)
 
 r_=transpose(transform*transpose(rper_));%Transform perifocal vectors
 v_=transpose(transform*transpose(vper_));
 
 r_=r_.*6378.1; %Convert from canonical to km
 v_=v_.*7.9053838; %Convert from canonical to km/s










