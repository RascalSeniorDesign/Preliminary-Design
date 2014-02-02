function [E,B,H] = thetaToAnomaly(e,theta)
%The thetaToAnomaly function takes in the eccentricity and true anamoly and
%finds the Eccentric, Parabolic, or Hyperbolic anamoly at that point.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      e         Eccentricity                  Scalar          Unitless
%      theta     True Anamoly                  Scalar             deg
%      E         Eccentric Anamoly             Scalar             deg
%      B         Parabolic Anamoly             Scalar             deg
%      H         Hyperbolic Anamoly            Scalar             deg
%==========================================================================
%Initial Release, thetaToAnamoly.m, Tom Moline, 2/01/2014

%Begin Code

%==========================================================================
%                       Find Eccentric Anamolies
%==========================================================================
if e<1.0
    E=asind((sind(theta)*sqrt(1-e^2))/(1+e*cosd(theta))); %Ellipse
elseif e==1.0
    B=tand(theta/2); %Parabola
else
    H=asinhd((sind(theta)*sqrt(e^2-1))/(1+e*cosd(theta))); %Hyperbola
end


