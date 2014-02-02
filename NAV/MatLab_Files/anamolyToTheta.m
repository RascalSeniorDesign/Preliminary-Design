function [theta] = anamolyToTheta(e,E)
%The anamolyToTheta function takes in the eccentricity and eccentric
%anamoly of a spacecraft and returns its true anamoly
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      e         Eccentricity                  Scalar          Unitless
%      theta     True Anamoly                  Scalar             deg
%      E         Eccentric Anamoly             Scalar             deg
%==========================================================================
%Initial Release, AnamolyToTheta.m, Tom Moline, 2/01/2014

%Begin Code

%==========================================================================
%                       Find Eccentric Anamolies
%==========================================================================
theta=asind((sind(E)*sqrt(1-e^2))/(1+e*cosd(E)));
