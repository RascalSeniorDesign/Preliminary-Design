function [E] = keplarEquationE(M,e)
%The keplarEquationE function inputs the Mean Anamoly and the eccentricity
%of a spacecraft in a particular orbit and oupts the eccentric anamoly at
%said point.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      e         Eccentricity                  Scalar          Unitless
%      M         Mean Anamoly                  Scalar             deg
%      E         Eccentric Anamoly             Scalar             deg
%==========================================================================
%Initial Release, keplarEquationE.m, Tom Moline, 2/01/2014

%Begin Code

%==========================================================================
%                       Find Eccentric Anamoly
%==========================================================================
M=M*(pi/180);
if M>-pi || M<pi
    En=M-e; %Initial Value
else
    En=M+e;
end

i=1;
Nmax=10^3; %Iteration Limit

while i<Nmax
    Entemp=En+(M-En+e*sin(En))/(1-e*cos(En));%New E value
    if abs(Entemp-En)<10^6 %Check for convergence
        break
    else
        En=Entemp; %Set up for next itration
        i=i+1;
    end
end
    
    
    
    