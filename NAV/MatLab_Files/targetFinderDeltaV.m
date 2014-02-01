%==========================================================================
%                 targetFinderDeltaV Script Description
%==========================================================================
%The targetFinderDeltaV script defines the initial positon and velocity of
%an intercepter and target statellite. From this, it calculates the
%deltaV required for various transfer orbits such that the two spacecraft
%end up in the same position at the end of said orbits. It then plots the
%total deltaV required for such a maneuver and the target/interceptor
%orbits

%Initial Release, targetFinderDeltaV.m, Tom Moline, 1/31/2014

%Begin Code

%==========================================================================
%             Define Target/Interceptor Positions/Velocities
%==========================================================================

vtgt_=[-1.962372 7.323674 0.000000]; %Target velocity, km/s
vint_=[-4.864779 5.816486 .240163]; %Interceptor velocity, km/s
rtgt_=[6697.4756 1794.5831 0.000]; %Target Position, km
rint_=[5328.7862 4436.1273 101.4720]; %Interceptor Position, km

%==========================================================================
%             Define Transfer Times and Pre-Allocate for Speed
%==========================================================================
tf=linspace(5,250,10000); %Transfer time, minutes
rtgtx=zeros(1,length(tf)); %X, Y, and Z components of target position
rtgty=zeros(1,length(tf));
rtgtz=zeros(1,length(tf));
rintx=zeros(1,length(tf));%X, Y, and Z components of interceptor position
rinty=zeros(1,length(tf));
rintz=zeros(1,length(tf));
deltaVax=zeros(1,length(tf));%X, Y, and Z components of initial deltaV
deltaVay=zeros(1,length(tf));
deltaVaz=zeros(1,length(tf));
deltaVbx=zeros(1,length(tf));%X, Y, and Z components of final deltaV
deltaVby=zeros(1,length(tf));
deltaVbz=zeros(1,length(tf));

%==========================================================================
%        Find Orbit Positions/Velocities and Calculate Transfer DeltaV
%==========================================================================
for i=1:length(tf)
    [deltaVatemp_,deltaVbtemp_] = targetFinder(rint_,rtgt_,vint_,vtgt_,tf(i));
    [rtgttemp_,vtgttemp_] = keplarSolver(rtgt_,vtgt_,tf(i));
    [rinttemp_,vinttemp_] = keplarSolver(rint_,vint_,tf(i));
    rtgtx(i)=rtgttemp_(1);
    rtgty(i)=rtgttemp_(2);
    rtgtz(i)=rtgttemp_(3);
    rintx(i)=rinttemp_(1);
    rinty(i)=rinttemp_(2);
    rintz(i)=rinttemp_(3);
    deltaVax(i)=deltaVatemp_(1);
    deltaVay(i)=deltaVatemp_(2);
    deltaVaz(i)=deltaVatemp_(3);
    deltaVbx(i)=deltaVbtemp_(1);
    deltaVby(i)=deltaVbtemp_(2);
    deltaVbz(i)=deltaVbtemp_(3);
end

%==========================================================================
%                        Find DeltaV Magnitudes
%==========================================================================

deltaVamag=(abs(deltaVax).^2+abs(deltaVay).^2+abs(deltaVaz).^2).^.5;
deltaVbmag=(abs(deltaVbx).^2+abs(deltaVby).^2+abs(deltaVbz).^2).^.5;
deltaVtotal=deltaVbmag+deltaVamag;

%==========================================================================
%                             Plot Results
%==========================================================================

figure(1)
plot(tf,deltaVamag,tf,deltaVtotal)
xlabel('Time (mins)')
ylabel('Delta V (km/s)')

figure(2)
plot3(rtgtx,rtgty,rtgtz,rintx,rinty,rintz)
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')




