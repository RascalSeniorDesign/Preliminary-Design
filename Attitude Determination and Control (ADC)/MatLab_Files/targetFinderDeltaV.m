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

clear
clc

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
tf=linspace(0,250,200); %Transfer time, minutes
td=linspace(0,250,40);
rtgtx=zeros(1,length(tf)); %X, Y, and Z components of target position
rtgty=zeros(1,length(tf));
rtgtz=zeros(1,length(tf));
rintx=zeros(1,length(tf));%X, Y, and Z components of interceptor position
rinty=zeros(1,length(tf));
rintz=zeros(1,length(tf));
deltaVax=zeros(length(td),length(tf));%X, Y, and Z components of initial deltaV
deltaVay=zeros(length(td),length(tf));
deltaVaz=zeros(length(td),length(tf));
deltaVbx=zeros(length(td),length(tf));%X, Y, and Z components of final deltaV
deltaVby=zeros(length(td),length(tf));
deltaVbz=zeros(length(td),length(tf));
rtransx=zeros(1,length(tf));
rtransy=zeros(1,length(tf));
rtransz=zeros(1,length(tf));
[deltaVatemp1_,deltaVbtemp2_] = targetFinder(rint_,rtgt_,vint_,vtgt_,tf(200),-1);

for i=1:length(tf)
    [rtgttemp_,vtgttemp_] = keplarSolver(rtgt_,vtgt_,tf(i));
    [rinttemp_,vinttemp_] = keplarSolver(rint_,vint_,tf(i));
    [rtranstemp_,vtranstemp_] = keplarSolver(rint_,(vint_+deltaVatemp1_),tf(i));
    rtgtx(i)=rtgttemp_(1); %Fill position/deltaV vector components
    rtgty(i)=rtgttemp_(2);
    rtgtz(i)=rtgttemp_(3);
    rintx(i)=rinttemp_(1);
    rinty(i)=rinttemp_(2);
    rintz(i)=rinttemp_(3);
    rtransx(i)=rtranstemp_(1);
    rtransy(i)=rtranstemp_(2);
    rtransz(i)=rtranstemp_(3);
end

%==========================================================================
%        Find Orbit Positions/Velocities and Calculate Transfer DeltaV
%==========================================================================
for j=1:length(td)
    [rtgt_,vtgt_] = keplarSolver(rtgt_,vtgt_,td(j));
    [rint_,vint_] = keplarSolver(rint_,vint_,td(j));
    for i=1:length(tf)
        tm=1;
        [deltaVatempshort_,deltaVbtempshort_] = targetFinder(rint_,rtgt_,vint_,vtgt_,tf(i),tm);
        tm=-1;
        [deltaVatemplong_,deltaVbtemplong_] = targetFinder(rint_,rtgt_,vint_,vtgt_,tf(i),tm);
        deltaVatempshortmag=sqrt(sum(abs(deltaVatempshort_).^2));
        deltaVatemplongmag=sqrt(sum(abs(deltaVatemplong_).^2));
        deltaVbtempshortmag=sqrt(sum(abs(deltaVbtempshort_).^2));
        deltaVbtemplongmag=sqrt(sum(abs(deltaVbtemplong_).^2));
        if deltaVatempshortmag<deltaVatemplongmag
            deltaVax(j,i)=deltaVatempshort_(1);
            deltaVay(j,i)=deltaVatempshort_(2);
            deltaVaz(j,i)=deltaVatempshort_(3);
            deltaVbx(j,i)=deltaVbtempshort_(1);
            deltaVby(j,i)=deltaVbtempshort_(2);
            deltaVbz(j,i)=deltaVbtempshort_(3);
        else
            deltaVax(j,i)=deltaVatemplong_(1);
            deltaVay(j,i)=deltaVatemplong_(2);
            deltaVaz(j,i)=deltaVatemplong_(3);
            deltaVbx(j,i)=deltaVbtemplong_(1);
            deltaVby(j,i)=deltaVbtemplong_(2);
            deltaVbz(j,i)=deltaVbtemplong_(3);
        end
    end
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
mesh(tf,td,deltaVtotal,'EdgeColor','black','FaceColor','None')
view(10,20)
set(gca,'GridLineStyle','-')
xlabel('Transfer Time (mins)')
ylabel('Delay Time (mins)')
zlabel('Total Delta V (km/s)')

figure(2)
[x,y,z]=sphere; %Creates unit sphere, to represent Earth
r=6378.137; %Earth radius, km
hsurface=surf(r*x,r*y,r*z); %fills in shpere with grid
set(hsurface,'FaceColor',[0 0 0], 'FaceAlpha', 0.5) %Black, transperency
hold on %Allows other data to be plotted
plot3(rtgtx,rtgty,rtgtz,rintx,rinty,rintz,rtransx,rtransy,rtransz)
axis equal %maintains constant aspect ration (avoids distortions)
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
hold off




