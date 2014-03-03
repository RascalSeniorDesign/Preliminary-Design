clear
clc

drfx=linspace(.001,10,3);
drfmatrix=[drfx;zeros(1,length(drfx));zeros(1,length(drfx))];
field1='deltaVtot';
dvx0=linspace(0,0.1,20);
dvmatrix0=[dvx0;zeros(1,length(dvx0));zeros(1,length(dvx0))];

dr0_=[0;0;0];
rtgt_=[6697.4756; 1794.5831; 0.0];
drf_=[0;1;0];
t=linspace(60*.1,60*300,200);

deltaV=struct(field1,{});

for k=1:length(drfx)
    for i=1:length(dvx0)
        for j=1:length(t)
            [deltaVi,deltaVf]=CWPrussingRendezvousSolver(dr0_,dvmatrix0(:,i),rtgt_,drfmatrix(:,k),t(j));
            deltaVvox(i,j)=sum(sqrt((abs(deltaVi)+abs(deltaVf)).^2));
        if deltaVvox(i,j)>=.15
            deltaVvox(i,j)=.15;
        end
        end
    end
    deltaV(k).deltaVtot=deltaVvox(:,:);
end
figure (2)

for i=1:length(drfx)
    hold on;
    surf(t./60,dvx0*1000,deltaV(i).deltaVtot*1000,'FaceAlpha',0.6)
    legendinfo{i}=['Final Relative Position (m): ' num2str(drfx(i)*1000)];
end
hold off
set(gca,'GridLineStyle','-')
xlabel('Transfer Time (min)')
ylabel('Initial Relative Velocity (m/s)')
zlabel('Total Delta V (m/s)')
legend(legendinfo)
view(13,38)




