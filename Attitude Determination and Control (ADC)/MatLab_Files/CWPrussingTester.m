clear
clc

dvx0=linspace(0,0.01,20);
dvmatrix0=[dvx0;zeros(1,length(dvx0));zeros(1,length(dvx0))];
dr0_=[0;0;0];
rtgt_=[6697.4756; 1794.5831; 0.0];
t=linspace(0,60*90,200);

for i=1:length(dvx0)
    for j=1:length(t)
        [dr_,dv_]=CWPrussing(dr0_,dvmatrix0(:,i),rtgt_,t(j));
        drx(i,j)=dr_(1);
        dry(i,j)=dr_(2);
        drz(i,j)=dr_(3);
        dvx(i,j)=dv_(1);
        dvy(i,j)=dv_(2);
        dvz(i,j)=dv_(3);
    end
end
figure (1)
mesh(t,dvx0,drx,'EdgeColor','black','FaceColor','None')
hold on
mesh(t,dvx0,dry,'EdgeColor','blue','FaceColor','None')
mesh(t,dvx0,drz,'EdgeColor','red','FaceColor','None')
set(gca, 'GridLineStyle','-')
view(-138,14)
ylabel('Initial Relative X Velocity (km/s)')
xlabel('Simulation Time (s)')
zlabel('Relative Position (km)')

figure (2)
mesh(t,dvx0,dvx,'EdgeColor','black','FaceColor','None')
hold on
mesh(t,dvx0,dvy,'EdgeColor','blue','FaceColor','None')
mesh(t,dvx0,dvz,'EdgeColor','red','FaceColor','None')
set(gca, 'GridLineStyle','-')
view(-138,14)
ylabel('Initial Relative X Velocity (km/s)')
xlabel('Simulation Time (s)')
zlabel('Relative Velocity (km/s)')