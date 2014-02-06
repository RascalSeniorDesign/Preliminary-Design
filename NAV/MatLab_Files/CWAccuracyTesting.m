clear
clc

vtgt_=[-1.962372; 7.323674; 0.000000]; %Target velocity, km/s
vint_=[-4.864779; 5.816486; .240163]; %Interceptor velocity, km/s
rtgt_=[6697.4756; 1794.5831; 0.000]; %Target Position, km
rint_=[5328.7862; 4436.1273; 101.4720]; %Interceptor Position, km

Xtgt_=[rtgt_ vtgt_];
Xint_=[rint_ vint_];

[x0,y0,z0,x0dot,y0dot,z0dot] = ECI2CW(Xtgt_,Xint_);

tf=linspace(0,400,100);

x=zeros(length(tf),length(tf));
y=zeros(length(tf),length(tf));
z=zeros(length(tf),length(tf));
xdot=zeros(length(tf),length(tf));
ydot=zeros(length(tf),length(tf));
zdot=zeros(length(tf),length(tf));
deltaV=zeros(length(tf),length(tf));

rtgt=sqrt(sum(abs(rtgt_).^2));

for i=1:length(tf)
    [x0dot(i),y0dot(i),z0dot(i)] = CWDocking(x0*1000,y0*1000,z0*1000,rtgt,tf(i));
    deltaV(i)=sqrt(abs(x0dot(i))^2+abs(y0dot(i))^2+abs(z0dot(i))^2);
    for j=1:length(tf)
        [x(i,j),y(i,j),z(i,j),xdot(i,j),ydot(i,j),zdot(i,j)] = CWSolver(x0*1000,y0*1000,z0*1000,x0dot(i)*1000,y0dot(i)*1000,z0dot(i)*1000,rtgt,tf(j));
    end
end

figure (1)
plot(tf,deltaV,tf,x0dot,tf,y0dot,tf,z0dot)
xlabel('Time of Transfer (min)')
ylabel('Delta V Required (min)')

figure(2)
subplot(2,1,1)
plot(tf,x(1,:),tf,y(1,:),tf,z(1,:))
xlabel('Time (min)')
ylabel('Relative Position (m)')
legend('x','y','z')
% 
subplot(2,1,2)
plot(tf,xdot(1,:),tf,ydot(1,:),tf,zdot(1,:))
xlabel('Time (min)')
ylabel('Relative Velocity (m/s)')
legend('dx','dy','dz')

figure (3)
plot3(x(1,:),y(1,:),z(1,:))
grid on
set(gca,'GridLineStyle','-')
xlabel('X Relative Position (m)')
ylabel('Y Relative Position (m)')

figure(4)
subplot(2,1,1)
mesh(tf,tf,x,'EdgeColor','black','FaceColor','None')
view(44,28)
grid on
set(gca,'GridLineStyle','-')
subplot(2,1,2)
mesh(tf,tf,y,'EdgeColor','black','FaceColor','None')
view(44,28)
grid on
set(gca,'GridLineStyle','-')

figure(5)
mesh(tf,tf,z,'EdgeColor','black','FaceColor','None')
view(44,28)
grid on
set(gca,'GridLineStyle','-')
xlabel('Transfer Time (m)')
ylabel('Plot Time (m)')
zlabel('Relative Y Position (m)')