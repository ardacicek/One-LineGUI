function [Ymid,xmid,maxRs,duration]=OneLine(H0,T,alfa0,m,d50,duration,dt,GroinSt,GroinFin)
%% A function to find shoreline changes in a defined duration,calculated in each timestep dt.
%% Initial conditions
%% Consant parameters
% duration = 24*3600*duration; %day to seconds
% dt = 24*3600*dt;
ro_s=2650; %sediment density kg/m3
ro_w=1025; %water density kg/m3
s=ro_s/ro_w; %sediment specific gravity
p=0.4; %sand porosity
d50 = d50/1000;
%% Breaking wave characteristics
[H0,L0,T,~,Hb,Lb, alfab]=MonochromaticBreaking(H0,T,alfa0,m);
Cgb=Lb/T;
K1=0.0001/d50; %https://etd.lib.metu.edu.tr/upload/12608582/index.pdf pg27
beta1=K1/(16*(s-1)*(1-p));
alfabs_init=alfab;

%% Setting up the grid
xmid=repmat(linspace(0,1500,500),round(duration/dt),1); %adjust dx to satisfy Rs
ymid=repmat(zeros(1,size(xmid,2)),size(xmid,1),1); 

Dc=H0*(2.28-10.9*H0/L0); %depth of closure

%% Initializing the simulation
Q_init=(Hb^2*Cgb)*beta1*sin(2*alfabs_init);
Q=zeros(round(duration/dt),size(xmid,2)+1); %edit 2000 manually according to dur/dt
Q(:,1)=Q_init; %Left boundary condition for Q
Q(:,size(xmid,2)+1)=Q_init; %Right boundary condition for Q
dx=xmid(1,2)-xmid(1,1); %step size of x

%% For the beach part before the groin
alfas=zeros(round(duration/dt),size(xmid,2)-1);
alfabs=zeros(round(duration/dt),size(xmid,2)-1);
Ymid=repmat(zeros(1,size(xmid,2)),size(xmid,1),1);
for i=1:round(duration/dt)
    for j=1:GroinSt-2
        alfas(i,j)=(atan((ymid(i,j+1)-ymid(i,j))/dx));
        alfabs(i,j)=alfab-alfas(i,j);
        Q(i,j+1)=(Hb^2*Cgb)*beta1*sin(2*alfabs(i,j));
    end
        Q(i,GroinSt:GroinFin)=0;
    for j=1:GroinSt
        Ymid(i,j)=ymid(i,j)+dt/(Dc*dx)*(Q(i,j)-Q(i,j+1));
        ymid(i+1,j)=Ymid(i,j);
    end
end
%% For the beach part after the groin

for i=1:round(duration/dt)
    for j=GroinFin:size(xmid,2)-1
        alfas(i,j)=(atan((ymid(i,j+1)-ymid(i,j))/dx));
        alfabs(i,j)=alfab-alfas(i,j);
        Q(i,j+1)=(Hb^2*Cgb)*beta1*sin(2*alfabs(i,j));
    end
    for j=GroinFin:size(xmid,2)
        Ymid(i,j)=ymid(i,j)+dt/(Dc*dx)*(Q(i,j)-Q(i,j+1));
        ymid(i+1,j)=Ymid(i,j);
    end
    
end   
%% Stability Check
Rs=Q./(alfab*Dc*dx^2)*dt;
maxRs=max(max(Rs));
% if maxRs<0.5
%     
%     cprintf('*green','Stability Condition is Satisfied.\n');
% else
%     cprintf('*red','!Check dt or dx to Satisfy Stability Condition!\n');
% end
%  duration=duration/3600/24; %seconds to day
%  cprintf('*blue','Duration=%f days\n',duration);
%% Plotting
% plot(xmid(1,:),Ymid(1,:));                              %Initial shape
% hold on 
% % plot(xmid(round(size(xmid,1)/64),:),Ymid(round(size(Ymid,1)/64),:));
% % plot(xmid(round(size(xmid,1)/16),:),Ymid(round(size(Ymid,1)/16),:));  %@t=duration/16 shape
% % plot(xmid(round(size(xmid,1)/8),:),Ymid(round(size(Ymid,1)/8),:));    %@t=duration/8 shape
% % plot(xmid(round(size(xmid,1)/4),:),Ymid(round(size(Ymid,1)/4),:));   %@t=duration/4 shape
% % plot(xmid(round(size(xmid,1)/2),:),Ymid(round(size(Ymid,1)/2),:));   %@t=duration/2 shape
% plot(xmid(size(xmid,1),:),Ymid(size(Ymid,1),:));        %@Final shape
% % legend('Initial Profile',['Duration (days)= ' num2str(round(duration/64))],['Duration (days)= ' num2str(round(duration/16))],['Duration (days)= ' num2str(round(duration/8))],['Duration (days)= ' num2str(round(duration/4))],['Duration (days)= ' num2str(round(duration/2))],['Duration (days)= ' num2str(round(duration))]);
% legend('Initial Profile','Final');

end

