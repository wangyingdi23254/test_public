clear
close all
tic
load p_knock.mat
p=p*1E5;
v=v*1E-9;
dp=centerdiff(ca,p,0.2);
dv=centerdiff(ca,v,0.2);
ca_ivc=-137;
ca_spark=-5;
compspan=ca_ivc:0.2:ca_spark;
compspan=compspan';
T_ivc=330;
p1=interp1(ca,p,compspan,'linear');
v1=interp1(ca,v,compspan,'linear');
T1=p1.*v1/(p1(1)*v1(1))*T_ivc;
T_spark=T1(end);
p_spark=p1(end);
V_spark=v1(end);
enginedata;
Tbb=Tadiabatic(p_spark,T_spark,phi,f,fueltype,airscheme);%以绝热火焰温度作为点火时已燃区温度
mbb=0.001*mass0;%假设着火时已燃区质量为总质量的0.1%
[~,~,vbb,~,~,~,~,~,~]=ecp(p_spark,Tbb,phi,fueltype,airscheme);%求已燃区的体积
Vbb=mbb*vbb;

thetaspanODE=[-5 30];
options = odeset('RelTol',1e-5,'AbsTol',[1e-7 1e-9 1e-8 1E-10 1e-3 1e-3]);
[T,Y]=ode45(@(theta,y) myunbruntemp(theta,y,ca,p,dp,dv), thetaspanODE, [mass0 mbb V_spark Vbb T_spark Tbb],options);
plot(T,Y(:,5),'g')
hold on
plot(T,Y(:,6),'r')
hold off
figure, plot(T,Y(:,3),'g')
hold on
plot(T,Y(:,4),'r')
hold off
figure, plot(T,Y(:,1),'g')
hold on
plot(T,Y(:,2),'r')
hold off
toc