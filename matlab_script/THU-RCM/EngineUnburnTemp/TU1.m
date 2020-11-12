function [ca_soc,Tu_ko,mu_ko]=TU1(p,v,intq_cal,ca,ca_cal,ca_ko)
%TU1和TU2的区别在与，TU1不计算最高压力点已燃区声速
p=p*1E5;
v=v*1E-9;
dp=centerdiff(ca,p,0.2);
dv=centerdiff(ca,v,0.2);
enginedata;
ca_ivc=-137;
ca_soc=findcax(ca_cal,intq_cal,1);
compspan=ca_ivc:0.2:ca_soc;
compspan=compspan';
p1=interp1(ca,p,compspan,'linear');
v1=interp1(ca,v,compspan,'linear');
T1=p1.*v1/(p1(1)*v1(1))*T_ivc;
T_soc=T1(end);
p_soc=p1(end);
V_soc=v1(end);

Tbb=Tadiabatic(p_soc,T_soc,phi,f,fueltype,airscheme);%以绝热火焰温度作为点火时已燃区温度
mbb=0.001*mass0;%假设着火时已燃区质量为总质量的0.1%
[~,~,vbb,~,~,~,~,~,~]=ecp(p_soc,Tbb,phi,fueltype,airscheme);%求已燃区的体积
Vbb=mbb*vbb;

thetaspanODE=[ca_soc ca_ko];
options = odeset('RelTol',1e-5,'AbsTol',[1e-7 1e-9 1e-8 1E-10 1e-3 1e-3]);
[t,Y]=ode45(@(theta,y) myunbruntemp(theta,y,ca,p,dp,dv), thetaspanODE, [mass0 mbb V_soc Vbb T_soc Tbb],options);
Tu_ko=Y(end,5);
mu_ko=Y(end,1)/Y(1,1);
% plot(t,Y(:,5),'g')
% hold on
% plot(t,Y(:,6),'r')
% hold off
% figure, plot(t,Y(:,3),'g')
% hold on
% plot(t,Y(:,4),'r')
% hold off
% figure, plot(t,Y(:,1),'g')
% hold on
% plot(t,Y(:,2),'r')
% hold off
