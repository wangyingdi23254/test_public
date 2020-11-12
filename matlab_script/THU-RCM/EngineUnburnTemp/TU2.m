function [ca_soc,Tu_ko,mu_ko,c_p_max,kappa_ko]=TU2(p,v,intq_cal,ca,ca_cal,ca_ko,af)
%TU1和TU2的区别在与，TU1不计算最高压力点已燃区声速
%ca_ko，爆震时刻
%ca_cal，放热率计算范围
%af，空燃比
%intq_cal，计算范围内的累计放热量
warning off all
p=p*1E5;%将压力由bar转换为pa
v=v*1E-9;%将容积由mm3转换为m3
dp=centerdiff(ca,p,0.2);%求压升率
dv=centerdiff(ca,v,0.2);%求容积变化率

enginedata;

p_ko=interp1(ca,p,ca_ko,'linear');%载入发动机基本参数
v_ko=interp1(ca,v,ca_ko,'linear');
ca_ivc=-137;%进气门关闭时刻
ca_soc=findcax(ca_cal,intq_cal,1);%以1%放热率作为燃烧始点
compspan=ca_ivc:0.2:ca_soc;%压缩行程的范围
compspan=compspan';
p1=interp1(ca,p,compspan,'linear');%压缩行程范围内的压力值
v1=interp1(ca,v,compspan,'linear');%压缩行程范围内的容积值
T1=p1.*v1/(p1(1)*v1(1))*T_ivc;%压缩行程范围内的温度值
T_soc=T1(end);%燃烧始点的温度
p_soc=p1(end);%燃烧始点的压力
V_soc=v1(end);%燃烧始点的容积

Tbb=Tadiabatic(p_soc,T_soc,phi,f,fueltype,airscheme);%以绝热火焰温度作为点火时已燃区温度
mbb=0.001*mass0;%假设着火时已燃区质量为总质量的0.1%
[~,~,vbb,~,~,~,~,~,~]=ecp(p_soc,Tbb,phi,fueltype,airscheme);%求已燃区的体积
Vbb=mbb*vbb;%求已燃区容积

thetaspanODE=[ca_soc ca_ko];%计算未燃区温度范围从燃烧始点开始到爆震时刻结束
options = odeset('RelTol',1e-5,'AbsTol',[1e-7 1e-9 1e-8 1E-10 1e-3 1e-3]);%计算精度条件

[t,Y]=ode45(@(theta,y) myunbruntemp(af,theta,y,ca,p,dp,dv), thetaspanODE, [mass0 mbb V_soc Vbb T_soc Tbb],options);%求解

Tu_ko=interp1(t,Y(:,5),ca_ko,'linear');%爆震时未燃区温度
Tb_ko=interp1(t,Y(:,6),ca_ko,'linear');%爆震时已燃区温度
% mu_ko=interp1(t,Y(:,1),ca_ko,'linear')/Y(1,1);
mu_ko=interp1(t,Y(:,1),ca_ko,'linear');%爆震时未燃质量
Tb_pmax=Y(end,6);%已燃温度最大值
[~,~,~,~,~,cp_ko,cv_ko,~,~]=ecp(p_ko,Tu_ko,phi,fueltype,airscheme);%求爆震时未燃区的定压和定容比热
kappa_ko=cp_ko/cv_ko;%未燃区的绝热指数
[~,~,~,~,~,cpb,cvb,~,~]=ecp(p_ko,Tb_ko,phi,fueltype,airscheme);%求爆震时已燃区的定压和定容比热
c_b=(cpb/cvb*(cpb-cvb)*Tb_ko)^0.5;%爆震时已燃区的声速，用来计算理论爆震一阶频率

p_cal=interp1(ca,p,ca_cal,'linear');
v_cal=interp1(ca,v,ca_cal,'linear');
p_lpf=avlfilter(ca_cal,p_cal,fc,fs,'low');
p_lpf_max=max(p_lpf);
v_lpf_max=v_cal(p_lpf==p_lpf_max);
c_p_max=(cpb/cvb*p_lpf_max(1)*v_lpf_max(1)/mass0)^0.5;


% p_max=max(p);
% v_p_max=v(p==p_max);
% c_p_max=(cpb/cvb*p_max(1)*v_p_max(1)/mass0)^0.5;

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