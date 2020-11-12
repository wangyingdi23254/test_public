% enginedata.m

% ***** engine geometry **********************************************
b=0.095; % engine bore (m)
stroke=0.115; % engine stroke (m)
r=stroke/2;
crl=215;%connection rod length
eps=r/crl; % half stroke to rod ratio, s/2l
cr=14; % compression ratio

ca_alpha=linspace(-2*pi,2*pi,7201);%曲轴转角
ca_beta=asin(eps*sin(ca_alpha));                                %连杆与气缸轴线的夹角
x=(r+crl)-(r*cos(ca_alpha)+crl*cos(ca_beta));                           %活塞位移
V=(stroke/(cr-1)+x)*pi*(b/2)^2;                                         %燃烧室容积
vs=stroke*pi*(b/2)^2;
% ***** engine thermofluids parameters *******************************
f=0.05; % residual fraction
fueltype='gasoline';
Hu=44000;
L0=14.7;
airscheme='GMcB';
phi=L0/14.7; % equivalence ratio
RPM=1200;
omega=RPM*pi/30; % engine speed in rad/s
fc=5000;
fs=36000;


tw=420; % engine surface temperature
% ***** initial conditions *******************************************
p_ivc=100e3;%进气门关闭时刻压力
T_ivc=330;%进气门关闭时刻温度
theta1=-pi;
ca_ivc=-137;%进气门关闭时刻;
V_ivc=interp1(ca_alpha,V,ca_ivc*pi/180,'linear');%进气门关闭时刻气缸容积
[~,~,v0,~,~,~,~,~]=farg(p_ivc,T_ivc,phi,f,fueltype,airscheme);%进气门关闭时刻气缸比容积
mass0=V_ivc/v0;%进气门关闭时刻气缸内工质质量