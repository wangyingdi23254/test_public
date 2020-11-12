% enginedata.m

% ***** engine geometry **********************************************
b=0.095; % engine bore (m)
stroke=0.115; % engine stroke (m)
r=stroke/2;
crl=215;%connection rod length
eps=r/crl; % half stroke to rod ratio, s/2l
cr=14; % compression ratio

ca_alpha=linspace(-2*pi,2*pi,7201);%����ת��
ca_beta=asin(eps*sin(ca_alpha));                                %�������������ߵļн�
x=(r+crl)-(r*cos(ca_alpha)+crl*cos(ca_beta));                           %����λ��
V=(stroke/(cr-1)+x)*pi*(b/2)^2;                                         %ȼ�����ݻ�
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
p_ivc=100e3;%�����Źر�ʱ��ѹ��
T_ivc=330;%�����Źر�ʱ���¶�
theta1=-pi;
ca_ivc=-137;%�����Źر�ʱ��;
V_ivc=interp1(ca_alpha,V,ca_ivc*pi/180,'linear');%�����Źر�ʱ�������ݻ�
[~,~,v0,~,~,~,~,~]=farg(p_ivc,T_ivc,phi,f,fueltype,airscheme);%�����Źر�ʱ�����ױ��ݻ�
mass0=V_ivc/v0;%�����Źر�ʱ�������ڹ�������