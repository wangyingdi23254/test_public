function [ca_soc,Tu_ko,mu_ko,c_p_max,kappa_ko]=TU2(p,v,intq_cal,ca,ca_cal,ca_ko,af)
%TU1��TU2���������룬TU1���������ѹ������ȼ������
%ca_ko������ʱ��
%ca_cal�������ʼ��㷶Χ
%af����ȼ��
%intq_cal�����㷶Χ�ڵ��ۼƷ�����
warning off all
p=p*1E5;%��ѹ����barת��Ϊpa
v=v*1E-9;%���ݻ���mm3ת��Ϊm3
dp=centerdiff(ca,p,0.2);%��ѹ����
dv=centerdiff(ca,v,0.2);%���ݻ��仯��

enginedata;

p_ko=interp1(ca,p,ca_ko,'linear');%���뷢������������
v_ko=interp1(ca,v,ca_ko,'linear');
ca_ivc=-137;%�����Źر�ʱ��
ca_soc=findcax(ca_cal,intq_cal,1);%��1%��������Ϊȼ��ʼ��
compspan=ca_ivc:0.2:ca_soc;%ѹ���г̵ķ�Χ
compspan=compspan';
p1=interp1(ca,p,compspan,'linear');%ѹ���г̷�Χ�ڵ�ѹ��ֵ
v1=interp1(ca,v,compspan,'linear');%ѹ���г̷�Χ�ڵ��ݻ�ֵ
T1=p1.*v1/(p1(1)*v1(1))*T_ivc;%ѹ���г̷�Χ�ڵ��¶�ֵ
T_soc=T1(end);%ȼ��ʼ����¶�
p_soc=p1(end);%ȼ��ʼ���ѹ��
V_soc=v1(end);%ȼ��ʼ����ݻ�

Tbb=Tadiabatic(p_soc,T_soc,phi,f,fueltype,airscheme);%�Ծ��Ȼ����¶���Ϊ���ʱ��ȼ���¶�
mbb=0.001*mass0;%�����Ż�ʱ��ȼ������Ϊ��������0.1%
[~,~,vbb,~,~,~,~,~,~]=ecp(p_soc,Tbb,phi,fueltype,airscheme);%����ȼ�������
Vbb=mbb*vbb;%����ȼ���ݻ�

thetaspanODE=[ca_soc ca_ko];%����δȼ���¶ȷ�Χ��ȼ��ʼ�㿪ʼ������ʱ�̽���
options = odeset('RelTol',1e-5,'AbsTol',[1e-7 1e-9 1e-8 1E-10 1e-3 1e-3]);%���㾫������

[t,Y]=ode45(@(theta,y) myunbruntemp(af,theta,y,ca,p,dp,dv), thetaspanODE, [mass0 mbb V_soc Vbb T_soc Tbb],options);%���

Tu_ko=interp1(t,Y(:,5),ca_ko,'linear');%����ʱδȼ���¶�
Tb_ko=interp1(t,Y(:,6),ca_ko,'linear');%����ʱ��ȼ���¶�
% mu_ko=interp1(t,Y(:,1),ca_ko,'linear')/Y(1,1);
mu_ko=interp1(t,Y(:,1),ca_ko,'linear');%����ʱδȼ����
Tb_pmax=Y(end,6);%��ȼ�¶����ֵ
[~,~,~,~,~,cp_ko,cv_ko,~,~]=ecp(p_ko,Tu_ko,phi,fueltype,airscheme);%����ʱδȼ���Ķ�ѹ�Ͷ��ݱ���
kappa_ko=cp_ko/cv_ko;%δȼ���ľ���ָ��
[~,~,~,~,~,cpb,cvb,~,~]=ecp(p_ko,Tb_ko,phi,fueltype,airscheme);%����ʱ��ȼ���Ķ�ѹ�Ͷ��ݱ���
c_b=(cpb/cvb*(cpb-cvb)*Tb_ko)^0.5;%����ʱ��ȼ�������٣������������۱���һ��Ƶ��

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