function [h,u,v,s,Y,cp,cv,dlvlT,dlvlp]=ecp(p,T,phi,fueltype,airscheme,Yguess)
%
% [h,u,v,s,Y,cp,dlvlT,dlvlp]=ecp(p,T,phi,fueltype,airscheme,Yguess)
%
% Routine to determine the equilibrium state of combustion products.
% Method closely follows that of:
% 1. Ferguson, C.R., 1986, "Internal Combustion Engines", Wiley, p122;
% which uses the method described by:
% 2. Olikara, C., and Borman, G.L., 1975, "A Computer Program for
% Calculating Properties of Equilibrium Combustion Products with
% Some Applications to I.C. Engines", SAE Paper 750468.
% ********************************************************************
% input:
% p,T,phi - pressure (Pa), temperature (K), and equivalence ratio
% fueltype - 'gasoline', 'diesel', etc - see fueldata.m for full list
% airscheme - 'GMcB' (Gordon and McBride) or 'Chemkin'
% Yguess - (optional) initial estimate for mole fractions of the
% species CO2 H2O N2 O2 CO H2 H O OH and NO
% output:
% h - enthalpy (J/kg), u - internal energy (J/kg),
% v - specific volume (m^3/kg), s - entropy (J/kgK),
% Y - mole fractions of 10 species, cp - specific heat (J/kgK),
% dlvlT - partial derivative of log(v) wrt log(T)
% dlvlp - partial derivative of log(v) wrt log(p)
% ********************************************************************
[alpha,beta,gamma,delta,Afuel]=fueldata(fueltype);
switch airscheme
    case 'GMcB'
        A0=airdata('GMcB_hi');
    case 'Chemkin'
        A0=airdata('Chemkin_hi');
end
% Equilibrium constant data from Olikara and Borman via Ferguson
Kp=[ 0.432168E+00 -0.112464E+05 0.267269E+01 -0.745744E-04 0.242484E-08
    0.310805E+00 -0.129540E+05 0.321779E+01 -0.738336E-04 0.344645E-08
    -0.141784E+00 -0.213308E+04 0.853461E+00 0.355015E-04 -0.310227E-08
    0.150879E-01 -0.470959E+04 0.646096E+00 0.272805E-05 -0.154444E-08
    -0.752364E+00 0.124210E+05 -0.260286E+01 0.259556E-03 -0.162687E-07
    -0.415302E-02 0.148627E+05 -0.475746E+01 0.124699E-03 -0.900227E-08];
MinMol=1e-25;
tol=3e-12;
Ru=8314.34; % J/kmol.K
M=[44.01 18.02 28.008 32.000 28.01 2.018 1.009 16 17.009 30.004]'; % kg/kmol
dcdT=zeros(4,1);
dcdp=zeros(4,1);
dfdT=zeros(4,1);
dfdp=zeros(4,1);
dYdT=zeros(10,1);
dYdp=zeros(10,1);
B=zeros(4,1);
% check if solid carbon will form
eps=0.210/(alpha+0.25*beta-0.5*gamma);
if phi>(0.210/eps/(0.5*alpha-0.5*gamma))
    error('phi too high - c(s) and other species will form');
end
if nargin==5 % no Yguess so estimate the composition using farg
    [h,u,v,s,Y,cp,cv,dlvlT,dlvlp]=farg(p,T,phi,1,fueltype,airscheme);
    Y(7:10)=ones(4,1)*MinMol; % since farg only returns first 6 species
else
    Y=Yguess;
end
% evaluate constants
patm=p/101.325e3; % convert Pa to atmospheres
TKp=[log(T/1000) 1/T 1 T T^2]';
K=10.^(Kp*TKp);
c=K.*[1/sqrt(patm) 1/sqrt(patm) 1 1 sqrt(patm) sqrt(patm)]';
d=[beta/alpha (gamma+0.42/eps/phi)/alpha (delta+1.58/eps/phi)/alpha]';
if abs(phi-1)<tol
    phi=phi*(1+tol*sign(phi-1));
end
i=find(Y<MinMol);
Y(i)=ones(length(i),1)*MinMol;
DY3to6=2*tol*ones(4,1);
MaxIter=500;
MaxVal=max(abs(DY3to6));
Iter=0;
DoneSome=0;
while (Iter<MaxIter)&&((MaxVal>tol)||(DoneSome<1))
    Iter=Iter+1;
    if Iter>2,
        DoneSome=1;
    end
    D76=0.5*c(1)/sqrt(Y(6));
    D84=0.5*c(2)/sqrt(Y(4));
    D94=0.5*c(3)*sqrt(Y(6)/Y(4));
    D96=0.5*c(3)*sqrt(Y(4)/Y(6));
    D103=0.5*c(4)*sqrt(Y(4)/Y(3));
    D104=0.5*c(4)*sqrt(Y(3)/Y(4));
    D24=0.5*c(5)*Y(6)/sqrt(Y(4));
    D26=c(5)*sqrt(Y(4));
    D14=0.5*c(6)*Y(5)/sqrt(Y(4));
    D15=c(6)*sqrt(Y(4));
    A(1,1)=1+D103;
    A(1,2)=D14+D24+1+D84+D104+D94;
    A(1,3)=D15+1;
    A(1,4)=D26+1+D76+D96;
    A(2,1)=0;
    A(2,2)=2*D24+D94-d(1)*D14;
    A(2,3)=-d(1)*D15-d(1);
    A(2,4)=2*D26+2+D76+D96;
    A(3,1)=D103;
    A(3,2)=2*D14+D24+2+D84+D94+D104-d(2)*D14;
    A(3,3)=2*D15+1-d(2)*D15-d(2);
    A(3,4)=D26+D96;
    A(4,1)=2+D103;
    A(4,2)=D104-d(3)*D14;
    A(4,3)=-d(3)*D15-d(3);
    A(4,4)=0;
    B(1)=-(sum(Y)-1);
    B(2)=-(2*Y(2)+2*Y(6)+Y(7)+Y(9)-d(1)*Y(1)-d(1)*Y(5));
    B(3)=-(2*Y(1)+Y(2)+2*Y(4)+Y(5)+Y(8)+Y(9)+Y(10)-d(2)*Y(1)-d(2)*Y(5));
    B(4)=-(2*Y(3)+Y(10)-d(3)*Y(1)-d(3)*Y(5));
%     invA=inv(A);
    DY3to6=A\B;
    MaxVal=max(abs(DY3to6));
    Y(3:6)=Y(3:6)+DY3to6/10;
    i=find(Y<MinMol);
    Y(i)=ones(length(i),1)*MinMol;
    Y(7)=c(1)*sqrt(Y(6));
    Y(8)=c(2)*sqrt(Y(4));
    Y(9)=c(3)*sqrt(Y(4)*Y(6));
    Y(10)=c(4)*sqrt(Y(4)*Y(3));
    Y(2)=c(5)*sqrt(Y(4))*Y(6);
    Y(1)=c(6)*sqrt(Y(4))*Y(5);
end
if Iter>=MaxIter
    warning('convergence failure in composition loop');
end
TdKdT=[1/T -1/T^2 1 2*T]';
dKdT=2.302585*K.*(Kp(:,[1 2 4 5])*TdKdT);
dcdT(1)=dKdT(1)/sqrt(patm);
dcdT(2)=dKdT(2)/sqrt(patm);
dcdT(3)=dKdT(3);
dcdT(4)=dKdT(4);
dcdT(5)=dKdT(5)*sqrt(patm);
dcdT(6)=dKdT(6)*sqrt(patm);
dcdp(1)=-0.5*c(1)/p;
dcdp(2)=-0.5*c(2)/p;
dcdp(5)=0.5*c(5)/p;
dcdp(6)=0.5*c(6)/p;
x1=Y(1)/c(6);
x2=Y(2)/c(5);
x7=Y(7)/c(1);
x8=Y(8)/c(2);
x9=Y(9)/c(3);
x10=Y(10)/c(4);
dfdT(1)=dcdT(6)*x1+dcdT(5)*x2+dcdT(1)*x7+dcdT(2)*x8+ ...
    dcdT(3)*x9+dcdT(4)*x10;
dfdT(2)=2*dcdT(5)*x2+dcdT(1)*x7+dcdT(3)*x9-d(1)*dcdT(6)*x1;
dfdT(3)=2*dcdT(6)*x1+dcdT(5)*x2+dcdT(2)*x8+dcdT(3)*x9+ ...
    dcdT(4)*x10-d(2)*dcdT(6)*x1;
dfdT(4)=dcdT(4)*x10-d(3)*dcdT(6)*x1;
dfdp(1)=dcdp(6)*x1+dcdp(5)*x2+dcdp(1)*x7+dcdp(2)*x8;
dfdp(2)=2*dcdp(5)*x2+dcdp(1)*x7-d(1)*dcdp(6)*x1;
dfdp(3)=2*dcdp(6)*x1+dcdp(5)*x2+dcdp(2)*x8-d(2)*dcdp(6)*x1;
dfdp(4)=-d(3)*dcdp(6)*x1;
B=-dfdT;
dYdT(3:6)=A\B;
dYdT(1)=sqrt(Y(4))*Y(5)*dcdT(6)+D14*dYdT(4)+D15*dYdT(5);
dYdT(2)=sqrt(Y(4))*Y(6)*dcdT(5)+D24*dYdT(4)+D26*dYdT(6);
dYdT(7)=sqrt(Y(6))*dcdT(1)+D76*dYdT(6);
dYdT(8)=sqrt(Y(4))*dcdT(2)+D84*dYdT(4);
dYdT(9)=sqrt(Y(4)*Y(6))*dcdT(3)+D94*dYdT(4)+D96*dYdT(6);
dYdT(10)=sqrt(Y(4)*Y(3))*dcdT(4)+D104*dYdT(4)+D103*dYdT(3);
B=-dfdp;
dYdp(3:6)=A\B;
dYdp(1)=sqrt(Y(4))*Y(5)*dcdp(6)+D14*dYdp(4)+D15*dYdp(5);
dYdp(2)=sqrt(Y(4))*Y(6)*dcdp(5)+D24*dYdp(4)+D26*dYdp(6);
dYdp(7)=sqrt(Y(6))*dcdp(1)+D76*dYdp(6);
dYdp(8)=sqrt(Y(4))*dcdp(2)+D84*dYdp(4);
dYdp(9)=D94*dYdp(4)+D96*dYdp(6);
dYdp(10)=D104*dYdp(4)+D103*dYdp(3);
% calculate thermodynamic properties
Tcp0=[1 T T^2 T^3 T^4]';
Th0=[1 T/2 T^2/3 T^3/4 T^4/5 1/T]';
Ts0=[log(T) T T^2/2 T^3/3 T^4/4 1]';
cp0=A0(:,1:5)*Tcp0;
h0=A0(:,1:6)*Th0;
s0=A0(:,[1:5 7])*Ts0;
% Y(1) and Y(2) reevaluated
Y(1)=(2*Y(3)+Y(10))/d(3)-Y(5);
Y(2)=(d(1)/d(3)*(2*Y(3)+Y(10))-2*Y(6)-Y(7)-Y(9))/2;
i=find(Y<MinMol);
Y(i)=ones(length(i),1)*MinMol;
% properties of mixture
h=sum(h0.*Y);
s=sum((s0-log(Y)).*Y);
cp=sum(Y.*cp0+h0.*dYdT*T);
MW=sum(Y.*M);
MT=sum(dYdT.*M);
Mp=sum(dYdp.*M);
R=Ru/MW;
v=R*T/p;
cp=R*(cp-h*T*MT/MW);
cv=cp-R;
dlvlT=1+max(-T*MT/MW,0);
dlvlp=-1-max(p*Mp/MW,0);
h=R*T*h;
s=R*(-log(patm)+s);
u=h-R*T;