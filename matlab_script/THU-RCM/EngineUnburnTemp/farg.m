function [h,u,v,s,Y,cp,cv,dlvlT,dlvlp]=farg(p,T,phi,f,fueltype,airscheme);
%
% [h,u,v,s,Y,cp,dlvlT,dlvlp]=farg(p,T,phi,f,fueltype,airscheme)
%
% Routine to determine the state of mixtures of fuel, air
% and residual combustion products at low temperatures.
% Method closely follows that of:
% 1. Ferguson, C.R., 1986, "Internal Combustion Engines", Wiley, p108;
% who uses the results of:
% 2. Hires, S.D., Ekchian, A., Heywood, J.B., Tabaczynski, R.J., and
% Wall, J.C., 1976, "Performance and NOx Emissions Modeling of a Jet
% Ignition Pre-Chamber Stratified Charge Engine", SAE Trans., Vol 85,
% Paper 760161.
% ********************************************************************
% input:
% p,T,phi - pressure (Pa), temperature (K), and equivalence ratio
% f - residual mass fraction; set f=0 if no combustion products
% are present and f=1 if only combustion products are present
% fueltype - 'gasoline', 'diesel', etc - see fueldata.m for full list
% airscheme - 'GMcB' (Gordon and McBride) or 'Chemkin'
% output:
% h - enthalpy (J/kg), u - internal energy (J/kg),
% v - specific volume (m^3/kg), s - entropy (J/kgK),
% Y - mole fractions of 6 species: CO2, H2O, N2, O2, CO, and H2,
% cp - specific heat (J/kgK),
% dlvlT - partial derivative of log(v) wrt log(T)
% dlvlp - partial derivative of log(v) wrt log(p)
% ********************************************************************
[alpha,beta,gamma,delta,Afuel]=fueldata(fueltype);
switch airscheme
    case 'GMcB'
        A=airdata('GMcB_low');
    case 'Chemkin'
        A=airdata('Chemkin_low');
end
Ru=8314.34; % J/kmolK
table=[-1 1 0 0 1 -1]';
M=[44.01 18.02 28.008 32.000 28.01 2.018]'; % kg/kmol
MinMol=1e-25;
dlvlT=1; dlvlp=-1;
eps=0.210/(alpha+0.25*beta-0.5*gamma);
if phi <= 1.0 % stoichiometric or lean
    nu=[alpha*phi*eps beta*phi*eps/2 0.79+delta*phi*eps/2 ...
        0.21*(1-phi) 0 0]';
    dcdT=0;
else % rich
    z=1000/T;
    K=exp(2.743+z*(-1.761+z*(-1.611+z*0.2803)));
    dKdT=-K*(-1.761+z*(-3.222+z*0.8409))/1000;
    a=1-K;
    b=0.42-phi*eps*(2*alpha-gamma)+K*(0.42*(phi-1)+alpha*phi*eps);
    c=-0.42*alpha*phi*eps*(phi-1)*K;
    nu5=(-b+sqrt(b^2-4*a*c))/2/a;
    dcdT=dKdT*(nu5^2-nu5*(0.42*(phi-1)+alpha*phi*eps)+ ...
        0.42*alpha*phi*eps*(phi-1))/(2*nu5*a+b);
    nu=[alpha*phi*eps-nu5 0.42-phi*eps*(2*alpha-gamma)+nu5 ...
        0.79+delta*phi*eps/2 0 nu5 0.42*(phi-1)-nu5]';
end
% mole fractions and molecular weight of residual
tmoles=sum(nu);
Y=nu/tmoles;
Mres=sum(Y.*M);
% mole fractions and molecular weight of fuel-air
fuel=eps*phi/(1+eps*phi);
o2=0.21/(1+eps*phi);
n2=0.79/(1+eps*phi);
Mfa=fuel*(12.01*alpha+1.008*beta+16*gamma+14.01*delta)+ ...
    32*o2+28.02*n2;
% mole fractions of fuel-air-residual gas
Yres=f/(f+Mres/Mfa*(1-f));
Y=Y*Yres;
Yfuel=fuel*(1-Yres);
Y(3)=Y(3)+n2*(1-Yres);
Y(4)=Y(4)+o2*(1-Yres);
% component properties
Tcp0=[1 T T^2 T^3 T^4]';
Th0=[1 T/2 T^2/3 T^3/4 T^4/5 1/T]';
Ts0=[log(T) T T^2/2 T^3/3 T^4/4 1]';
cp0=A(1:6,1:5)*Tcp0;
h0=A(1:6,1:6)*Th0;
s0=A(1:6,[1:5 7])*Ts0;
Mfuel=12.01*alpha+1.008*beta+16.000*gamma+14.01*delta;
a0=Afuel(1); b0=Afuel(2); c0=Afuel(3); d0=Afuel(6); e0=Afuel(7);
cpfuel=Afuel(1:5)*[1 T T^2 T^3 1/T^2]';
hfuel=Afuel(1:6)*[1 T/2 T^2/3 T^3/4 -1/T^2 1/T]';
s0fuel=Afuel([1:5 7])*[log(T) T T^2/2 T^3/3 -1/T^2/2 1]';
% set min value of composition so log calculations work
if Yfuel<MinMol
    Yfuel=MinMol;
end
i=find(Y<MinMol);
Y(i)=ones(length(i),1)*MinMol;
% properties of mixture
h=hfuel*Yfuel+sum(h0.*Y);
s=(s0fuel-log(Yfuel))*Yfuel+sum((s0-log(Y)).*Y);
cp=cpfuel*Yfuel+sum(cp0.*Y)+sum(h0.*table*T*dcdT*Yres/tmoles);
MW=Mfuel*Yfuel+sum(Y.*M);
R=Ru/MW;
h=R*T*h;
u=h-R*T;
v=R*T/p;
s=R*(-log(p/101.325e3)+s);
cp=R*cp;
cv=cp-R;