clear all

%% Initial conditions

XFuel=0.016529;
XO2=0.206612;
XN2=0.776859;
Tini=1000;               % Temperature in K
Pini=40;               % Pressure in atm


XEND=1000.0;            % Length of the domain in cm

%% Calculation of Di,ThermiMax,uCJ,tau1,tau2,T1,T2

run ZNDCalc

%% Calculation of reduced activation energy Theta

Theta=(1/Tvn)*((log(tau2)-log(tau1))/((1/T2)-(1/T1)));

%% Calculation of chi

chi=Theta*Di*(ThermiMax/uCJ);

%% Calculation of cell size

% Correlation parameters

A0=30.465860763763;
a1=89.55438805808153;
a2=-130.792822369483;
a3=42.02450507117405;
b1=-0.02929128383850;
b2=1.026325073064710E-5;
b3=-1.031921244571857E-9;

% Calcul cell size

lambda=1E3*Di*(A0+(a1/chi)+(a2/chi^2)+(a3/chi^3)+b1*chi+b2*chi^2+b3*chi^3);

