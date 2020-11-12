function [alpha,beta,gamma,delta,Afuel]=fueldata(fuel)
%
% [alpha,beta,gamma,delta,Afuel]=fueldata(fuel)
%
% Routine to specify the thermodynamic properties of a fuel.
% Data taken from:
% 1. Ferguson, C.R., 1986, "Internal Combustion Engines", Wiley;
% 2. Heywood, J.B., 1988, "Internal Combustion Engine Fundamentals",
% McGraw-Hill; and
% 3. Raine, R. R., 2000, "ISIS_319 User Manual", Oxford Engine Group.
% ********************************************************************
% input:
% fuel switch
% from Ferguson: 'gasoline', 'diesel', 'methane', 'methanol',
% 'nitromethane', 'benzene';
% from Heywood: 'methane_h', 'propane', 'hexane', 'isooctane_h',
% 'methanol_h', 'ethanol', 'gasoline_h1', gasoline_h2', 'diesel_h';
% from Raine: 'toluene', 'isooctane'.
% output:
% alpha, beta, gamma, delta - number of C, H, 0, and N atoms
% Afuel - vector of polynomial coefficients for cp/R, h/RT, and s/R
% of the form h/RT=a1+a2*T/2+a3*T^2/3+a4*T^3/4-a5/T^2+a6/T (for
% example) where T is expressed in K.
% ********************************************************************
% Set values for conversion of Heywood data to nondimensional format
% with T expressed in K
SVal=4.184e3/8.31434;
SVec=SVal*[1e-3 1e-6 1e-9 1e-12 1e3 1 1];
switch fuel
    case 'gasoline' % Ferguson
        alpha=7; beta=17; gamma=0; delta=0;
        Afuel=[4.0652 6.0977E-02 -1.8801E-05 0 0 -3.5880E+04 15.45];
    case 'diesel' % Ferguson
        alpha=14.4; beta=24.9; gamma=0; delta=0;
        Afuel=[7.9710 1.1954E-01 -3.6858E-05 0 0 -1.9385E+04 -1.7879];
    case 'methane' % Ferguson
        alpha=1; beta=4; gamma=0; delta=0;
        Afuel=[1.971324 7.871586E-03 -1.048592E-06 0 0 -9.930422E+03 8.873728];
    case 'methanol' % Ferguson
        alpha=1; beta=4; gamma=1; delta=0;
        Afuel=[1.779819 1.262503E-02 -3.624890E-06 0 0 -2.525420E+04 1.50884E+01];
    case 'nitromethane' % Ferguson
        alpha=1; beta=3; gamma=2; delta=1;
        Afuel=[1.412633 2.087101E-02 -8.142134E-06 0 0 -1.026351E+04 1.917126E+01];
    case 'benzene' % Ferguson
        alpha=6; beta=6; gamma=0; delta=0;
        Afuel=[-2.545087 4.79554E-02 -2.030765E-05 0 0 8.782234E+03 3.348825E+01];
    case 'toluene' % Raine
        alpha=7; beta=8; gamma=0; delta=0;
        Afuel=[-2.09053 5.654331e-2 -2.350992e-5 0 0 4331.441411 34.55418257];
    case 'isooctane' % Raine
        alpha=8; beta=18; gamma=0; delta=0;
        Afuel=[6.678E-1 8.398E-2 -3.334E-5 0 0 -3.058E+4 2.351E+1];
    case 'methane_h' % Heywood
        alpha=1; beta=4; gamma=0; delta=0;
        Afuel=[-0.29149 26.327 -10.610 1.5656 0.16573 -18.331 19.9887/SVal].*SVec;
    case 'propane' % Heywood
        alpha=3; beta=8; gamma=0; delta=0;
        Afuel=[-1.4867 74.339 -39.065 8.0543 0.01219 -27.313 26.4796/SVal].*SVec;
    case 'hexane' % Heywood
        alpha=6; beta=14; gamma=0; delta=0;
        Afuel=[-20.777 210.48 -164.125 52.832 0.56635 -39.836 79.5542/SVal].*SVec;
    case 'isooctane_h' % Heywood
        alpha=8; beta=18; gamma=0; delta=0;
        Afuel=[-0.55313 181.62 -97.787 20.402 -0.03095 -60.751 27.2162/SVal].*SVec;
    case 'methanol_h' % Heywood
        alpha=1; beta=4; gamma=1; delta=0;
        Afuel=[-2.7059 44.168 -27.501 7.2193 0.20299 -48.288 31.1406/SVal].*SVec;
    case 'ethanol' % Heywood
        alpha=2; beta=6; gamma=1; delta=0;
        Afuel=[6.990 39.741 -11.926 0 0 -60.214 8.01623/SVal].*SVec;
    case 'gasoline_h1' % Heywood
        alpha=8.26; beta=15.5; gamma=0; delta=0;
        Afuel=[-24.078 256.63 -201.68 64.750 0.5808 -27.562 NaN].*SVec;
    case 'gasoline_h2' % Heywood
        alpha=7.76; beta=13.1; gamma=0; delta=0;
        Afuel=[-22.501 227.99 -177.26 56.048 0.4845 -17.578 NaN].*SVec;
    case 'diesel_h' % Heywood
        alpha=10.8; beta=18.7; gamma=0; delta=0;
        Afuel=[-9.1063 246.97 -143.74 32.329 0.0518 -50.128 NaN].*SVec;
end