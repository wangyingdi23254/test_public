function Tb=Tadiabatic(p,Tu,phi,f,fueltype,airscheme);
%
% Tb=Tadiabatic(p,Tu,phi,f,fueltype,airscheme)
%
% Routine for calculating the adiabatic flame temperature.
% Method involves iteratively selecting flame temperatures until
% the enthalpy of the combustion products (in equilibrium) matches
% the enthalpy of the initial gas mixture.
% farg.m is used to determine the enthalpy of the unburned mixture,
% and ecp.m is used to determine the enthalpy of the burned gas.
% ********************************************************************
% input:
% p - pressure (Pa)
% Tu - temperature of the unburned mixture (K)
% phi - equivalence ratio
% f - residual mass fraction; set f=0 if no combustion products
% are present and f=1 if only combustion products are present
% fueltype - 'gasoline', 'diesel', etc - see fueldata.m for full list
% airscheme - 'GMcB' (Gordon and McBride) or 'Chemkin'
% output:
% Tb - temperature of the burned gas (K) - adiabatic flame temperature
% ********************************************************************
MaxIter=50;
Tol=0.00001; % 0.001% allowable error in temperature calculation
Tb=2000; % initial estimate
DeltaT=2*Tol*Tb; % something big
Iter=0;
[hu,u,v,s,Y,cp,dlvlT,dlvlp]=farg(p,Tu,phi,f,fueltype,airscheme);
while (Iter<MaxIter)&(abs(DeltaT/Tb)>Tol)
Iter=Iter+1;
[hb,u,v,s,Y,cp,dlvlT,dlvlp]=ecp(p,Tb,phi,fueltype,airscheme);
DeltaT=(hu-hb)/cp;
Tb=Tb+DeltaT;
end
if Iter>=MaxIter
warning('convergence failure in adiabatic flame temperature loop');
end