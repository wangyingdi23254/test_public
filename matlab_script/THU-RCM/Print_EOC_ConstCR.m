function Print_EOC_ConstCR()
clear; clc;
%%
T_0=273.15;
T_initial = T_0+25 ;% T_initial=T_initial;  %从这里输入初始温度
phi = 1.0;  %燃空当量比
Innert = 3.76;% Innert / oxygen
atom = [1 4 1];% C H O
Oxy2fuel = atom(1)*1 +atom(2)/4 - atom(3)/2;% Oxygen / fuel

filename = [ 'methanol_VCR' '.txt']; %输出文件名按时间命名，记得更改
fid = fopen(filename,'w');

for CR = 10:0.5:16
    fprintf( fid,'CR = %.2f\n',CR );
    for Rat = 0:0.5:0 %Rat = Mole(Inert/N2)
        fprintf( fid,'Ratio = %.2f\n',Rat );
        for P_initial = 1:0.2:2
            %for T_initial = 273.15 + 60:20:273.15+80.5  %可能会用到加热，初始温度变动范围
            [T_eff, P_eff ] = Calc_T_P(CR, T_initial, P_initial, phi, Rat,Innert,Oxy2fuel);
            fprintf( fid,'P_initial = %.2f\t T_initial = %.2f\t P_eff = %.2f\t T_eff = %.2f\t\n',...
            P_initial,T_initial,P_eff, T_eff );
            %end
        end
    end
end
fclose(fid);

function [T_eff, P_eff] = Calc_T_P(CR, T_initial, P_initial, phi, Rat,Innert,Oxy2fuel)
% 根据不同组分、初始温度、压缩比，由绝热假设来计算压缩终止温度。
% 用于实验前初步计算终止温度。
V_initial = 1;
dT=0.1; 
T_temp = T_initial;
Molar_sum = phi+Oxy2fuel*(1+Innert);
Molar = struct(...
    'O2', Oxy2fuel./Molar_sum,'Ar',(Oxy2fuel*Innert*Rat/(Rat+1))./Molar_sum,...%配气组分
    'N2',(Oxy2fuel*Innert/(Rat+1))./Molar_sum ,'H2',0,'CO2',0, 'CO',0,...
    'H2O',0,'nHeptane',0,'CH4',0,'Toluene',0,'Isooctane',0,...
    'Methylcyclohexane',0, 'Isobutanol',0,'CH3OH',phi./Molar_sum ,'nButane',0 ,'DMM3',0);
calc=0;
while calc<log(CR)
   calc = calc + 1./ (gamma_all(T_temp,Molar.O2, Molar.N2, Molar.Ar, Molar.H2, Molar.CO2, ...
       Molar.CO,Molar.CH4, Molar.H2O, Molar.nHeptane, Molar.Toluene, Molar.Isooctane, ...
       Molar.Isobutanol, Molar.CH3OH, Molar.nButane, Molar.DMM3)-1)./T_temp.*dT;
   % CALL function gamma
   T_temp = T_temp +dT;
end
V_eff = V_initial/CR;
T_eff = T_temp;
P_eff = P_initial*V_initial*T_eff/( V_eff*T_initial );
end

function [gamma]=gamma_all(T, Molar_O2, Molar_N2, Molar_Ar, Molar_H2, Molar_CO2, ...
    Molar_CO,Molar_CH4, Molar_H2O, Molar_nHeptane, Molar_Toluene, Molar_Isooctane, ...
    Molar_Isobutanol, Molar_CH3OH, Molar_nButane, Molar_DMM3)
% global Molar_H2 Molar_O2 Molar_Ar Molar_N2 Molar_CO2 Molar_H2O Molar_CO Molar_isooctane
% 1  2   3   4   5   6    7   8    9         10       11         12
% T, O2, N2, Ar, H2, CO2, CO, H2O, nHeptane, Toluene, Isooctane, nButane
% T=x;  %Mixture temperature
%Gas molar fraction

%Universal gas constant
R=8.31451;        %J/(mol-K)
%Cp/R of isooctane  MW=114.23092
if(T>1000)    
    %Cp/R of N2  MW=28.01348
    N2=2.95257626+1.39690057E-3.*T-4.92631691E-7.*T.^2+7.86010367E-11.*T.^3-4.60755321E-15.*T.^4;
    %Cp/R of O2   MW=31.9988
    O2=3.66096083+6.56365523E-4.*T-1.41149485E-7.*T.^2+2.05797658E-11.*T.^3-1.29913248E-15.*T.^4;
    CO2=4.63659493+2.74131991E-3.*T-9.95828531E-7.*T.^2+1.60373011E-10.*T.^3-9.16103468E-15.*T.^4;
    H2O=2.67703787+2.97318329E-3.*T-7.73769690E-7.*T.^2+9.44336689E-11.*T.^3-4.26900959E-15.*T.^4;
    CO=3.04848583+1.35172818E-3.*T-4.85794075E-7.*T.^2+7.88536486E-11.*T.^3-4.69807489E-15.*T.^4;
    H2=3.33727920E+00-4.94024731E-05.*T+4.99456778E-07.*T.^2-1.79566394E-10.*T.^3+2.00255376E-14.*T.^4;
    CH4=7.48514950E-02+ 1.33909467E-02.*T-5.73285809E-06.*T^2 + 1.22292535E-09.*T^3 -1.01815230E-13.*T^4;
    CH3OH= 3.52726795E+00 + 1.03178783E-02.*T -3.62892944E-06.*T^2 + 5.77448016E-10.*T^3 -3.42182632E-14.*T^4;
else
   
    N2=3.53100528-1.23660987E-4.*T-5.02999437E-7.*T.^2+2.43530612E-09.*T.^3-1.40881235E-12.*T.^4;
    O2=3.78245636-2.99673415E-3.*T+9.84730200E-6.*T.^2-9.68129508E-09.*T.^3+3.24372836E-12.*T.^4;
    CO2=2.35677352+8.98459677E-3.*T-7.12356269E-6.*T.^2+2.45919022E-09.*T.^3-1.43699548E-13.*T.^4;
    H2O=4.19864056-2.03643410E-3.*T+6.52040211E-6.*T.^2-5.48797062E-09.*T.^3+1.77197817E-12.*T.^4;
    CO=3.57953347-6.10353680E-4.*T+1.01681433E-6.*T.^2+9.07005884E-10.*T.^3-9.04424499E-13.*T.^4;
    H2=2.34433112E+00+7.98052075E-03.*T-1.94781510E-05.*T.^2+2.01572094E-08.*T.^3-7.37611761E-12.*T.^4;
    CH4=5.14987613E+00 -1.36709788E-02.*T + 4.91800599E-05.*T^2 -4.84743026E-08.*T^3+ 1.66693956E-11.*T^4;
    CH3OH= 5.65851051E+00 -1.62983419E-02.*T +6.91938156E-05.*T^2 -7.58372926E-08.*T^3 +2.80427550E-11.*T^4;  
end

if (T>1391)
    n_Heptane = 2.22148969E+01+3.47675750E-02.*T-1.18407129E-05.*T.^2+1.83298478E-09.*T.^3-1.06130266E-13.*T.^4; %last number should be 4.
else
    n_Heptane =-1.26836187E+00+8.54355820E-02.*T-5.25346786E-05.*T.^2+1.62945721E-08.*T.^3-2.02394925E-12.*T.^4;
end
if (T>1396)
    Toluene =   1.69989044e+01+2.20052825e-02.*T-7.58020314e-06.*T.^2+1.18267678e-09.*T.^3-6.88594262e-14.*T.^4;
    Isooctane = 2.71373590E+01+3.79004890E-02.*T-1.29437358E-05.*T.^2+2.00760372E-09.*T.^3-1.16400580E-13.*T.^4;
else
    Toluene   =-5.45944164e+00+7.71089443e-02.*T-5.99305765e-05.*T.^2+2.40404364e-08.*T.^3-3.92116250e-12.*T.^4;
    Isooctane =-4.20868893E+00+1.11440581E-01.*T-7.91346582E-05.*T.^2+2.92406242E-08.*T.^3-4.43743191E-12.*T.^4;
end
if (T>1394)
    Isobutanol = 1.47423131e+01 + 2.19843767e-02.*T -7.51584192e-06.*T.^2 +1.16633393e-09.*T.^3 -6.76421638e-14.*T.^4;
else
    Isobutanol =-8.37465362e-01 + 5.76520639e-02.*T - 3.90215462e-05.*T.^2 +1.40131231e-08.*T.^3-2.11159111e-12.*T.^4;
end
if (T>1392)
    n_Butane=1.24940183E+01+2.17726258E-02.*T-7.44272215E-06.*T.^2+1.15487023E-09.*T.^3-6.69712949E-14.*T.^4; 
else
    n_Butane=-4.55756824E-01+4.80323389E-02.*T-2.65497552E-05.*T.^2+6.92544700E-09.*T.^3-6.38317504E-13.*T.^4;
end
if (T>1000)
    DMM3=1.9899747e+01+3.9912973e-02.*T-1.9118063e-05.*T^2+4.4929191e-09.*T^3-4.2155413e-13.*T^4;
else
    DMM3=-1.6150180e+01+2.3867562e-01.*T-4.2010835e-04.*T^2+3.5798712e-07.*T^3-1.1574661e-10.*T^4;
end
%Cp/R of Ar   MW=40
Ar=2.5;

Cp=(H2*Molar_H2+N2*Molar_N2+O2*Molar_O2+Ar*Molar_Ar+CO2*Molar_CO2+...
    H2O*Molar_H2O+CO*Molar_CO+CH4*Molar_CH4+n_Heptane*Molar_nHeptane+...
    Toluene * Molar_Toluene + Isooctane*Molar_Isooctane+...
    Isobutanol*Molar_Isobutanol+n_Butane*Molar_nButane+...
    CH3OH*Molar_CH3OH+DMM3*Molar_DMM3)*R;
Cv=Cp-R;
gamma=Cp/Cv;
end

end