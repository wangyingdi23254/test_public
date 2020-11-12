%% ZND calculation for cell size estimation

%% Make chem.bin

system INTERP.exe

%% Find Vcj

DET=200000;
TLIM=2000;
DELX=0.1;

nom='zndkin.inp';
fid = fopen(nom,'w');
fprintf(fid,'%s\n','DCJ');
fprintf(fid,'%s','DET', ' ');
fprintf(fid,'%i\n',DET);
fprintf(fid,'%s','PRES', ' ');
fprintf(fid,'%i\n',Pini);
fprintf(fid,'%s','TEMP', ' ');
fprintf(fid,'%i\n',Tini);
fprintf(fid,'%s','REAC IC8H18', ' ');
fprintf(fid,'%i\n',XFuel);
fprintf(fid,'%s','REAC O2', ' ');
fprintf(fid,'%i\n',XO2);
fprintf(fid,'%s','REAC N2', ' ');
fprintf(fid,'%i\n',XN2);
fprintf(fid,'%s','TLIM', ' ');
fprintf(fid,'%i\n',TLIM);
fprintf(fid,'%s','XEND', ' ');
fprintf(fid,'%i\n',XEND);
fprintf(fid,'%s','DELX', ' ');
fprintf(fid,'%i\n',DELX);
fprintf(fid,'%s','END');
fclose(fid)

system zndkin.exe

M1=importdata('zndkin_Xs.dat',' ');
Headline=M1.textdata;
Matrix=M1.data;

Vcj=Matrix(1,3);

%% Find other parameters

nom='001zndkin.inp';
fid = fopen(nom,'w');
fprintf(fid,'%s\n','DUSR');
fprintf(fid,'%s','DET', ' ');
fprintf(fid,'%i\n',100*Vcj);
fprintf(fid,'%s','PRES', ' ');
fprintf(fid,'%i\n',Pini);
fprintf(fid,'%s','TEMP', ' ');
fprintf(fid,'%i\n',Tini);
fprintf(fid,'%s','REAC IC8H18', ' ');
fprintf(fid,'%i\n',XFuel);
fprintf(fid,'%s','REAC O2', ' ');
fprintf(fid,'%i\n',XO2);
fprintf(fid,'%s','REAC N2', ' ');
fprintf(fid,'%i\n',XN2);
fprintf(fid,'%s','TLIM', ' ');
fprintf(fid,'%i\n',TLIM);
fprintf(fid,'%s','XEND', ' ');
fprintf(fid,'%i\n',XEND);
fprintf(fid,'%s','DELX', ' ');
fprintf(fid,'%i\n',DELX);
fprintf(fid,'%s','END');
fclose(fid)

nom='002zndkin.inp';
fid = fopen(nom,'w');
fprintf(fid,'%s\n','DUSR');
fprintf(fid,'%s','DET', ' ');
fprintf(fid,'%i\n',100*0.99*Vcj);
fprintf(fid,'%s','PRES', ' ');
fprintf(fid,'%i\n',Pini);
fprintf(fid,'%s','TEMP', ' ');
fprintf(fid,'%i\n',Tini);
fprintf(fid,'%s','REAC IC8H18', ' ');
fprintf(fid,'%i\n',XFuel);
fprintf(fid,'%s','REAC O2', ' ');
fprintf(fid,'%i\n',XO2);
fprintf(fid,'%s','REAC N2', ' ');
fprintf(fid,'%i\n',XN2);
fprintf(fid,'%s','TLIM', ' ');
fprintf(fid,'%i\n',TLIM);
fprintf(fid,'%s','XEND', ' ');
fprintf(fid,'%i\n',XEND);
fprintf(fid,'%s','DELX', ' ');
fprintf(fid,'%i\n',DELX);
fprintf(fid,'%s','END');
fclose(fid)

nom='003zndkin.inp';
fid = fopen(nom,'w');
fprintf(fid,'%s\n','DUSR');
fprintf(fid,'%s','DET', ' ');
fprintf(fid,'%i\n',100*1.01*Vcj);
fprintf(fid,'%s','PRES', ' ');
fprintf(fid,'%i\n',Pini);
fprintf(fid,'%s','TEMP', ' ');
fprintf(fid,'%i\n',Tini);
fprintf(fid,'%s','REAC IC8H18', ' ');
fprintf(fid,'%i\n',XFuel);
fprintf(fid,'%s','REAC O2', ' ');
fprintf(fid,'%i\n',XO2);
fprintf(fid,'%s','REAC N2', ' ');
fprintf(fid,'%i\n',XN2);
fprintf(fid,'%s','TLIM', ' ');
fprintf(fid,'%i\n',TLIM);
fprintf(fid,'%s','XEND', ' ');
fprintf(fid,'%i\n',XEND);
fprintf(fid,'%s','DELX', ' ');
fprintf(fid,'%i\n',DELX);
fprintf(fid,'%s','END');
fclose(fid)

system ZNDCalc.bat

M3=importdata('001zndkin_Xs.dat',' ');
Headline3=M3.textdata;
Matrix3=M3.data;

uCJ=Matrix3(2,3);
Thermi=Matrix3(:,54);
ThermiMax=max(Matrix3(:,54));
x=Matrix3(:,1);
Tvn=Matrix3(3,5);
posiMaxXThermi=find(Thermi>=ThermiMax);
Di=(x(posiMaxXThermi(1),1));

figure(1)
semilogx(x,Thermi,'r','LineWidth',2)
xlabel('Distance (m)','FontSize',16,'FontWeight','bold')
ylabel('Thermicity (1/s)','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
%xmin=1E-2;
%xmax=1E1;
%xlim([xmin xmax]);

M4=importdata('002zndkin_Xs.dat',' ');
Headline4=M4.textdata;
Matrix4=M4.data;

Thermi_1=Matrix4(:,54);
ThermiMax_1=max(Matrix4(:,54));
x_1=Matrix4(:,1);
dx_1=diff(x_1);
Vel_1=Matrix4(:,3);
T1=Matrix4(3,5);
posiMaxXThermi_1=find(Thermi_1>=ThermiMax_1);

for i=1:posiMaxXThermi_1
    t1(i)=dx_1(i)/Vel_1(i);
end

tau1=1E6*sum(t1(:));

M5=importdata('003zndkin_Xs.dat',' ');
Headline5=M5.textdata;
Matrix5=M5.data;

Thermi_2=Matrix5(:,54);
ThermiMax_2=max(Matrix5(:,54));
x_2=Matrix5(:,1);
dx_2=diff(x_2);
Vel_2=Matrix5(:,3);
T2=Matrix5(3,5);
posiMaxXThermi_2=find(Thermi_2>=ThermiMax_2);

for i=1:posiMaxXThermi_2
    t2(i)=dx_2(i)/Vel_2(i);
end

tau2=1E6*sum(t2(:));