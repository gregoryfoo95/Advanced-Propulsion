%% Advanced Propulsion Coursework %%
close all
clear
clc
%% Professionally Crafted By:
%% Gregory Foo
%% Jonathan Ting
%% Ngiam Yen Kui

%% Definitions of Geometry
%(1) - Intake(Freestream)
%(C1) - Inlet Throat
%(x) - Upstream of Shock
%(y) - Downstream of Shock
%(2) - Beginning of Burner
%(b) - End of Burner (Burn Complete)
%(C2) - Nozzle Throat
%(4) - Engine Exhaust

%% Code Organisation Parameters
font=13;
titl=15;
%% Notations
%Stagnation Temperature - T0
%Static Temperature - T
%Stagnation Pressure - P0
%Static Pressure - P
%Mach Number - M
%Area Ratio - A
%Reference Throat Area - A*
%% Assumptions for Design

%(1) - Exhaust has exactly the same properties as air
%(2) - Gamma = 1.4 and is constant throughout
%(3) - Compression and Expansion are Isentropic
%(4) - All heat is added between the combustion chamber, (2) and (b), and P2=Pb
%(5) - Engine is Adiabatic
%(6) - Shock is infinitely thin, Ax = Ay, hence As encapsulates both
%(7) - Ideal Gas Law is obeyed P=rho*R*T
%(8) - Isentropic everywhere except for shock

%% Miscellaneous Equations Utilized in Derivations
%(1) - U = M*sqrt(gamma*R*T)

%% Declaration of input parameters 
P1 = 70; %freestream pressure (kPa)
T1 = 210; %freestream temperature (K)
M1 = 3.24; %freestream mach number
Mx = 1.2; % Mach number at shock
M2 = 0.3; %Mach number at burner entry
Tb = 1400;%burner temperature (K)
F = 10; %engine thrust (kN)
gamma = 1.4; %Specific Heat Ratio of Air
R = 287; %Ideal Gas Constant
%% Computation of Engine Design Cross-sectional Areas for Selected Inputted Parameters from Tutorial 3
%-Inlet Area A1fixed
%-Inlet Throat Area AC1fixed
%-Burner Entry Area A2fixed
%-Burner Exit Area Abfixed
%-Nozzle Throat Area AC2fixed
%-Exit Area A4fixed
%-Propulsive Efficiency npfixed
%-Thermodynamic Efficiency ncyclefixed
%-Total Efficiency ntfixed
[ntfixed,ncyclefixed,npfixed,~,~,A1fixed,AC1fixed,AC2fixed,A2fixed,Abfixed,A4fixed] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R); 
fprintf('The calculated values of areas are shown below in m^2:\n A1=%f\n AC1=%f\n AC2=%f\n A2=%f\n Ab=%f\n A4=%f\n',A1fixed,AC1fixed,AC2fixed,A2fixed,Abfixed,A4fixed);
fprintf('The calculated values of efficiencies are shown below:\n np=%f\n ncycle=%f\n nt=%f\n\n',npfixed,ncyclefixed,ntfixed);
%% For Clean Code Viewing
disp('This program allows the variation of parameters of:');
disp('1.Pressure');
disp('2.Temperature');
disp('3.Mach Number');
disp('4.Normal Shock Strength');
disp('5.Burner Entry Mach Number'); 
disp('6.Burner Temperature'); 
disp('7.Thrust');
cases = input('Please enter the corresponding number of the parameter which you wish to vary: ');
switch (cases)
case 1
%% ----- (1) Freestream pressure P1 variation ----- %%
initialP1 = P1*10^3; %Initial Freestream Pressure for Variation (Pa)
increP1 = 1000; %Pressure increment (Pa)
finalP1 = 400*10^3; %Final Freestream Pressure for Variation (Pa)
counterP1 = (finalP1-initialP1)/increP1; %Number of test variations of pressure

%-Utilizing the function to compute the varying efficiencies and putting them into matrices corresponding to their freestream pressure-%
%-Preallocating Matrices for Speed
P1mat=zeros(1,counterP1+1);
ntP=zeros(1,counterP1+1);
ncycleP=zeros(1,counterP1+1);
npP=zeros(1,counterP1+1);
for i = 1:1:counterP1+1
    P1 = initialP1 + (i-1)*increP1;
    [ntP(i), ncycleP(i), npP(i)] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R);
    P1mat(i) = initialP1 + (i-1)*increP1;
end

%-Plotting the variation of the efficiencies wrt Freestream Pressure-%
figure
hold on
plot (P1mat/1000,npP,'k-','LineWidth',1.2);
title ('Variation of Efficiencies with Freestream Pressure P_1','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Efficiencies');
Hx=xlabel ('Freestream Pressure P_1 (kPa)');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
plot (P1mat/1000,ncycleP,'k-.','LineWidth',1.2);
hold on
plot (P1mat/1000,ntP,'k:','LineWidth',1.2);
legend('Propulsive efficiency, n_p','Thermodynamic efficiency, n_c_y_c_l_e','Total efficiency, n_t_o_t_a_l','Location','northeastoutside');

 case 2
%% ----- (2) Freestream temperature T1 variation ----- %%
initialT1 = 100; %Initial Freestream Temperature for Variation (K)
increT1 = 1; %Temperature increment (K)
finalT1 = 500; %Final Freestream Temperature for Variation (K)
counterT1 = (finalT1-initialT1)/increT1; %Number of test variations of Temperature

%-Utilizing the function to compute the varying efficiencies and putting them into matrices corresponding to their freestream temperature-%
T1mat=zeros(1,counterT1+1);
for i = 1:1:counterT1+1
    T1 = initialT1 + (i-1)*increT1;
    [ntT(i), ncycleT(i), npT(i), T2T(i),Mbnew(i)] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R);
    T1mat(i) = initialT1 + (i-1)*increT1;
end

%-Plotting the variation of the efficiencies wrt Freestream Temperature-%
figure
hold on
plot (T1mat,npT,'k-','LineWidth',1.2);
title ('Variation of Efficiencies with Freestream Temperature T_1','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Efficiencies');
Hx=xlabel ('Freestream Temperature T_1 (K)');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
plot (T1mat,ncycleT,'k-.','LineWidth',1.2);
hold on
plot (T1mat,ntT,'k:','LineWidth',1.2);
xlim([150,500]);
ylim([0,1.5]);

%-Plotting of Zone of 'Physical Impossibility'-%
for i=1:0.02:2.5 %for Propulsive efficiency higher than 1
    hold on
    yline(i,'r-');
end
for j=100:2:183 %for temperature ratio condition violated, (T2/Tb)*(M2+1/(gamma*M2))^2-(4/gamma) > 0
    hold on
    xline(j,'r-');
end
for k=455:2:500 %Burner Entry temperature would exceed Burner Exit Temperature
    xline(k,'r-');
end
text(200,1.1,'$Physically\ Impossible$','interpreter','latex','FontName','Times New Roman','FontSize',14);
h = text(160,0.9,'$Physically\ Impossible$','interpreter','latex','FontName','Times New Roman','FontSize',14);
set(h,'Rotation',90); %Rotate Text 
h3 = text(475,0.1,'$Unrealistic$','interpreter','latex','FontName','Times New Roman','FontSize',14);
set(h3,'Rotation',90); %Rotate Text 
xline(210); %Chosen Operating Condition
xline(290); %Sea Level Condition
h1 = text(285,0.1,'$Sea\ Level\ Condition$','interpreter','latex','FontName','Times New Roman','FontSize',10);
set(h1,'Rotation',90); %Rotate Text 
h2 = text(205,0.1,'$Chosen\ Operating\ Condition$','interpreter','latex','FontName','Times New Roman','FontSize',10);
set(h2,'Rotation',90); %Rotate Text
legend('Propulsive efficiency, n_p','Thermodynamic efficiency, n_c_y_c_l_e','Total efficiency, n_t','Location','North');

%-Check to ensure engine does not melt under excessive burner entry
%temperature-%
figure
hold on
title ('Variation of Burner Entry Temperature, T_2 with Freestream Temperature T_1','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Burner Entry Temperature T_2 (K)');
Hx=xlabel ('Freestream Temperature T_1 (K)');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
plot(T1mat,T2T);
ylim([0 2100]);
yline(1500,'k:');
text(350,1550,'$\downarrow\ Melting\ Point\ Lower\ limit$','interpreter','latex','FontName','Times New Roman');
yline(2000,'k-.');
text(350,1950,'$\uparrow\ Melting\ Point\ Higher\ limit$','interpreter','latex','FontName','Times New Roman');
yline(1400,'k--');
text(300,1350,'$\uparrow\ Defined\ Burner\ Exit\ Temperature$','interpreter','latex','FontName','Times New Roman');
xline(210); %Chosen Operating Condition
xline(290); %Sea Level Condition
h1 = text(285,200,'$Sea\ Level\ Condition$','interpreter','latex','FontName','Times New Roman','FontSize',10);
set(h1,'Rotation',90); %Rotate Text 
h2 = text(205,200,'$Chosen\ Operating\ Condition$','interpreter','latex','FontName','Times New Roman','FontSize',10);
set(h2,'Rotation',90); %Rotate Text 
legend('Variation of T_2 with T_1','Lower limit of metal melting point','Upper limit of metal melting point','Location','Southeast');

case 3
%% ----- (3) Freestream Mach Number M1 variation ----- %%
%-Notes: Based on Lecture 3, from the effect of burner temperature slides,
%it is understood that the engine would melt at M1=6 and beyond, therefore,
%the maximum limit of the freestream mach number variation would be 6.5 and
%the minimum limit is M1 = 1 as it is a ramjet and is built for supersonic
%flight

initialM1 = 1.0;
increM1 = 0.05;
finalM1 = 6.5;
counterM1 = (finalM1-initialM1)/increM1;

%-Utilizing the function to compute the varying efficiencies and putting them into matrices corresponding to their freestream mach number-%
M1mat=zeros(1,counterM1+1);
for i = 1:1:counterM1+1
    M1 = initialM1 + (i-1)*increM1;
    % '~' is used to skip the outputs which are unnecessary for this analysis
    [ntM(i), ncycleM(i), npM(i), T2M(i),~,A1M(i),~,~,~,~,~,~,~,~,T4M(i)] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R); 
    M1mat(i) = initialM1 + (i-1)*increM1;
end

%-Plotting the variation of the efficiencies wrt Freestream Mach Number-%
figure
hold on
plot (M1mat,npM,'k-','LineWidth',1.2);
title ('Variation of Efficiencies with Freestream Mach Number','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Efficiencies');
Hx=xlabel ('Freestream Mach Number ');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
plot (M1mat,ncycleM,'k-.','LineWidth',1.2);
hold on
plot (M1mat,ntM,'k:','LineWidth',1.2);
%-Plotting of Zone of 'Physical Impossibility'-%
for i=1:0.02:1.2 %for temperature ratio condition violated, (T2/Tb)*(M2+1/(gamma*M2))^2-(4/gamma) > 0
    hold on
    yline(i,'r-');
end
for i=1:0.05:2.9 %for propulsive efficiency higher than 1
    hold on
    xline(i,'r-');
end
for i=5.3:0.05:7 %Burner Entry temperature would exceed Burner Exit Temperature
    hold on
    xline(i,'r-');
end
legend('Propulsive efficiency, n_p','Thermodynamic efficiency, n_c_y_c_l_e','Total efficiency, n_t','Location','south');
text(3,1.1,'$Physically\ Impossible$','interpreter','latex','FontName','Times New Roman','FontSize',14);
h = text(1.5,0.7,'$Physically\ Impossible$','interpreter','latex','FontName','Times New Roman','FontSize',14);
set(h,'Rotation',90); %Rotate Text 
h1 = text(6.5,0.05,'$Physically\ Impossible$','interpreter','latex','FontName','Times New Roman','FontSize',14);
set(h1,'Rotation',90); %Rotate Text 
figure
hold on
plot(M1mat,T2M,'k-');
title ('Burner Entry Temperature T_2 with Freestream Mach Number','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Burner Entry Temperature T_2');
Hx=xlabel ('Freestream Mach Number ');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
yline(1500,'k:');
text(2,1550,'$\downarrow\ Melting\ Point\ Lower\ limit$','interpreter','latex','FontName','Times New Roman');
yline(2000,'k-.');
text(2,1950,'$\uparrow\ Melting\ Point\ Higher\ limit$','interpreter','latex','FontName','Times New Roman');
yline(1400,'k--');
text(2,1350,'$\uparrow\ Defined\ Burner\ Exit\ Temperature$','interpreter','latex','FontName','Times New Roman');
figure

%-Variation of Inlet Area with Freestream Mach Number
plot (M1mat,A1M,'k-');
title ('Variation of Inlet Area A_{1}','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Inlet Area  A_{1} (m^2)');
Hx=xlabel ('Freestream Mach Number M_1');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
xlim([3, 5]);
%-Variation of Temperature ratio T4/T1 with freestream mach number
figure
plot (M1mat,T4M/T1,'k-');
title ('Variation of Temperature Ratio T_4/T_1','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Temperature Ratio T_4/T_1');
Hx=xlabel ('Freestream Mach Number M_1');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
xlim([3, 5]);

case 4
%% ----- (4) Normal Shock Strength Mx Variation ----- %%
%- A Normal shock is required to bring the flow to a subsonic speed in the
%- combustion chamber. It is also favourable to have a weak shock so that
%- there is minimum stagnation pressure loss
%- Vary Normal Shock Strength from 1.0 to 4.0
initialMx = 1.0;
increMx = 0.05;
finalMx = 4.0;
counterMx = (finalMx-initialMx)/increMx;

%-Utilizing the function to compute the varying efficiencies and putting them into matrices corresponding to their Normal Shock Strength-%
Mxmat=zeros(1,counterMx+1);
for i = 1:1:counterMx+1
    Mx = initialMx + (i-1)*increMx;
    [ntMx(i), ncycleMx(i), npMx(i), T2Mx(i),~,A1Mx(i),~,~,~,~,~,M4Mx(i),U4Mx(i),U1Mx(i)] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R);
    Mxmat(i) = initialMx + (i-1)*increMx;
end

%-Plotting the variation of the efficiencies wrt Normal Shock Strength-%
figure
hold on
plot (Mxmat,npMx,'k-','LineWidth',1.2);
title ('Variation of Efficiencies with Normal Shock Strength','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Efficiencies');
Hx=xlabel ('Normal Shock Strength, M_x');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
plot (Mxmat,ncycleMx,'k-.','LineWidth',1.2);
hold on
plot (Mxmat,ntMx,'k:','LineWidth',1.2);
legend('Propulsive efficiency, n_p','Thermodynamic efficiency, n_c_y_c_l_e','Total efficiency, n_t_o_t_a_l','Location','East');

%Plot to find the max shock strength that can be physically possible
figure
hold on
plot(Mxmat,U4Mx.^2-U1Mx.^2,'k-');
title ('(U_4)^2-(U_1)^2 against Normal Shock Strength M_x','fontname','times new roman','fontsize',titl);
Hy=ylabel ('(U_4)^2-(U_1)^2');
Hx=xlabel ('Normal Shock Strength, M_x ');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
yline(0,'k--');
xline(3.375,'k--');
case 5
%% ----- (5) Burner Entry Mach Number M2 Variation ----- %%        
%- The flow within the combustion chamber has to be subsonic, hence the
%- variation was performed for 0.1 to 1.0

initialM2 = 0.1;
increM2 = 0.01;
finalM2 = 1.0;
counterM2 = (finalM2-initialM2)/increM2;

%-Utilizing the function to compute the varying efficiencies and putting them into matrices corresponding to their Burner Entry Mach Number-%
M2mat=zeros(1,counterM2+1);
for i = 1:1:counterM2+1
    M2 = initialM2 + (i-1)*increM2;
    [ntM2(i), ncycleM2(i), npM2(i),T2M2(i),~,~,~,~,~,~,~,~,~,~,T4M2(i)] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R);
    M2mat(i) = initialM2 + (i-1)*increM2;
end

%-Plotting the variation of the efficiencies wrt Burner Entry mach Number%

figure
hold on
plot (M2mat,npM2,'k-','LineWidth',1.2);
title ('Variation of Efficiencies with Burner Entry Mach Number, M_2','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Efficiencies');
Hx=xlabel ('Burner Entry Mach Number, M_2');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
plot (M2mat,ncycleM2,'k-.','LineWidth',1.2);
hold on
plot (M2mat,ntM2,'k:','LineWidth',1.2);
%-Plotting of Zone of 'Physical Impossibility'-%
for i=0.33:0.0035:1 %Burner Entry temperature would exceed Burner Exit Temperature
    hold on
    xline(i,'r-');
end
h1 = text(0.355,0.55,'$Physically\ Impossible$','interpreter','latex','FontName','Times New Roman','FontSize',14);
set(h1,'Rotation',90); %Rotate Text 
xlim([0.1, 0.4]);
legend('Propulsive efficiency, n_p','Thermodynamic efficiency, n_c_y_c_l_e','Total efficiency, n_t_o_t_a_l','Location','NorthEastOutside');

figure
hold on
plot(M2mat,T4M2,'k-');
title ('Variation of T_4 with Burner Entry Mach Number, M_2','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Temperature (K)');
Hx=xlabel ('Burner Entry Mach Number, M_2');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
legend('Temperature at Exhaust , T_4','Location', 'NorthEast');
xlim([0.1,0.35]);
figure
hold on
plot(M2mat,T2M2,'k-');
xlim([0.1,0.35]);
title ('Variation of T_2 with Burner Entry Mach Number, M_2','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Temperature (K)');
Hx=xlabel ('Burner Entry Mach Number, M_2');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
legend('Temperature at Exhaust , T_2','Location', 'NorthEast');
case 6
%% ----- (6) Burner Exit Temperature, Tb Variation ----- %%        
%- Temperature range to observe effect after 1610K
initialTb = 1200;
increTb = 1;
finalTb = 1700;
counterTb = (finalTb-initialTb)/increTb;

%-Utilizing the function to compute the varying efficiencies and putting them into matrices corresponding to their Burner Temperature-%
Tbmat=zeros(1,counterTb+1);
for i = 1:1:counterTb+1
    Tb = initialTb + (i-1)*increTb;
    [ntTb(i), ncycleTb(i), npTb(i), T2Tb(i),Mbnewb(i),~,~,~,~,~,~,~,~,~,T4Tb(i)] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R);
    Tbmat(i) = initialTb + (i-1)*increTb;
end

%-Plotting the variation of the efficiencies wrt Burner Exit Temperature-%
figure
hold on
plot (Tbmat,npTb,'k-','LineWidth',1.2);
title ('Variation of Efficiencies with Burner Exit Temperature, T_b ','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Efficiencies');
Hx=xlabel ('Burner Temperature, T_b');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
plot (Tbmat,ncycleTb,'k-.','LineWidth',1.2);
hold on
plot (Tbmat,ntTb,'k:','LineWidth',1.2);
%-Plotting of Zone of 'Physical Impossibility'-%
for i=1609:2:1700 %Imaginary Roots 
    hold on
    xline(i,'r-');
end    
h1 = text(1650,0.5,'$Physically\ Impossible$','interpreter','latex','FontName','Times New Roman','FontSize',14);
set(h1,'Rotation',90); %Rotate Text     
legend('Propulsive efficiency, n_p','Thermodynamic efficiency, n_c_y_c_l_e','Total efficiency, n_t','Location','southwest');

case 7
%% ----- (7) Thrust, F Variation ----- %%          
%- Thrust was allowed to vary from 1kN to 200kN

initialF = 1*10^3;
increF = 1;
finalF = 200*10^3;
counterF = (finalF-initialF)/increF;

%-Utilizing the function to compute the varying efficiencies and putting them into matrices corresponding to their Burner Temperature-%
Fmat=zeros(1,counterF+1);
for i = 1:1:counterF+1
    F = initialF + (i-1)*increF;
    [ntF(i), ncycleF(i), npF(i), T2F(i)] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R);
    Fmat(i) = initialF + (i-1)*increF;
end

%-Plotting the variation of the efficiencies wrt Thrust-%
figure
hold on
plot (Fmat,npF,'k-','LineWidth',1.2);
title ('Variation of Efficiencies with Thrust','fontname','times new roman','fontsize',titl);
Hy=ylabel ('Efficiencies');
Hx=xlabel ('Thrust, F');
set(Hy, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
set(Hx, 'fontname', 'times new roman', 'fontsize', font, 'rotation', 0);
hold on
plot (Fmat,ncycleF,'k-.','LineWidth',1.2);
hold on
plot (Fmat,ntF,'k:','LineWidth',1.2);
legend('Propulsive efficiency, n_p','Thermodynamic efficiency, n_c_y_c_l_e','Total efficiency, n_t_o_t_a_l','Location','northeastoutside');

end

%% Design Code
%%- Able to provide Output of -%
%%- Total Efficiency, Thermodynamic Efficiency, Propulsive Efficiency,
%%- Burner Exit Mach, Inlet Area, Inlet Throat Area, Nozzle Throat Area, 
%%- Burner Entry Area, Burner Exit Area, Exhaust Area, Exhaust Mach,
%%- Velocity at Exhaust, Velocity at Inlet, Temperature at Exhaust
%%- Respectively
 function [nt,ncycle,np,T2,Mbnew,A1,AC1,AC2,A2,Ab,A4,M4,U4,U1,T4] = design(P1,T1,M1,Mx,M2,Tb,F,gamma,R)

%% -----Analysis at Inlet Station (1)----- %%
%---Calculation of Isentropic Properties---%
%-Freestream Stagnation Temperature/Freestream Static Temperature, T01/T1-%
tempratio1 = 1 + ((gamma-1)/2).*(M1.^2); 

%-Freestream Stagnation Pressure/Freestream Static Pressure, P01/P1-%
presratio1 = tempratio1.^(gamma/(gamma-1)); 

%-Area of Inlet/Area of Inlet Throat, A1/AC1-%
arearatio1 = (1/M1)*((2/(gamma+1))*(tempratio1)).^((gamma+1)/...
    (2*(gamma-1))); 

%% -----Analysis at Inlet Throat Station (C1)----- %%
%-Area of Inlet Throat/Area of Inlet, AC1/A1-%
arearatioC1 = 1 / arearatio1; 

%% -----Analysis at Normal Shock Station (x,y)----- %%
%Note: Stagnation Conditions change across the shock; Non-isentropic
%-Mach Number immediately downstream of Shock, My-%
My = sqrt(((gamma-1)*Mx.^2 + 2)/(2*gamma*Mx.^2 - (gamma-1)));

%-Static Pressure After Shock/Static Pressure Before Shock, Py/Px-%
presratioshock = (2*gamma*Mx.^2 - (gamma-1))/(gamma+1);

%-Density After Shock/Density Before Shock, rhoy/rhox-%
densityratioshock = ((gamma+1)*Mx.^2)/((gamma-1)*Mx.^2 + 2);

%%-Stagnation Conditions at (x)-%%
%-Stagnation Temperature/Static Temperature, T0x/Tx-%
tempratiox = 1+((gamma-1)/2)*Mx.^2;
%-Stagnation Pressure/Static Pressure, P0x/Px-%
presratiox = (tempratiox)^(gamma/(gamma-1));
%-Area of Shock/Area of Inlet Throat, Ax/AC1 = As/AC1-%
arearatiox = (1/Mx)*((2/(gamma+1))*tempratiox)^((gamma+1)/(2*(gamma-1)));

%%-Stagnation Conditions at (y)-%%
%-Stagnation Temperature/Static Temperature, T0y/Ty-%
tempratioy = 1+((gamma-1)/2)*My.^2;
%-Stagnation Pressure/Static Pressure, P0y/Py-%
presratioy = (tempratioy)^(gamma/(gamma-1));
%-Area of Shock/Area of Inlet Throat, Ay/AC1 = As/AC1 = As/Ay*-%
arearatioy = (1/My)*((2/(gamma+1))*tempratioy)^((gamma+1)/(2*(gamma-1)));

%% -----Analysis at Burner Entry (2)----- %%
%-Stagnation Temperature/Static Temperature, T02/T2-%
tempratio2 = 1+((gamma-1)/2)*M2^2;
%-Stagnation Pressure/Static Pressure, P02/P2-%
presratio2 = (tempratio2)^(gamma/(gamma-1));
%-Area of Burner Entry/Area of Inlet Throat, A2/A2* = A2/Ay*-%
arearatio2 = (1/M2)*((2/(gamma+1))*tempratio2)^((gamma+1)/(2*(gamma-1)));

%-Burner Entry Area/Inlet Area, A2/A1-%
%-Note: Isentropic, A2*/Ay* = A1*/Ax* = 1
%-Note: Due to infinitely thin shock, Ay/Ax = 1
%-A2/A1 = (A2/A2*)*(A2*/Ay*)*(Ay*/Ay)*(Ay/Ax)*(Ax/Ax*)*(Ax*/A1*)*(A1*/A1)
arearatio21 = (arearatio2)*(arearatioy^(-1))*(arearatiox)*(arearatio1^(-1));

%-Burner Entry Temperature/Inlet Temperature, T2/T1-%
%-Note: Isentropic, T02/T0y = T01/T0x
%-Note: Across Shock assumed to be adiabatic, T0y=T0x
%-Burner Entry Temperature/Inlet Temperature, T2/T1-%
tempratio21 = ((tempratio2)^-1)*tempratio1;

%Computation of Burner Entry Temperature, T2 for subsequent analysis
T2 = tempratio21*T1;

%% -----Analysis from Burner Entry to End of Burner (2,b)----- %%
%Note: Apply Conservation of Mass and Conservation of Momentum
%-Conservation of Mass, rho2*U2*A2 = rhob*Ub*Ab-%
%(P2*A2)/(Pb*Ab) = (Mb/M2)*sqrt(T2/Tb)

%-Conservation of Momentum, rho2*U2^2*A2 + P2A2 = rhob*Ub^2*Ab + PbAb
%(P2*A2)/(Pb*Ab) = (gamma*Mb^2+1)/(gamma*M2^2)+1

%-Combination of both expressions-%
%(Mb/M2)*sqrt(T2/Tb) = (gamma*MB^2+1)/(gamma*M2^2+1)
coeffMbsq = gamma*M2*sqrt(Tb/T2);%(gamma)/(gamma*M2^2+1);
coeffMb = -(gamma*M2^2+1);%-sqrt(T2/Tb)/M2;
coeffconstant = M2*sqrt(Tb/T2); %1/(gamma*M2^2+1);

eqn = [coeffMbsq coeffMb coeffconstant];
Mb = transpose(roots(eqn));
%-Supersonic solution is undesired, take subsonic Mb-%
  for i=1:1:2
      if Mb(i) < 1      
          Mbnew =Mb(i);
      end
  end

%-Since the burning process is isobaric, Pb=P2-%
%-(P2*A2)/(Pb/Ab) = A2/Ab = (gamma*Mb^2+1)/(gamma*M2^2+1)-%
arearatio2b = (gamma*Mbnew^2+1)/(gamma*M2^2+1);

%-To calculate Ab/A1, Ab/A1 = (Ab/A2)*(A2/A1)-%
arearatiob1 = (arearatio2b)^(-1)*arearatio21;

%-New Stagnation Conditions at Burner Exit (b)-%
%-Stagnation Temperature at (b)/Static Temperature at (b), T0b/Tb-%
tempratiob = (1+((gamma-1)/2)*Mbnew^2);
%-Stagnation Pressure at (b)/Static Pressure at (b), P0b/Pb-%
presratiob = tempratiob^(gamma/(gamma-1));
%-Area of Burner Exit/Area of Throat, Ab/Ab* = Ab/AC2-%
arearatiob = (1/Mbnew)*((2/(gamma+1))*(tempratiob))^((gamma+1)/(2*(gamma-1)));

%-Nozzle Throat Area/Area of Inlet, Ab*/A1 = AC2/A1-%
%-Ab*/A1 = (Ab*/Ab)*(Ab/A1)-%
arearatiobstar1 = arearatiob^(-1)*arearatiob1;

%% -----Analysis from Burner Exit to Exhaust (b,4)----- %%
%-Stagnation Pressure at (4)/Static Pressure at (4), P04/P4-%
%- P04/P4 =
%- (P04/P0b)*(P0b/Pb)*(Pb/P2)*(P2/P02)*(P02/P0y)*(P0y/Py)*(Py/Px)*(Px/P0x)*(P0x/P01)*(P01/P1)*(P1/P4)-%
%-Note: Isentropic Flow, P04/P0b = P02/P0y = P0x/P01 = 1
%-Note: By Design, Pb/P2 = P1/P4 = 1
presratio4 = presratiob*((presratio2)^(-1))*presratioy*presratioshock*((presratiox)^(-1))*presratio1;

%-Mach Number at Exit, M4 for P04/P4-%
M4 = sqrt(((presratio4^((gamma-1)/gamma))-1)/((gamma-1)/2));
%-Stagnation Conditions at Exit (4)-%

%-Stagnation Temperature at (4)/Static Temperature at (4),T04/T4-%
tempratio4 = 1 + ((gamma-1)/2)*M4^2;

%-Stagnation Pressure at (4)/Static Pressure at (4), P04/P4-%
presratio4calc = tempratio4^(gamma/(gamma-1)); 

%-Area of Exit/Area of Nozzle Throat, A4/AC2 = A4/A4*-%
arearatio4 = (1/M4)*((2/(gamma+1))*(tempratio4))^((gamma+1)/(2*(gamma-1)));

%-Area of Exhaust/Area of Inlet, A4/A1-%
%-A4/A1 = (A4/A4*)*(A4*/Ab*)*(Ab*/A1)-%
%-Note: Isentropic flow from (b) to (4) hence, A4*/Ab* = 1
arearatio41 = arearatio4*arearatiobstar1;


%% -----Thrust and Cross-sectional Area Calculation----- %%
%- F ~ rho4*U4^2*A4 - rho1*U1^2*A1 + P4*A4 - P1*A1
%- P4*A4 - P1*A1 is 0 due to control volume encompassing the whole engine-%
%- F/(P1*A1) = gamma*M1^2*((M4/M1)^2*(A4/A1)-1)
thrustPA = gamma*M1^2*((M4/M1)^2*arearatio41-1);

%-Inlet Area, A1-%
A1 = F/(thrustPA*P1);

%-Inlet Throat Area, AC1-%
AC1 = arearatioC1 * A1;

%-Burner Entry Area, A2-%
A2 = arearatio21 * A1;

%-Burner Exit Area, Ab-%
Ab = arearatiob1 * A1;

%-Nozzle Throat Area, AC2-%
AC2 = arearatiobstar1 * A1;

%-Exit Area, A4-%
A4 = arearatio41 * A1;

%% -----Efficiency Calculation----- %%
%Temperature at Exit (4)/Temperature at Inlet (1), T4/T1-%
tempratio41 = tempratio4^(-1)*tempratiob*(Tb/T1);
%-Propulsive Efficiency-%
np = (F/(P1*A1))*(2/gamma)*(1/((tempratio41*M4^2)-M1^2));
%-Thermodynamic Efficiency-%
%-Note:Assuming perfect efficiency of compressor and turbine
ncycle = 1 - (tempratio41*T1 - T1)/(Tb-T2);%1- tempratio21^(-1);
%-Total Efficiency-%
nt = ncycle .* np;

%-Calculation of Parameters for Post-Processing
T4 = tempratio41*T1; %temperature at exhaust
a4 = sqrt(gamma*R*T4); %speed of sound at exhaust
a1 = sqrt(gamma*R*T1); %speed of sound at inlet
U4 = (M4^2)*a4; %velocity at exhaust
U1 = (M1^2)*a1; %velocity at inlet
end
