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
font=18;
titl=20;
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
Mx = 1.3; % Mach number at shock
M2 = 0.4; %Mach number at burner entry
Tb = 1400;%burner temperature (K)
F = 10; %engine thrust (kN)
gamma = 1.4; %Specific Heat Ratio of Air



%% Design Code
 %function [nt,np,ncycle,T2] = design(P1,T1,M1,Mx,M2,Tb,F,gamma)

%% -----Analysis at Inlet Station (1)----- %%
%---Calculation of Isentropic Properties---%
%-Freestream Stagnation Temperature/Freestream Static Temperature, T01/T1-%
tempratio1 = 1 + ((gamma-1)/2)*(M1^2); 

%-Freestream Stagnation Pressure/Freestream Static Pressure, P01/P1-%
presratio1 = tempratio1^(gamma/(gamma-1)); 

%-Area of Inlet/Area of Inlet Throat, A1/AC1-%
arearatio1 = (1/M1)*((2/(gamma+1))*(tempratio1))^((gamma+1)/...
    (2*(gamma-1))); 

%% -----Analysis at Inlet Throat Station (C1)----- %%
%-Area of Inlet Throat/Area of Inlet, AC1/A1-%
arearatioC1 = 1 / arearatio1; 

%% -----Analysis at Normal Shock Station (x,y)----- %%
%Note: Stagnation Conditions change across the shock; Non-isentropic
%-Mach Number immediately downstream of Shock, My-%
My = sqrt(((gamma-1)*Mx^2 + 2)/(2*gamma*Mx^2 - (gamma-1)));

%-Static Pressure After Shock/Static Pressure Before Shock, Py/Px-%
presratioshock = (2*gamma*Mx^2 - (gamma-1))/(gamma+1);

%-Density After Shock/Density Before Shock, rhoy/rhox-%
densityratioshock = ((gamma+1)*Mx^2)/((gamma-1)*Mx^2 + 2);

%%-Stagnation Conditions at (x)-%%
%-Stagnation Temperature/Static Temperature, T0x/Tx-%
tempratiox = 1+((gamma-1)/2)*Mx^2;
%-Stagnation Pressure/Static Pressure, P0x/Px-%
presratiox = (tempratiox)^(gamma/(gamma-1));
%-Area of Shock/Area of Inlet Throat, Ax/AC1 = As/AC1-%
arearatiox = (1/Mx)*((2/(gamma+1))*tempratiox)^((gamma+1)/(2*(gamma-1)));

%%-Stagnation Conditions at (y)-%%
%-Stagnation Temperature/Static Temperature, T0y/Ty-%
tempratioy = 1+((gamma-1)/2)*My^2;
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
%Note: Isentropic, A2*/Ay* = A1*/Ax* = 1...
%Note: Due to infinitely thin shock, Ay/Ax = 1
%A2/A1 = (A2/A2*)*(A2*/Ay*)*(Ay*/Ay)*(Ay/Ax)*(Ax/Ax*)*(Ax*/A1*)*(A1*/A1)
arearatio21 = (arearatio2)*(arearatioy^(-1))*(arearatiox)*(arearatio1^(-1));

%-Burner Entry Temperature/Inlet Temperature, T2/T1-%
%Note: Isentropic, T02/T0y = T01/T0x
%Note: Across Shock assumed to be adiabatic, T0y=T0x
%Burner Entry Temperature/Inlet Temperature, T2/T1-%
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
%-Supersonic solution is impossible, take subsonic Mb-%
  for i=1:1:2
      if Mb(i) < 1      
          Mbnew = real(Mb(i));
      end
  end
%Mb(1)=[];
%Mb=real(Mb);

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


%% -----Thrust Calculation----- %%
%- F ~ rho4*U4^2*A4 - rho1*U1^2*A1 + P4*A4 - P1*A1
%- P4*A4 - P1*A1 is 0 due to control volume encompassing the whole engine-%
%- F/(P1*A1) = gamma*M1^2*((M4/M1)^2*(A4/A1)-1)
thrustPA = gamma*M1^2*((M4/M1)^2*arearatio41-1);

%-Inlet Area, A1-%
A1 = F/(thrustPA*P1); %m^2

%-Inlet Throat Area, AC1-%
AC1 = arearatioC1 * A1; %m^2

%-Burner Entry Area, A2-%
A2 = arearatio21 * A1; %m^2

%-Burner Exit Area, Ab-%
Ab = arearatiob1 * A1; %m^2
 
%-Nozzle Throat Area, AC2-%
AC2 = arearatiobstar1 * A1; %m^2

%-Exit Area, A4-%
A4 = arearatio41 * A1; %m^2

%% -----Efficiency Calculation----- %%
%Temperature at Exit (4)/Temperature at Inlet (1), T4/T1-%
tempratio41 = tempratio4^(-1)*tempratiob*(Tb/T1);
%-Propulsive Efficiency-%
np = (F/(P1*A1))*(2/gamma)*(1/((tempratio41*M4^2)-M1^2));
%-Thermodynamic Efficiency-%
%-Note:Assuming perfect efficiency of compressor and turbine
ncycle = 1 - (tempratio41*T1 - T1)/(Tb-T2); %1- tempratio21^(-1);
%-Total Efficiency-%
nt = ncycle.* np;
% end


