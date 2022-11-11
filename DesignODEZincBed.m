function [Zinc_ODE] = DesignODEZincBed(V,Y)
%% Reactor Volume Parameters zinc oxide bed
% Use the kinetics to kind of decide on a space time or amount of volume of
% the reactor
% this comes first because everything is going to be done on a per tube basis
reactorLength =10; %meters, CHANGE ME
 %kg/hr
mDotC = 10000; %kg/hr, not sure what this is

%% Other things that change
T = Y(9); %Kelvin
T_celcius = T-273;
T_0 = 483; %Kelvin, initial temp
P = Y(11); %bar
P_0 = 35; %bar
Tc = Y(10); %Kelvin
%ignore the kinetics of ZnO reacting with HCl and just consider ZnO and H2S
%the HCl is negligible???
F_ZnO = Y(1); %zinc oxide (solid)
F_H2S = Y(2); %hydrogen sulfide (gas)
F_ZnS = Y(3); %ZnS (solid)
F_H2O = Y(4); %LIQUID - The vaporization temp at 35 bar is 515.75
F_HCl= Y(5); %gas
F_ZnCl2 =Y(6); %solid
F_tot_gas = F_H2S+F_HCl;

%% Rate Related Constants - This is updated for the zinc bed
gasConst = .008314; %kj/mol-K
%assume ideal gas law and that partial pressures are valid for use of
%concentration
k1 = 179571*exp(-19.26/gasConst/T); %units of 1/hr
rxn1 = k1*F_H2O/F_tot_gas; %units of concentration of H2O/hr

%% Ergun Related Equations - THIS HAS NOT BEEN FIXED
F_tot_0_overall = 970240 + 2910750 + 970240 + 5820; %[mol/hr], constant
F_tot_0 = F_tot_0_overall/numTubes; %[mol/hr] per tube, constant
rho_0 = 18.002; %kg/m3, constant
rho = rho_0*(P/P_0)*(T_0/T)*(F_tot_0/F_tot); %kg/m3, changing
volumetricFlowRate_tot_0 = 166430.8*1000; %L/hr, calculated in writeup, overall, constant
volumetricFlowRate_0 = volumetricFlowRate_tot_0/numTubes; %L/hr, per tube, constant
Tau = 3.6; %space time in seconds, given 
reactorVol = Tau/3600*volumetricFlowRate_tot_0; %3600 to convert from s to hr: [L]
reactorVolm3 = reactorVol/1000; %[m3]
Ac = reactorVolm3/(numTubes*reactorLength); %[m]
tubeDiameter=2*sqrt(Ac/pi); %cross sectional area [m2]
%here, Matlab calculates a reactor length
% intentionally without a semicolon so it comes out in our command window
% as a result
superFicVelocity = volumetricFlowRate_0/1000/Ac/3600; % m/s
G = rho*superFicVelocity; % superficial mass velocity [kg/m2-s]
phi = 0.50; %void fraction/porosity [unitless], constant
mu = 0.02306/1000; %viscosity of gas [Pa-s] or [kg/m-s], constant
gc = 1; %unit conversion [metric], constant
particleDiameter = tubeDiameter/8; %based on heuristic [m], constant
beta_0 = G*(1-phi)/(rho_0*gc*particleDiameter*phi^3)*(150*(1-phi)*mu/particleDiameter+1.75*G)/1000; %to get units of kPa/m

%% Energy Balance Related Equations
C_ZnO = 0.05285;
Cp_ZnO = C_ZnO;

%constants for heat capacity of H2S, where T is in CELSIUS
C_H2S = [0.03351, 1.547*10^(-5), 0.3012*10^(-8), -3.292*10^(-12)];
%Cp has units of kJ/mol DEGREE C
Cp_H2S = C_H2S(1)+C_H2S(2)*T_celcius+C_H2S(3)*T_celcius^2+C_H2S(4)*T_celcius^3;

%constants to calculate for heat capacity of ZnS, Where T is in KELVIN
C_ZnS = 0.045715306;
%Cp has units of kJ/mol degree K
Cp_ZnS = C_ZnS; 

%this is for liquid water, DEGREES C
C_H2O = 0.0754;
Cp_H2O = C_H2O;

%constants to calculate for heat capacity of HCl, where T is in CELCIUS
C_HCl = [0.02913, -0.1341*10^(-5), 0.9715*10^(-8),-4.335*10^(-12)];
Cp_HCl = C_HCl(1)+C_HCl(2)*T_celcius+C_HCl(3)*T_celcius^2+C_HCl(4)*T_celcius^3;


C_ZnCl2 = 0.1007889484;
Cp_ZnCl2 = C_ZnCl2;


%should these be neg or pos
hrxn1 = -239.111; %kj/mol
Ua = 300*tubeDiameter*pi*reactorLength/(reactorVolm3/numTubes)*0.0036; % [kJ/(hr L-cat K) Heat capacity*surface area of heat transfer/volume
%% The ODEs
% mass balance
rC2H4 = -(-rxn1-rxn3);
rProd = -(-rxn2 + rxn1);
rHCl = -(-2*rxn1 - rxn2 - 2*rxn4);
rO2 = -(-0.5*rxn1-0.5*rxn2-3*rxn3-0.5*rxn4);
rCO2 = -(2*rxn3);
rH2O = -(rxn1+rxn2+2*rxn3+rxn4);
rCl3Eth = -rxn2;
rCl2 = -rxn4;
%Ergun equation
rP = -beta_0/(1-phi)/Ac*(P_0/P)*(T/T_0)*(F_tot/F_tot_0); %units of kPa/m/m2, everything else is unitless
%Thermal equation
chekpointA = -rxn1*hrxn1-rxn2*hrxn2-rxn3*hrxn3-rxn4*hrxn4;
checkpointb=Ua*(T-Tc); %goal have the check points be about the same
numerator = (rxn1*hrxn1+rxn2*hrxn2+rxn3*hrxn3+rxn4*hrxn4)-(Ua*(T-Tc));
denominator = (F_C2H4*Cp_C2H4+F_Prod*Cp_Prod+F_HCl*Cp_HCl+F_O2*Cp_O2+F_CO2*Cp_CO2+F_H2O*Cp_H2O+F_Cl3Eth*Cp_Cl3Eth+F_Cl2*Cp_Cl2)*(1-phi);
rT = numerator/denominator;
%cooling thermal balance
rTc = Ua*(T-Tc)/(mDotC/numTubes)/Cp_dowtherm;
%% Convert to original script
dC2H4 = rC2H4;
dProd = rProd;
dHCl = rHCl;
dO2 = rO2;
dCO2 = rCO2;
dH2O = rH2O;
dCl3Eth = rCl3Eth;
dCl2 = rCl2;
dP = rP;
dT = rT;
dTc = rTc;
Zinc_ODE = [dProd; dC2H4; dHCl; dO2; dCO2; dH2O; dCl3Eth; dCl2; dT; dTc; dP];
% %Goal
% F_C2H4_fin = 0.04*970240/numTubes;
% if F_C2H4 <= F_C2H4_fin 
%     return;
% end
end 