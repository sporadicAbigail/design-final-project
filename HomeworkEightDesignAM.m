% This code is for CHEN 4520 Chemical Process Synthesis (Fall 2022)
% Homework 8
% Written by: Aimee Maravi, Abigail Hutabarat, Mady Murphy, 
% Sean Oishi-Holder, and Mohammed Alahmed

clear
clc

T_out = 533; %outlet Temperature in kelvin 
Cbd = 1975; %kg cat/m^3 this is catalyst bulk density 
Void = 0.50; %catalyst voidage
%heurisitc to account for later: catalyst particle diameter is an eigth of
%the tube diameter

%% Initial Conditions
%T0 = 1; %what units 

%this is molar flow rates [mol/hr] 
F_Product_in = 0; 
F_C2H4_i= 970240;
F_HCl_i = 2910750;
F_O2_i= 970240;
F_CO2_i = 0;
F_H2O_i = 0;
F_Cl3Eth_i = 0;
F_Cl2_i = 5.820*1000;
%initial molar flow rate  is constant
molFlow_i = F_C2H4_i + F_HCl_i + F_O2_i + F_CO2_i + F_H2O_i + F_Cl3Eth_i + F_Cl2_i; 

%Guess T initial
Tube_temp_in= 350; %assuming the reaction happens at 533K
Coolant_temp_in = 278;
Pressure_i = 2000; %inlet pressure in kPa
IC = [F_Product_in, F_C2H4_i, F_HCl_i, F_O2_i, F_CO2_i, F_H2O_i, F_Cl3Eth_i, F_Cl2_i, Tube_temp_in, Coolant_temp_in, Pressure_i];

%% this is my Domain
V_I = 0;
%guessing the length is of the reactor is z? do we also just guess and
%check this value???????????????/
V_Domain = [V_I 3875]; %Define the temperature domain
%% SOLVE ODE 
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[Vsol, Ysol] = ode15s('FinalTakeOde', V_Domain, IC,opts);

%% Data handling 
% Extraction Solutions for each variable 
Product_sol = Ysol (:,1);
C2H4_sol = Ysol(:,2); 
HCl_sol = Ysol(:,3);
O2_sol = Ysol(:,4);
CO2_sol = Ysol(:,5);
H2O_sol = Ysol(:,6);
Cl3Eth_sol = Ysol(:,7); 
Cl2_sol = Ysol(:,8);
TubeTemp_sol = Ysol(:,9);
CoolantTemp_sol = Ysol(:,10);
Pressure_sol = Ysol(:,11);

%% plotting 

%I realized we need a reactor volume rate

fig1 = figure(1);
hold on
set(fig1,'name', 'Reactor Temperature Profile')
xlabel('Reactor Volume (L)');
ylabel('FlowRates');
plot(Vsol,C2H4_sol,'--')
plot(Vsol,Product_sol,'o')
plot(Vsol,HCl_sol,'b--o')
plot(Vsol,O2_sol,':')
plot(Vsol,H2O_sol)
plot(Vsol,Cl3Eth_sol)
plot(Vsol,Cl2_sol)
plot(Vsol,CO2_sol,'.')
legend('C2H2','Product','HCl','O2','H2O','Cl3eth','Cl2','CO2')
hold off 

fig2 = figure(2);
set(fig2,'name', 'Reactor Temperature Profile')
plot(Vsol,TubeTemp_sol);
xlabel('Tube Volume');
ylabel('Temperature ');

fig3 = figure(3); 
set(fig2,'name', 'Temp')
plot(Vsol,CoolantTemp_sol);
xlabel('Reactor Volume');
ylabel('Coolant Temperature');

fig4 = figure(4);
set(fig4, 'name','pressure')
plot(Vsol,Pressure_sol);
xlabel('Reactor Volume (L)');
ylabel('pressure (kPa)');

Flow_p = mean(Product_sol);
Flow_C2H4 = mean(C2H4_sol);
Flow_HCl = mean(HCl_sol);
Flow_O2 = mean(O2_sol);
Flow_H2O = mean(H2O_sol);
Flow_Cl3Eth = mean(Cl3Eth_sol);
Flow_Cl2 = mean(Cl2_sol);
Flow_CO2 = mean(CO2_sol);

%
%points of confusion for me 1) what the heck is the units of Ysol?
%Ysol is solution vector for all the ODEs
%Xsol will be in units of volume, everything in terms of volume
%are we integrating by 2 variables? No - technically more than one