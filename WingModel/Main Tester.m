% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

%Flight Condition
clear;clc;
velocity = 62.6;%102m/s -> 0.3M
%velocity = 83.9; %cessna max
%velocity = 205.8; %globalhawk cruise. make sure to increase span
AoA = 5; %degrees

MainFlightCondition = FlightCondition(velocity, AoA);
airfoil_path = "JustAirfoil-NACA0012.xlsx";


%inputs = [5.0 1.5 1.5 0.01 .01 .01]'; %rect wing analysis geometry
%inputs = [20 1.5 0.3 0.008 .005 .02]';
inputs = [8 1.5 0.4 0.0028 .01 .01]'; %rect wing analysis geometry
%doesn't work for root chord of less than 1.5. investigate why later

%%
i=1;
tic;
[total_mass(i), max_disp(i), KS_stress(i), CD(:,i), CL(:,i), CL_approx(i), CD_approx(i)] = WingModeler(MainFlightCondition, airfoil_path, inputs(:,i), 0);
toc;
delete('patran.ses.*');
disp('0-AoA test generated! 1')
%%
AoA = 5; %degrees
MainFlightCondition = FlightCondition(velocity, AoA);
airfoil_path = "JustAirfoil-NACA0012.xlsx";

i=1;
tic;
[total_mass(i), max_disp(i), KS_stress(i), CD(:,i), CL(:,i), CL_approx(i), CD_approx(i)] = WingModeler(MainFlightCondition, airfoil_path, inputs(:,1), i);
toc;
delete('patran.ses.*');
disp('5-AoA test generated! 2')

%%
AoA = 0;
MainFlightCondition = FlightCondition(velocity, AoA);
airfoil_path = "JustAirfoil-NACA2412.xlsx";

i=3;
tic;
[total_mass(i), max_disp(i), KS_stress(i), CD(:,i), CL(:,i), CL_approx(i), CD_approx(i)] = WingModeler(MainFlightCondition, airfoil_path, inputs(:,1));
toc;
delete('patran.ses.*');
disp('0-AoA, NACA2412 test generated! 3')