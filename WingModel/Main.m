% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

clear; clc; close all; 

fprintf('Aerospace_MDO  Copyright (C) 2020  Adam Benabbou\n')
fprintf('This program comes with ABSOLUTELY NO WARRANTY\n')
fprintf('This is free software, and you are welcome to redistribute it\n')
fprintf('under certain conditions\n')


%62.6 m/s for cessna cruise
%83.9 m/s for cessna max

velocity = 62.6;%102m/s -> 0.3M
AoA = 5; %degrees

MainFlightCondition = FlightCondition(velocity, AoA);

%INPUTS  -- sample
span = 6; %meters - const.
chord_root = 1.0; % meters - VAR.
chord_tip = .40; % meters - VAR.
panel_t1 = 0.01; % meters - VAR.

rib_t = .0051; % meters - const.lb
spar_t = .0251; % meters - const.

%variable order: span,chord_root, chord_tip, panel_t1, rib_t, spar_t

lb = [ 5.0   1.5   0.3   0.003  0.005  0.015 ]';
ub = [ 12.0   1.5   1.5   0.02  0.005  0.015]';

N = 1000;
n = length(lb);
Q = floor(N/4);
xn = lhsdesign(n,N,'criterion','maximin','iterations',100); %normalized LHS design
inputs = bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb))); %LHS on the problem design space domain

%airfoil path 
airfoil_path = "JustAirfoil-NACA0012.xlsx";

%preallocation
total_mass = zeros(1,length(inputs));
max_disp = zeros(1,length(inputs));
KS_stress = zeros(1,length(inputs));
CD = zeros(2,length(inputs));
CL = zeros(2,length(inputs));
CD_approx = zeros(1,length(inputs));
CL_approx = zeros(1,length(inputs));

%% Generate Samples
tic
parfor i=1:1000
[total_mass(i), max_disp(i), KS_stress(i), CD(:,i), CL(:,i), CL_approx(i), CD_approx(i)] = WingModeler(MainFlightCondition, airfoil_path, inputs(:,i),i);
end
toc

disp('==========( DONE )==========')
delete('patran.ses.*');

%% Generate Data Matrix to Visualize Sampling and Design Space
figure
data = [inputs(1,:)' inputs(3,:)' inputs(4,:)' total_mass(:) max_disp(:) KS_stress(:) CL(2,:)' CD(2,:)'];
labels = {'Span(m)' 'Tip Chord(m)' 'Panel Thickn.(m)' 'Total Mass(kg)' 'Tip Disp.(m)' 'KS Stress(Pa)' 'CL' 'CD'};
[S, AX, BigAx, H,HAx] = plotmatrix(data, '*');                        % create a 4 x 4 matrix of plots
for i = 1:8                                       % label the plots
  xlabel(AX(8,i), labels{i}, 'FontSize', 14)
  ylabel(AX(i,1), labels{i}, 'FontSize', 10)
end

for i=1:3
    for j = 1:3
        S(i,j).Color = [.4 0.247 0.741];
    end
end


for i=4:8
    for j = 4:8
        S(i,j).Color = [.4 0.747 0.741];
    end
end

for i=4:8
    for j = 1:3
        S(i,j).Color = [.7 0.247 0.341];
    end
end
title('DataMatrix of Sample Designs');

%% Test Different Spread Values for RBN
goal_mse_mass = (1e-4*range(total_mass));
goal_mse_disp = (1e-4*range(max_disp));
goal_mse_stress = (1e-4*range(KS_stress));
goal_mse_CL = (1e-4*range(CL(2,:)));
goal_mse_CD = (1e-4*range(CD(2,:)));

goal_mse_CL_approx = (1e-4*range(CL_approx));
 
spreads = logspace(-2,1,20);

parfor i = 1:length(spreads)
     net_mass = newrb(inputs,total_mass,goal_mse_mass,spreads(i),Q,1);
     y_sim_mass(i,:) = sim(net_mass,inputs);
     mse_mass(i) = mse(net_mass,total_mass,y_sim_mass(i,:));
end

parfor i = 1:length(spreads)
     net_disp = newrb(inputs,max_disp,goal_mse_disp,spreads(i),Q,1);
     y_sim_disp(i,:) = sim(net_disp,inputs);
     mse_disp(i) = mse(net_disp,max_disp,y_sim_disp(i,:));
end

parfor i = 1:length(spreads)
     net_stress = newrb(inputs,KS_stress,goal_mse_stress,spreads(i),Q,1);
     y_sim_stress(i,:) = sim(net_stress,inputs);
     mse_stress(i) = mse(net_stress,KS_stress,y_sim_stress(i,:));
end

parfor i = 1:length(spreads)
     net_CL_approx = newrb(inputs,CL_approx,goal_mse_CL_approx,spreads(i),Q,1);
     y_sim_CL_approx(i,:) = sim(net_CL_approx,inputs);
     mse_CL_approx(i) = mse(net_CL_approx,CL_approx,y_sim_CL_approx(i,:));
end

parfor i = 1:length(spreads)
     net_CL = newrb(inputs,CL(2,:),goal_mse_CL,spreads(i),Q,1);
     y_sim_CL(i,:) = sim(net_CL,inputs);
     mse_CL(i) = mse(net_CL,CL(2,:),y_sim_CL(i,:));
end

parfor i = 1:length(spreads)
     net_CD = newrb(inputs,CD(2,:),goal_mse_CD,spreads(i),Q,1);
     y_sim_CD(i,:) = sim(net_CD,inputs);
     mse_CD(i) = mse(net_CD,CD(2,:),y_sim_CD(i,:));
end
%% Visualize best Spread for RBN
figure (1)
hold on
plot(spreads, mse_mass(1:length(spreads)), 'o-')
title(['Mean Square Error with 250 Neurons or less (Total Mass)'])
xlabel(['Spread Value'])
ylabel(['Mean Square Error'])
grid on
hold off

figure(2)
hold on
plot(spreads, mse_disp(1:length(spreads)), 'o-')
title(['Mean Square Error with 250 Neurons or less (Max Disp)'])
xlabel(['Spread Value'])
ylabel(['Mean Square Error'])
grid on
hold off

figure (3)
hold on
plot(spreads, mse_stress(1:length(spreads)), 'o-')
title(['Mean Square Error with 250 Neurons or less (KS Stress)'])
xlabel(['Spread Value'])
ylabel(['Mean Square Error'])
grid on
hold off

figure(4)
hold on
plot(spreads, mse_CL(1:length(spreads)), 'o-')
title(['Mean Square Error with 250 Neurons or less(CL)'])
xlabel(['Spread Value'])
ylabel(['Mean Square Error'])
grid on
hold off

figure(5)
hold on
plot(spreads, mse_CD(1:length(spreads)), 'o-')
title(['Mean Square Error with 250 Neurons or less(CD)'])
xlabel(['Spread Value'])
ylabel(['Mean Square Error'])
grid on
hold off

figure(6)
hold on
plot(spreads, mse_CL_approx(1:length(spreads)), 'o-')
title(['Mean Square Error with 250 Neurons or less(CL Approx.)'])
xlabel(['Spread Value'])
ylabel(['Mean Square Error'])
grid on
hold off

%% Generate Surrogate Models 
[mse_mass, index] = min(mse_mass);
spread_mass = spreads(index);
net_mass = newrb(inputs,total_mass,goal_mse_mass, spread_mass,Q,1);
y_sim_mass(~index,:) = []; %delete all others

[mse_disp, index] = min(mse_disp);
spread_disp = spreads(index);
net_disp = newrb(inputs,max_disp,goal_mse_disp,spread_disp,Q,1);
y_sim_disp(~index,:) = []; %delete all others

[mse_stress, index] = min(mse_stress);
spread_stress = spreads(index);
net_stress = newrb(inputs,KS_stress,goal_mse_stress,spread_stress,Q,1);
y_sim_stress(~index,:) = []; %delete all others

[mse_CL, index] = min(mse_CL);
spread_CL = spreads(index);
net_CL = newrb(inputs,CL(2,:),goal_mse_CL, spread_CL,Q,1);
y_sim_CL(~index,:) = []; %delete all others

[mse_CD, index] = min(mse_CD);
spread_CD = spreads(index);
net_CD = newrb(inputs,CD(2,:),goal_mse_CD, spread_CD,Q,1);
y_sim_CD(~index,:) = []; %delete all others

[mse_CL_approx, index] = min(mse_CL_approx);
spread_CL_approx = spreads(index);
net_CL_approx = newrb(inputs,CL_approx,goal_mse_CL_approx, spread_CL_approx,Q,1);
y_sim_CL_approx(~index,:) = []; %delete all others

surr_mass       = @(x) sim(net_mass,x);
surr_disp       = @(x) sim(net_disp,x); %1/10 of span
surr_stress     = @(x) sim(net_stress,x); %switch to 1 from 100 
surr_CD         = @(x) sim(net_CD, x);
surr_CL         = @(x) sim(net_CL,x);
surr_CL_approx  = @(x) sim(net_CL_approx,x); 

%% Optimize
LW_ratio = 1.3;
%f = @(x) -surr_CL(x)/surr_CD(x);
f = @(x) surr_mass(x);
g = @(x) [  (surr_disp(x)- x(1)/10);
            (surr_stress(x) - 1);
            -(surr_CL_approx(x).*((1/2).*MainFlightCondition.rho.*(MainFlightCondition.velocity^2).*((1/2).*(x(2)+x(3)).*x(1))) - LW_ratio.*((1/0.15).*surr_mass(x)+ surr_mass(x)))];
h = @(x) [];
rng('shuffle'); %seeds random number generator based on time
population_init = 50;
population_max = 200;
population_top= (1/5)*population_init;

x_guess = (ub-lb).*rand(length(lb),1) +lb; %Random guess. literally.

% optimization 
f_min = []; sensitivity = []; x_min = []; iterations = [];
tic
parfor k = 1:10
    x_guess = (ub-lb).*rand(length(lb),1) +lb; %Random guess. literally.
    [f_min(k), sensitivity(:,k), x_min(:,k), iterations(k)] = ALM_GA_Func(f, g, h, lb, ub, x_guess, population_init, population_max, population_top);
end
toc
timelapsed_ga = toc;
disp('OPTIMIZATION COMPLETE')

% Verification
[f_min_result, index] = min(f_min);
input_opt = x_min(:,index) ;

%% Model optimum design using PATRAN/NASTRAN
[total_mass_opt, max_disp_opt, KS_stress_opt, CD_opt, CL_opt, CL_approx_opt, CD_approx_opt] = WingModeler(MainFlightCondition, airfoil_path, input_opt, 0);

%% Percent error calculations at optimum
total_mass_s = surr_mass(input_opt);
max_disp_s = surr_disp(input_opt);
KS_stress_s = surr_stress(input_opt);
CL_s = surr_CL(input_opt);
CD_s = surr_CD(input_opt);
CL_approx_s = surr_CL_approx(input_opt);

accuracy_mass   = total_mass_opt - total_mass_s;
accuracy_disp   = max_disp_opt - max_disp_s;
accuracy_stress = KS_stress_opt - KS_stress_s;
accuracy_CL = CL_opt(2,:) - CL_s;
accuracy_CD = CD_opt(2,:) - CD_s;
accuracy_CL_approx = CL_approx_opt - CL_approx_s;

pe_mass = accuracy_mass/total_mass_opt;
pe_disp = accuracy_disp/max_disp_opt;
pe_stress = accuracy_stress/KS_stress_opt;
pe_CL = accuracy_CL/CL_opt(2,:);
pe_CD = accuracy_CD/CD_opt(2,:);
pe_CL_approx = accuracy_CL_approx/CL_approx_opt;


%% for ease of copy paste
thematrix = [total_mass_s max_disp_s KS_stress_s CL_s CD_s CL_approx_s;
    total_mass_opt max_disp_opt KS_stress_opt CL_opt(2,:) CD_opt(2,:) CL_approx_opt;
    accuracy_mass accuracy_disp accuracy_stress accuracy_CL accuracy_CD accuracy_CL_approx;
    pe_mass pe_disp pe_stress pe_CL pe_CD pe_CL_approx;]

%% plot
DataVisualizer();
