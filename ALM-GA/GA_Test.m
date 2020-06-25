close all; clear; clc;
rng('shuffle'); %seeds random number generator based on time
lb = [-2 -2];
ub = [2 2];
A = [];
b = [];
fcount = 0;
nvars = 2;
strlength = 8;


f = @(x) 10*2 + (x(1)^2 - 10*cos(2*pi*x(1))) + (x(2)^2 - 10*cos(2*pi*x(2)));

%some B.S. constraints to get started. must change later
g1 = @(x) ((x(1)^2/9) + ((x(2)^2)/4) -1); %less than or equal to zero
g2 = @(x) x(1) - 1.5;

g = @(x) [g1(x) g2(x)];
h = @(x) [];

v = 1; %initial guess
r = 1; %initial guess
u = [1 1]; %initial guess. FML

LGFunc = @(x) LaGrangianFunction(f(x), h(x), g(x), x, v, r, u)



k = 0;
population_init = 200;
population_max = 1000;
population_top= 20;

global ALLMEMBERS;


[fmin, x_min ] = Simple_GA_Func(LGFunc, lb, ub, strlength, population_init, population_max, population_top );

%%
MemberTable = zeros(population_max,3);

for i = 1:population_max
    MemberTable(i,1:2) = ALLMEMBERS(i).Get_designvariables();
    MemberTable(i,3) = ALLMEMBERS(i).Get_functionvalue();
end

figure(1);
scatter3(MemberTable(:,1),MemberTable(:,2),MemberTable(:,3))
hold on;
scatter3(x_min(1),x_min(2), fmin, 'r', 'filled');
hold off;


