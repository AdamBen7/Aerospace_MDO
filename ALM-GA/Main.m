%% Calling ALM_GA Function
close all; clear; clc; 
f = @(x) x(1).^4 - 2*(x(1).^2).*x(2) + x(1).^2 + x(1).*x(2).^2 -2*x(1) + 4; %cost function
h = @(x) (1/2)*x(1).^2 + (1/2)*x(2).^2 - 1; %equality constraint func. %normalized cuz I read it somewhere...
g1 = @(x) 0.25*x(1).^2 + 0.75*x(2).^2 - 1; %inequality constraint func.
g2 = @(x) -x(1);
g3 = @(x) x(1) - 4;
g4 = @(x) -x(2);
g5 = @(x) x(2) - 4;
g = @(x) [g1(x) g2(x) g3(x) g4(x) g5(x)]';

rng('shuffle'); %seeds random number generator based on time
population_init = 200;
population_max = 750;
population_top= (1/5)*population_init;

global ALLMEMBERS;

lb = [-5 -5]';
ub = [5 5]';

x_guess = [3 2]';

parfor i = 1:4
[f_min(i), sensitivity(:,i), x_min(:,i)] = ALM_GA_Func(f, g, h, lb, ub, x_guess, population_init, population_max, population_top);
end

%% EggHolder Function
%ans: f(512,404.2319) = -959.6407
close all; clear; clc; 
f = @(x) -(x(2) + 47)*sin(sqrt(abs(x(1)/2+(x(2)+47)))) - x(1)*sin(sqrt(abs(x(1)-(x(2)+47)))); %cost function
h = @(x) [];%equality constraint func. %normalized cuz I read it somewhere...
g = @(x) [];

rng('shuffle'); %seeds random number generator based on time
population_init = 200;
population_max = 1000;
population_top= (1/5)*population_init;

global ALLMEMBERS;

lb = [-512 -512]';
ub = [512 512]';

x_guess = [300 200]';
N = 10;

f_min = zeros(1,N);
x_min = zeros(2,N);

tic
for k = 1:N
[f_min(k), sensitivity(:,k), x_min(:,k), iterations(k)] = ALM_GA_Func(f, g, h, lb, ub, x_guess, population_init, population_max, population_top);
end
toc

%% Bird Function - constrained
%ans: f(-3.1302468,-1.5821422) = -106.7645367
close all; clear;clc;
%f = @(x) -20*exp(-0.2*sqrt(0.5*(x(1)^2 + x(2)^2))) - exp(0.5*(cos(2*pi*x(1)) + cos(2*pi*x(2)))) + exp(1) + 20; % Ackley Function
f = @(x) sin(x(2))*exp((1-cos(x(1)))^2) + cos(x(1))*exp((1-sin(x(2)))^2) + (x(1)-x(2))^2;
h = @(x) [];
g = @(x) (x(1)+5)^2 + (x(2) + 5)^2 - 25;
lb = [-10 -6.5]';
ub = [0 0]';
x_guess = [-5 -5]';

%infeasible bound. for algorithm bug check
%lb = [-1 -1]';
%5ub = [3 3]';
%x_guess = [0,0]';

rng('shuffle'); %seeds random number generator based on time
population_init = 50;
population_max = 200;
population_top= (1/5)*population_init;


N =4;
f_min = zeros(1,N);
x_min = zeros(2,N);

tic
for k = 1:N
[f_min(k), sensitivity(:,k) x_min(:,k),iterations(k)] = ALM_GA_Func(f, g, h, lb, ub, x_guess, population_init, population_max, population_top);
end
toc

%% some function that takes a while
close all; clear; clc; 
f = @(x) x(1).^4 - 2*(x(1).^2).*x(2) + x(1).^2 + x(1).*x(2).^2 -2*x(1) + 4; %cost function
h = @(x) (1/2)*x(1).^2 + (1/2)*x(2).^2 - 1; %equality constraint func. %normalized cuz I read it somewhere...
g1 = @(x) 0.25*x(1).^2 + 0.75*x(2).^2 - 1; %inequality constraint func.
g2 = @(x) -x(1);
g3 = @(x) x(1) - 4;
g4 = @(x) -x(2);
g5 = @(x) x(2) - 4;
g = @(x) [g1(x) g2(x) g3(x) g4(x) g5(x)]';
lb = [-5 -5]';
ub = [5 5]';

rng('shuffle'); %seeds random number generator based on time
population_init = 50;
population_max = 750;
population_top= (1/5)*population_init;

N =4;
f_min = zeros(1,N);
x_min = zeros(2,N);
x_guess = [3 2]';  %initial guess

tic
for k = 1:N
[f_min(k), sensitivity(:,k) x_min(:,k),iterations(k)] = ALM_GA_Func(f, g, h, lb, ub, x_guess, population_init, population_max, population_top);
end
toc


%% testing with matlab built-in
ConstraintFunc = @constraintfunctions;
tic 
for k = 1:50
    [x_ga(:,k), f_ga(k)] = ga(f,2, [],[],[],[],lb, ub,ConstraintFunc);
end
toc



%% More plotting/testing  - Eggholder Function
% Lecture 30, Slide 42/43
N=500; %sample size parameter
x1 = linspace(-512,512,N);
x2 = linspace(-512,512,N);
[X1,X2] = meshgrid(x1,x2);
the_f =  -(X2 + 47).*sin(sqrt(abs((X1/2)+(X2+47)))) - X1.*sin(sqrt(abs(X1-(X2+47))));

hold on
figure(1);
the_surf = contour3(X1,X2,the_f,100);
%the_surf = surf(X1,X2,the_f);
title('Surface Plot of Eggholder Function')
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');
grid on;
view(127.5,20);
%shading faceted;



plot3(512,404.2319,-959.5407,'d','MarkerEdgeColor','g','MarkerFaceColor','y','MarkerSize',6);
scatter3(x_min(1,:), x_min(2,:),f_min, 20,'*b');
scatter3(x_ga(1,:), x_ga(2,:),f_ga, 50,'+r');
legend('Original Surface Contours', 'Best Known Optimum', 'ALM-GA', 'Matlab-GA')

hold off;

set(gca,'FontSize',14,'FontName','Times New Roman');



%% More Plotting/Testing Mishra's Bird Function
N=500; %sample size parameter
x1 = linspace(-10,0,N);
x2 = linspace(-6.5,0,N);
[X1,X2] = meshgrid(x1,x2);
g = (X1+5).^2 + (X2 + 5).^2 >= 25;
X1(g) = NaN;
X2(g) = NaN;

the_f = sin(X2).*exp((1-cos(X1)).^2) + cos(X1).*exp((1-sin(X2)).^2) + (X1-X2).^2;


hold on
figure(1);
the_surf = surfc(X1,X2,the_f, 'FaceAlpha',0.5);

%the_surf = surf(X1,X2,the_f);
title('Surface Plot of Bird Function')
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');
grid on;
view(-37.5,30);
shading interp;

plot3(-3.1302468, -1.5821422, -106.7645367,'d','MarkerEdgeColor','g','MarkerFaceColor','y','MarkerSize',6);
scatter3(x_min(1,:), x_min(2,:),f_min, 20,'*b');
scatter3(x_ga(1,:), x_ga(2,:),f_ga, 50,'+r');
legend('Original Surface', 'Original Contour', 'Best Known Optimum', 'ALM-GA', 'Matlab-GA')

hold off;

set(gca,'FontSize',14,'FontName','Times New Roman');

%%
function [c_ineq, c_eq] = constraintfunctions(x)
    c_ineq = [((x(1)+5)^2 + (x(2) + 5)^2 - 25)];
    c_eq = [];
end

function f = objectivefunction(x)
x1 = x(:,1);% I DON"T GET IT!!!! WHYY does separating this work??
x2 = x(:,2);
f = -(x2 + 47).*sin(sqrt(abs(x1/2+(x2+47)))) - x1.*sin(sqrt(abs(x1-(x2+47)))); %cost function

end


function f = mishrafunc(x)
X1 = x(:,1);% I DON"T GET IT!!!! WHYY does separating this work??
X2 = x(:,2);
f = sin(X2).*exp((1-cos(X1)).^2) + cos(X1).*exp((1-sin(X2)).^2) + (X1-X2).^2;
end

