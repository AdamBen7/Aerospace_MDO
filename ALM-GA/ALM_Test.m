%% Using Augmented Lagrange Multiplier Method with GA
close all; clear; clc; 
f = @(x) x(1).^4 - 2*(x(1).^2).*x(2) + x(1).^2 + x(1).*x(2).^2 -2*x(1) + 4; %cost function
h = @(x) (1/2)*x(1).^2 + (1/2)*x(2).^2 - 1; %equality constraint func. %normalized cuz I read it somewhere...
g1 = @(x) 0.25*x(1).^2 + 0.75*x(2).^2 - 1; %inequality constraint func.
g2 = @(x) -x(1);
g3 = @(x) x(1) - 4;
g4 = @(x) -x(2);
g5 = @(x) x(2) - 4;
g = @(x) [g1(x) g2(x) g3(x) g4(x) g5(x)]';
perturb = 10e-4;

rng('shuffle'); %seeds random number generator based on time
population_init = 200;
population_max = 750;
population_top= (1/5)*population_init;

global ALLMEMBERS;

lb = [-5 -5]';
ub = [5 5]';
strlength = 16;

k=1;
x(:,k) = [3 2]';  %initial guess

v(k) = 1; %initial guess
r(k) = 1; %initial guess
u(:,k) = ones(length(g(x(:,k))),1); %initial guess. FML

K = inf;
alfa = 2; beta = 2;  %initialized guesses
epsilon1 = 10e-6; epsilon2 = 10e-6;

tic
%LaGrangianFunction(f, h, g, x(k,:), v(k), r(k), u(k,:))
while true
    LGFunc = @(the_x) LaGrangianFunction(f, h, g, the_x, v(k), r(k), u(:,k));

    %Step 2
    k = k + 1;
    %Step 3
    [fmin, sentivity, x_min] = Simple_GA_Func(LGFunc, lb, ub, strlength, population_init, population_max, population_top, x(:,k-1));
    x(:,k) = x_min;   
%    lb = [(x(k,1)-5/k),(x(k,2)-5/k) ];
%    ub = [(x(k,1)+5/k),(x(k,2)+5/k) ];

    f_ALM(k) = f(x(:,k));
    %Step 4
    ConstraintVioParams(:,k) = [abs(h(x(:,k))); abs(max([g(x(:,k-1)); -u(:,k-1)/r(k-1)]))];   
    K_bar(k) = max(ConstraintVioParams(:,k));
    
    %incase some don't change
    v(k) = v(k-1);
    u(:,k) = u(:,k-1);
    r(k) = r(k-1);
    
    if (K_bar(k) <= epsilon1)|| (norm(MyGradient(LGFunc,x(:,k),perturb)) <= epsilon2*max(1,norm(x(:,k))))
        break; 
    end
    %Step 5
    if K_bar(k) >= K 
        r(k) = beta*r(k-1);
        continue;
    end
    %Step 6
    v(k) = v(k-1) + r(k-1)*h(x(:,k));
    u(:,k) = u(:,k-1) + r(k-1)*(max(g(x(:,k)), -u(:,k-1)/r(k-1))); %hopefully it's correct
    
    if K_bar(k) <= K/alfa
        K = K_bar(k);
        continue;
    end
    %Step 7
    r(k) = beta*r(k-1);
    K = K_bar(k);
end

x_ALM = x;

%useful for checking constraints are getting minimized
Max_Constraint = max(ConstraintVioParams, [],2);

toc

%% Using Augmented Lagrange Multiplier Method
clear f h g1 g2 g3 g4 g5 v r u 
f = @(x) x(1).^4 - 2*(x(1).^2).*x(2) + x(1).^2 + x(1).*x(2).^2 -2*x(1) + 4; %cost function
h = @(x) (1/2)*x(1).^2 + (1/2)*x(2).^2 - 1; %equality constraint func. %normalized cuz I read it somewhere...
g1 = @(x) 0.25*x(1).^2 + 0.75*x(2).^2 - 1; %inequality constraint func.
g2 = @(x) -x(1);
g3 = @(x) x(1) - 4;
g4 = @(x) -x(2);
g5 = @(x) x(2) - 4;
g = @(x) [g1(x) g2(x) g3(x) g4(x) g5(x)];

k=1;
x2(k,:) = [3 2];  %initial guess

v(k) = 1; %initial guess
r(k) = 1; %initial guess
u(k,:) = ones(1,length(g(x2(k,:)))); %initial guess. FML

K = inf;
alfa = 2; beta = 2;  %initialized guesses
epsilon1 = 10e-6; epsilon2 = 10e-6;

%LaGrangianFunction(f, h, g, x(k,:), v(k), r(k), u(k,:))
while true
    LGFunc2 = @(the_x) LaGrangianFunction(f, h, g, the_x, v(k), r(k), u(k,:));
    
    %Step 2
    k = k + 1;
    %Step 3
    x2(k,:) = fminunc(LGFunc2, x2(k-1,:));   
    
    
    f_ALM2(k) = f(x2(k,:));
    %Step 4
    ConstraintVioParams2(k,:) = [abs(h(x2(k,:))) abs(max(g(x2(k,:)), -u(k-1,:)/r(k-1)))];   
    K_bar(k) = max(ConstraintVioParams2(k,:));
    
    %incase some don't change
    v(k) = v(k-1);
    u(k,:) = u(k-1,:);
    r(k) = r(k-1);
    
    if (K_bar(k) <= epsilon1)|| (norm(gradient(x2(k,:))) <= epsilon2*max(1,norm(x2(k,:))))
        break;
         
    end
    %Step 5
    if K_bar(k) >= K 
        r(k) = beta*r(k-1);
        continue;
    end
    %Step 6
    v(k) = v(k-1) + r(k-1)*h(x2(k,:));
    u(k,:) = u(k-1,:) + r(k-1)*(max(g(x2(k,:)), -u(k-1,:)/r(k-1))); %hopefully it's correct
    
    if K_bar(k) <= K/alfa
        K = K_bar(k);
        continue;
    end
    %Step 7
    r(k) = beta*r(k-1);
    K = K_bar(k);
end
x_ALM2 = x2;



%% Marching plot
f = @(x1,x2) x1.^4 - 2*(x1.^2).*x2 + x1.^2 + x1.*x2.^2 -2*x1 + 4; %cost function
h = @(x1,x2) (1/2)*x1.^2 + (1/2)*x2.^2 - 1; %equality constraint func. %normalized cuz I read it somewhere...
g1 = @(x1,x2) 0.25*x1.^2 + 0.75*x2.^2 - 1; %inequality constraint func.
g2 = @(x1) -x1;
g3 = @(x1) x1 - 4;
g4 = @(x2) -x2;
g5 = @(x2) x2 - 4;

[x1, x2] = meshgrid(-3:0.05:4, -3:.05:4);
func = f(x1, x2);
contour(x1,x2,func,50)
xlim([-.1 3.125]);
ylim([-0.1 2.5]);
hold on;

plot(x_ALM(:,1),x_ALM(:,2),'r*-')
plot(x_ALM2(:,1),x_ALM2(:,2),'g*-')

cv1=[0:0.1:7];
const1 = contour(x1,x2,g1(x1,x2),cv1,'y');
cv1=[0 0.001]
const1 = contour(x1,x2,g1(x1,x2),cv1,'k');
cv2=[0:0.05:7];
const2 = contour(x1,x2,g2(x1),cv2,'y');
cv2=[0 0.001]
const2 = contour(x1,x2,g2(x1),cv2,'k');
cv3=[0:0.05:7];
const3 = contour(x1,x2,g3(x1),cv3,'y');
cv3=[0 0.001]
const3 = contour(x1,x2,g3(x1),cv3,'k');
cv4=[0:0.05:7];
const4 = contour(x1,x2,g4(x2),cv4,'y');
cv4=[0 0.001]
const4 = contour(x1,x2,g4(x2),cv4,'k');
cv5=[0:.1:7];
const5 = contour(x1,x2,g5(x2),cv5,'y');
cv5=[0 0.001]
const5 = contour(x1,x2,g5(x2),cv5,'k');

cv6=[0 0.003]
const6 = contour(x1,x2,h(x1,x2),cv6,'b');

text(.25,1.25,'Feasible Line', 'Color','blue')

plot(1,1,'o', 'MarkerSize',12, 'MarkerEdgeColor','black','MarkerFaceColor','yellow');
text(1.1,1.05,'Optimal Point') 
legend ('Cost Function Contour','ALM-GA', 'ALM-Fminunc')

