% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------
function [f_min, sensitivity, x_min, k] = ALM_GA_Func(f, g, h, lb, ub, x_guess, population_init, population_max, population_top )
rng('shuffle'); %seeds random number generator based on time

%x_guess can be randomly generated instead?
k = 1;
x(:,k) = x_guess;

%In case constraints are non-existent
if isempty(h(x(:,k)))
    h = @(x) 0;
end
if isempty(g(x(:,k)))
    g = @(x) 0;
end

%optimization parameters
K = inf;
alfa = 2; beta = 2;  %initialized guesses
epsilon1 = 10e-6; epsilon2 = 10e-6;
strlength = 16;
perturb = 10e-5;

v(:,k) = ones(1,length(h(x(:,k)))); %initial guess
r(:,k) = ones(1,1); %initial guess
u(:,k) = ones(1,length(g(x(:,k)))); %initial guess. FML

max_iters = 100; %100 is more than enough. Just give up past that
while k < max_iters
    LGFunc = @(the_x) LaGrangianFunction(f, h, g, the_x, v(:,k), r(k), u(:,k));
    
    %Step 2
    k = k + 1;
    
    %Step 3 Calls GA Function
    [f_min, sensitivity, x_min] = Simple_GA_Func(LGFunc, lb, ub, strlength, population_init, population_max, population_top, x(:,k-1));
    x(:,k) = x_min;   
    f_ALM(k) = f(x(:,k));
    
    %Step 4
    ConstraintVioParams(:,k) = [abs(h(x(:,k))); abs(max([g(x(:,k-1)); -u(:,k-1)/r(k-1)]))];   
    K_bar(k) = max(ConstraintVioParams(:,k));
    
    %incase some don't change
    v(:,k) = v(:,k-1);
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
    v(:,k) = v(:,k-1) + r(k-1)*h(x(:,k));
    u(:,k) = u(:,k-1) + r(k-1)*(max(g(x(:,k)), -u(:,k-1)/r(k-1))); %hopefully it's correct
    if K_bar(k) <= K/alfa
        K = K_bar(k);
        continue;
    end
    %Step 7
    r(k) = beta*r(k-1);
    K = K_bar(k);
end
    if k >= max_iters %experimental portion
        %set outputs as NaN if not converged
        f_min = f_min*NaN;
        x_min = x_min*NaN;
        sensitivity = sensitivity*NaN;
    end 
        
end