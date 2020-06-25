function MyGradients = MyGradient(f_func, x, perturb)
    
    %perturbMat = (eye(length(x))*perturb) + x;
    %by using theperturb, I no longer need to define perturbs. Remove later
    %once I verify that this finite difference gradient method is reliable/robust
    theperturb = 0.01*x; % 1 percent should be small enough based on Dr. Su
    perturbMat = (eye(length(x)).*theperturb) + x;
    MyGradients = zeros(length(x),length(f_func));
    for i = 1:length (x)
        MyGradients(i) = (f_func(perturbMat(:,i))/theperturb(i)) - (f_func(x)/theperturb(i));
        
    end

    
end
