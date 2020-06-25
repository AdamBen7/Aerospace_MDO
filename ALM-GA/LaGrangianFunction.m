function LaGrangian = LaGrangianFunction(f_func, h_func, g_func, x, v, r, u)
    LG_Eq = f_func(x) + v*h_func(x) + (1/2)*r*h_func(x)^2;
    g = g_func(x);
%do similar for h_func if h_func size is >1
    LaGrangian = LG_Eq;
    for j = 1:length(g)
        if (g(j) + (u(j)/r)) >= 0
            LG_InEq = u(j)*g(j) + (1/2)*r*g(j).^2;
        else %(g(x) + (u/r)) < 0
            LG_InEq = -(1/(2*r))*u(j).^2;
        end
        LaGrangian = LaGrangian + LG_InEq;
    end
end