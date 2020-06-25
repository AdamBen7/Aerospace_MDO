%% Design Variable to Binary then back proof of concept
clear;clc;
nvars = 2;
lb = [-5 -5];
ub = [5 5];
strlength = 32;

%x_min = [1.26684047660484 0.628586207895676];
x_min = [3 2];

g = (x_min - lb).*(2^strlength - 1)./(ub - lb);

for dvpointer = 1:nvars
    for bitpos = 1:strlength
        if 2^(strlength-bitpos) <= g(dvpointer)
            gene_mat(dvpointer,bitpos) = 1;
            g(dvpointer) = g(dvpointer) - 2^(strlength-bitpos);
        else
            gene_mat(dvpointer,bitpos) = 0;
        end
    end
end
for i = 1:nvars
    for j = 1:strlength
        gene(i,j) = num2str(gene_mat(i,j));
    end
end

designvar = zeros(nvars,1);
for dvpointer = 1:nvars
    designvar(dvpointer) = 0;
    for bitpos = 1:strlength
        designvar(dvpointer) = designvar(dvpointer)+ (str2double(gene(dvpointer,(bitpos)))*2^((strlength)-bitpos));%cleaner??
    end
    designvar(dvpointer) = designvar(dvpointer)* (ub(dvpointer) - lb(dvpointer))/((2^strlength)-1) + lb(dvpointer);
end