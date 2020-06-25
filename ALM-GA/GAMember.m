%isn't the most beautiful but it works. Proud first prototype of a GA code
classdef GAMember
    properties (Access = public)
        nvars;
        strlength;
        lb;
        ub;
        perturb;
    end
    properties (Access = private) %only class members have access
        f; 
        gene;
        designvariables;
        functionvalue;
        sensitivity;
    end
    methods
        %constructor
        function obj = GAMember(f, nvars, strlength, lb, ub, x_init)
            obj.f = f;
            obj.nvars = nvars;
            obj.strlength = strlength;
            obj.lb = lb;
            obj.ub = ub;
            if nargin >5
                obj.designvariables = x_init;
                g = (x_init - lb).*(2^strlength - 1)./(ub - lb);
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
                obj.gene = gene;
                obj.functionvalue = obj.f(x_init);  
                                
            else
                randombinary = round(rand(nvars,strlength)); 
                %gene = strings(nvars,strlength); %dunno why but pre-allocating
                %slows my code
                for i = 1:nvars
                    for j = 1:strlength
                        gene(i,j) = num2str(randombinary(i,j));
                    end
                end
                obj.gene = gene;
                designvar = zeros(nvars,1);
                for dvpointer = 1:nvars
                    designvar(dvpointer) = 0;
                    for bitpos = 1:strlength
                        designvar(dvpointer) = designvar(dvpointer)+ (str2double(gene(dvpointer,(bitpos)))*2^((strlength)-bitpos));%cleaner??
                    end                
                    designvar(dvpointer) = designvar(dvpointer)* (ub(dvpointer) - lb(dvpointer))/((2^strlength)-1) + lb(dvpointer);
                end
                obj.designvariables=designvar;
                obj.functionvalue = obj.f(obj.designvariables);
            end
            obj.perturb = 10e-6; %link with ALM's perturb later
            obj.sensitivity = MyGradient(obj.f, obj.designvariables,obj.perturb);
        end
        
        
        function designvar = binary2decimal(obj)
            designvar = zeros(obj.nvars,1);
            for dvpointer = 1:obj.nvars
                designvar(dvpointer) = 0;
                for bitpos = 1:obj.strlength
                    designvar(dvpointer) = designvar(dvpointer)+ (str2double(obj.gene(dvpointer,(bitpos)))*2^((obj.strlength)-bitpos));%cleaner??
                end                
                designvar(dvpointer) = designvar(dvpointer)* (obj.ub(dvpointer) - obj.lb(dvpointer))/((2^obj.strlength)-1) + obj.lb(dvpointer);
            end
        end
        
        function obj = UpdateMember(obj,new_genes) %use this to clean up redundancy of code
            obj.gene = new_genes;
            obj.designvariables = binary2decimal(obj);
            obj.functionvalue = obj.f(obj.designvariables);
        end
        
        %Get Methods
        function functionvalue = Get_functionvalue(obj)
            functionvalue = obj.functionvalue;
        end
        function sensitivity = Get_sensitivity(obj)
            sensitivity = obj.sensitivity;
        end
        function gene = Get_gene(obj)
            gene = obj.gene;
        end
        function designvars = Get_designvariables(obj)
            designvars = obj.designvariables;
        end
    end
end
