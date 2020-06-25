function [fmin, sensitivity, x_min ] = Simple_GA_Func(f, lb, ub, strlength, population_init, population_max, population_top, x_init)
rng('shuffle'); %seeds random number generator based on time
fcount = 0;

nvars = length(lb);
%epsilon = 2*10e-7;

k=0;
population = population_init;

%create an initial member
member(1) = GAMember(f,nvars,strlength,lb,ub,x_init);
populationtable(1,:) = [1 member(1).Get_functionvalue()];

%initialize population
for mem = 2:population
    member(mem) = GAMember(f,nvars,strlength,lb,ub);
    populationtable(mem,:) = [mem member(mem).Get_functionvalue()];
end

populationtable = sortrows(populationtable, 2);
while population < population_max
    k = k + 1;

    crossover_indexset = randperm(population_top); %top 20 
    for i = 1:2:population_top
        index_1 = populationtable(crossover_indexset(i),1);
        index_2 = populationtable(crossover_indexset(i+1),1);
        [new_gene1, new_gene2] = CrossOver(member(index_1),member(index_2));
        %creates new members
        member(population+i) = member(index_1).UpdateMember(new_gene1); 
        member(population+(i+1)) = member(index_2).UpdateMember(new_gene2);
        populationtable(population+i,:) = [(population+i) member((population+i)).Get_functionvalue()]; 
        populationtable(population+(i+1),:) = [(population+(i+1)) member((population+(i+1))).Get_functionvalue()]; 

    end
    population = length(member);
    populationtable = sortrows(populationtable, 2);

    %mutation
    for i = 1:ceil(0.05*population)
        mem = ceil(population*rand());
        if mem == populationtable(1,1)
            continue
        end
        new_gene = Mutate(member(mem));
        member(mem) = member(mem).UpdateMember(new_gene);
        index_mem = find(populationtable(:,1)== mem);
        populationtable(index_mem,2) = member(mem).Get_functionvalue();%update table value 
    end
    
    populationtable = sortrows(populationtable, 2);
    
end

    global ALLMEMBERS;
    ALLMEMBERS = member;
    
    fmin = member(populationtable(1,1)).Get_functionvalue();
    x_min = member(populationtable(1,1)).Get_designvariables();
    sensitivity = member(populationtable(1,1)).Get_sensitivity();

end

function [gene1,gene2] =  CrossOver(parent1,parent2)
    strlength = parent1.strlength;  %either parent works
    nvars = parent1.nvars;          %either parent works
    ends = randomgenestrip(strlength,nvars); %section which will crossover
    
    gene1 = parent1.Get_gene();
    gene2 = parent2.Get_gene();
    
    for i = 1:nvars
        for j = ends(i,1):ends(i,2) %thecrossover
            genestrip1(i,j) = gene1(i,j);
            gene1(i,j) = gene2(i,j);
            gene2(i,j) = genestrip1(i,j);        
        end
    end
end
function new_gene = Mutate(obj) %very simple mutation
    strlength = obj.strlength;
    nvars = obj.nvars;
    ends = randomgenestrip(strlength,nvars); %section which will crossover
    randombinary = round(rand(nvars,strlength)); %figure out how to change seeds
    new_gene = obj.Get_gene();
    for i = 1:nvars
        for j = ends(i,1):ends(i,2)
            new_gene(i,j) = num2str(randombinary(i,j));
        end
    end
end
        

function ends = randomgenestrip(strlength,nvars)
    ends = zeros(nvars,2);
    for i = 1:nvars
        while true
            ends(i,:) = ceil(strlength*rand(1,2)); %using ceil instead of round() to prevent zeros
            ends(i,:) = sort(ends(i,:));
            if ((ends(i,2)-ends(i,1)) < strlength-1)
                break
            end
        end
    end
end
