% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

function [totalmass, MaxDisp, KS_Stress, CD, CL, Lift, Drag] = F06DataExtractor_144(VersionName, FlightCondition)
f06file_name = sprintf('%s-fullmodel.f06',VersionName);
cleaner_f06_name = sprintf('%s-CleanerF06.txt',VersionName);
fid = fopen(f06file_name,'r');
tline = strtrim(fgetl(fid));

lineCounter = 1;
A = {};

marker_ref = [];
j = 1;

while ischar(tline)
    try
        tline = strtrim(fgetl(fid));
        lineCounter = lineCounter +1;
        if isempty(tline)
            continue
        else
            A{j}= tline;
            j = j + 1;
        end
    catch
        break;%EOF
    end

end
fclose(fid);

%for visualization/debugging purposes.
fid2 = fopen(cleaner_f06_name,'w');
for i = 1:numel(A)
        fprintf(fid2,'%s \n',A{i});
end
fclose(fid2);

k = 0;
data_loc = [];
data_start_index = [];
data_end_index = [];
table_type= {};

k=1;
data_start_index(k) = 1;
for i = 1:length(A)
   try
       if contains(A{i}(118:122),'PAGE') == true
           table_type{k} = A{data_start_index(k)+2}; 
           data_end_index(k) = i;
           k = k + 1;
           data_start_index(k) = i+1;
       end
   catch
       continue
   end
end


%populate tables
DisplacementTable = []; %{} depending on how you code it 
QuadElementTable = []; %tensile surface?
QuadElementTable_C = []; %bending
CBARTable= [];
CRODTable = [];
PBARTable = [];
PSHELLTable = [];
PRODTable = [];
VehicleAeroTable = [];

for k = 1:length(table_type) 
    condition = table_type{k};
    switch condition
        case '0                                                  MAXIMUM  APPLIED LOADS'
            i = (data_start_index(k)+5);
            totalloadline = split(A{i});
            total_loads = [];
            for j = 1:length(totalloadline)
                total_loads = [total_loads str2double(totalloadline{j})];
            end
        case 'F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T'
            i = data_end_index(k) - 1;
            pointer = split(A{i});
            SPC_Force = [str2double(pointer(1))];
            for j = 3:length(pointer)
                SPC_Force= [SPC_Force str2double(pointer(j))];
            end
                
       case 'D I S P L A C E M E N T   V E C T O R'
            for i = (data_start_index(k)+4):(data_end_index(k)-1) 
                %DisplacementTable{end+1} = split(A{i})';
                pointer = split(A{i});
                data = [];
                for j = 3:length(pointer)
                    data = [data str2double(pointer(j))];
                end
                DisplacementTable = [DisplacementTable; data];
            end
        case 'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN'
            for i = (data_start_index(k)+5):10:(data_end_index(k)-1) 
                pointer1 = split(A{i});
                pointer2 = split(A{i+1}); %for 2nd row, fiber distance is positive
                data = [];
                data2 = [];
                for j = 4:length(pointer1)
                    data = [data str2double(pointer1(j))];
                    data2 = [data2 str2double(pointer2(j-3))];
                end
                QuadElementTable = [QuadElementTable; data];
                QuadElementTable_C = [QuadElementTable_C; data2];
            end

        case 'E L E M E N T   P R O P E R T Y   S U M M A R Y     (BY ELEMENT TYPE / ID)'
            for i = (data_start_index(k)+4):(data_end_index(k)-1)
                if (strfind(A{i}(1:42),'TOTAL MASS FOR ALL SUPPORTED ELEMENT TYPES') == true)
                    totalmassline = split(A{i}(43:end));
                    totalmasses = [];
                    for j = 1:length(totalmassline)
                        totalmasses = [totalmasses str2double(totalmassline{j})];
                    end
                    break;
                else
                    continue
                end
            end
            
        case 'E L E M E N T   P R O P E R T Y   S U M M A R Y     (BY PROPERTY TYPE / ID)'
            tabletype_str = split(A{data_start_index(k)+3});
            if ismember(tabletype_str(4), 'PSHELL,')
                Prop_ID = str2double(tabletype_str(7));
                j = data_start_index(k)+5;
                pointer1 = split(A{j});
                for i = (data_start_index(k)+5):(data_end_index(k)-1)
                    pointer1 = split(A{i});
                    data = [];
                    if ismember(pointer1(4), 'PSHELL,')
                        Prop_ID = str2double(tabletype_str(7));
                        continue;
                    elseif contains(pointer1(2),'QUAD4') == true
                        data = [str2double(pointer1(1)) Prop_ID];
                        for j = 4:length(pointer1)
                            data = [data str2double(pointer1(j))];
                        end
                        PSHELLTable = [PSHELLTable; data];
                        
                    else
                        continue %goto next line
                    end
                end
            else
                continue %elseif... for other types of properties
            end
            
            
        case 'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S'
            if contains(A{data_start_index(k)+12}, 'TOTAL VEHICLE') == true
                for i =(data_start_index(k)+15):(data_end_index(k)-2)
                    pointer = split(A{i});
                    if length(pointer) ~= 8
                        continue
                    end
                    data = [];
                    for j = 3:length(pointer)
                        data = [data str2double(pointer(j))];
                    end
                    VehicleAeroTable = [VehicleAeroTable; data];
                end
                %do the 2nd table once i understand these tables
                
            else
                continue
            end
                
        otherwise %safely ignores crap
            continue
    end                   
end


% Calculating data and outputting to optimizer
[~, indexMaxDisp] = max(abs(DisplacementTable(:,3)));%change to magnitude of translation of each element in 3dof?
MaxDisp = DisplacementTable(indexMaxDisp,3); %accounts for sign
%PShell: 4 - Area
%QuadElement: 8 - SigmaVM
rho = 50; %150 from Mulani paper. 50 from Martins paper 
KS_Stress = 0;
part_of_KS = 0;
sigma_y = 276*10^6;% reminder that this is yield stress. cuz i'm an idiot for many reasons.

for i = 1:length(QuadElementTable)
    nextVMStress = max(QuadElementTable(i,8), QuadElementTable_C(i,8));
    %KS_Stress = KS_Stress + (1/rho)*log((1/PSHELLTable(i,4))*PSHELLTable(i,4)*exp(rho*nextVMStress/sigma_y)); %Kreisselmeier and Steinhauser criteria
    part_of_KS = part_of_KS + PSHELLTable(i,4)*exp(rho*nextVMStress/sigma_y); 
end
KS_Stress = (1/rho)*log((1/sum(PSHELLTable(:,4)))*(part_of_KS));

totalmass = totalmasses(end); %could be slightly more efficient
CD = [NaN NaN]';
CL = [NaN NaN]';

%0Load is gravitaitonal load in our case so just calculate it instead of
%iterating through file. Sacrificing robustness for speed
Lift = (9.81*totalmass) + (-SPC_Force(4));
%Drag = Lift*tan(FlightCondition.AoA_rad); %induced drag
Drag = (-SPC_Force(2)); %WRONG!

CD(1) = VehicleAeroTable(2,6);
CD(2) = VehicleAeroTable(14,6);
CL(1) = VehicleAeroTable(6,6);
CL(2) = VehicleAeroTable(18,6);
%total_load = sqrt(total_loads(3)^2 + total_loads(4)^2 + total_loads(5)^2);


movefile(cleaner_f06_name, VersionName);

end

