% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

function Target_BDF = BDF_Assembler(MainFlightCondition, VersionName)

Target_BDF = sprintf('%s-fullmodel.bdf',VersionName);
Target_ID = fopen(Target_BDF,'w');


fID=fopen('./1_Analysis_Header_144_15AoA.txt','r');

i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(Target_ID,'%s\n',newtline{j});  
end  

modelbdf_filepath = sprintf('./2_MainModelDeck-%s.bdf',VersionName);

fID=fopen(modelbdf_filepath,'r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(Target_ID,'%s\n',newtline{j});  
end

fprintf(Target_ID, '%s\n', '$ Loads for Load Case : Default');
fprintf(Target_ID,'%s\n','SPCADD   2       1');
fprintf(Target_ID,'%s\n','LOAD     2      1.      1.       1');

fprintf(Target_ID, '%s\n','$ Aeroelastic Model Parameters');
fprintf(Target_ID,'%s\n','PARAM   AUNITS  1.');

aerobdf_filepath = sprintf('./3_PartialBDF-%s.txt',VersionName);

fID=fopen(aerobdf_filepath,'r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(Target_ID,'%s\n',newtline{j});  
end  

%Z-TrimVariables
fID=fopen('./4_TrimVariables.txt','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(Target_ID,'%s\n',newtline{j});  
end  

fprintf(Target_ID, 'TRIM    1       %2.3f   %4.1f  ANGLEA  %1.4f  URDD1   0.\n',MainFlightCondition.Mach, MainFlightCondition.dynamic_pressure, MainFlightCondition.AoA_rad);
fprintf(Target_ID, '        URDD2   0.      URDD3   0.      URDD4   0.      URDD5   0.\n');
fprintf(Target_ID, '        URDD6   0.\n');
fprintf(Target_ID,'$ Referenced Coordinate Frames\n');

fprintf(Target_ID, '%s\n', 'ENDDATA 687e903c');

fclose(Target_ID);
end
