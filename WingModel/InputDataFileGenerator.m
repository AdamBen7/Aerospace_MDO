% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

function InputDataFileGenerator(VersionName,FlightCondition, inputs)
DBPath_Short =  sprintf('./%s', VersionName);
filepath_dat = sprintf('./%s/%s.dat', DBPath_Short,VersionName);

fileID = fopen(filepath_dat,'w'); 

fprintf(fileID,'global real span = %0.6f\n',inputs(1));
fprintf(fileID,'global real root_chord = %0.6f\n',inputs(2));
fprintf(fileID,'global real tip_chord = %0.6f\n',inputs(3));
fprintf(fileID,'global real panel_t1 = %0.6f\n',inputs(4));
fprintf(fileID,'global real rib_t = %0.6f\n',inputs(5));
fprintf(fileID,'global real spar1_t = %0.6f\n',inputs(6));

fprintf(fileID,'global real AoA = %0.6f\n',FlightCondition.AoA_rad);
fprintf(fileID,'global real velocity = %0.6f\n',FlightCondition.velocity);
fprintf(fileID,'global real rho = %0.3f\n',FlightCondition.rho); %0.3f based on what I saw in patran
fprintf(fileID,'global real Mach = %0.6f\n',FlightCondition.Mach);
fprintf(fileID,'global real dynamic_pressure = %0.6f\n',FlightCondition.dynamic_pressure);


spar1_location = 0.25; %25% away from LE
spar2_location = 0.6; %come up with a formula later.
fprintf(fileID,'global real spar1_location = %0.6f\n',spar1_location);
fprintf(fileID,'global real spar2_location = %0.6f\n',spar2_location);

%restricted parameters that could be additional dofs for higher fidelity
spar2_t = 0.05*inputs(6);
panel_t2 = 0.5*inputs(4);
panel_t3 = 0.75*panel_t2;
fprintf(fileID,'global real spar2_t = %0.6f\n',spar2_t);
fprintf(fileID,'global real panel_t2 = %0.6f\n',panel_t2);
fprintf(fileID,'global real panel_t3 = %0.6f\n',panel_t3);

mean_chord = (inputs(2) + inputs(3))/2;
wing_area = inputs(1)*mean_chord;

%physically dependent
fprintf(fileID,'global real mean_chord = %0.6f\n',mean_chord);
fprintf(fileID,'global real wing_area = %0.6f;\n',wing_area);

fclose(fileID);
