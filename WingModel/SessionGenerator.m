% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

function SessionGenerator(VersionName, airfoil_path, FlightCondition, inputs)
DBPath_Long = strcat(pwd,"\",VersionName);
DBPath_Short =  sprintf('./%s', VersionName);
SessionPath_Short = sprintf('%s/Full_Session-%s.ses', DBPath_Short, VersionName);


span        = inputs(1);
chord_root  = inputs(2);
chord_tip   = inputs(3);
panel_t1    = inputs(4);
rib_t       = inputs(5);
spar_t      = inputs(6);

panel_t2 = panel_t1/2; % dependent var
panel_t3 = panel_t2/2; % dependent var

rho = FlightCondition.rho;

mean_chord = (chord_root+chord_tip)/2;
area = mean_chord*span;
qchord_root = chord_root/4; %quarter chord location for spar/root bc placement


fileID = fopen(SessionPath_Short,'w'); 

fprintf(fileID,'$# ================== Session File Starts =====================================\n\n');
fprintf(fileID,'uil_file_new.go( "C:\\MSC.Software\\Patran_x64\\20190/template.db", "%s\\%s.db")\n', DBPath_Long,VersionName);
fprintf(fileID,'set_current_dir( "  %s")\n', DBPath_Long);

fprintf(fileID,'!!INPUT "%s\\%s.dat"\n',DBPath_Long, VersionName); %Attaching .dat file

%Generating All Groups
fprintf(fileID,'\n\n$# ================== Generating All Groups ===================================\n\n');

fID=fopen('./SessionFiles/1_FirstSession.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 

for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end  

%insert an airfoil generator here that exports as xlsx or csv for
%higher fidelity.

airfoil_original = readmatrix(airfoil_path);

index_bot = find(isnan(airfoil_original(:,1)));

%load airfoil_root top and bottom spline points
airfoil_top = zeros(index_bot-1,2);
airfoil_bot = zeros((size(airfoil_original,1) - index_bot),2);

for i = 1:index_bot-1
    airfoil_top(i,1) = airfoil_original(i,1);
    airfoil_top(i,2) = airfoil_original(i,3);
end

j = 1;
for i = index_bot+1:size(airfoil_original,1)
    airfoil_bot(j,1) = airfoil_original(i,1);
    airfoil_bot(j,2) = airfoil_original(i,3);
    j = j+1;
end

airfoil_root_top = airfoil_top*chord_root;
airfoil_tip_top = airfoil_top*chord_tip;
airfoil_root_bot = airfoil_bot*chord_root;
airfoil_tip_bot = airfoil_bot*chord_tip;

for i = 1:size(airfoil_top,1)
    s_airfoil_root_top{i} = sprintf("[%0.6f,0,%0.6f]",airfoil_root_top(i,1),airfoil_root_top(i,2));
    s_airfoil_tip_top{i} = sprintf("[%0.6f,%0.5f,%0.6f]",airfoil_tip_top(i,1),span, airfoil_tip_top(i,2));
end

%+1 that way, I don't start from coordinate (0,0,0)
n = size(airfoil_bot,1) -1;
for i = 1:n
    s_airfoil_root_bot{i} = sprintf("[%0.6f,0,%0.6f]",airfoil_root_bot(i+1,1),airfoil_root_bot(i+1,2));
    s_airfoil_tip_bot{i} = sprintf("[%0.6f,%0.6f,%0.6f]",airfoil_tip_bot(i+1,1),span, airfoil_tip_bot(i+1,2));
end
    %fprintf(fileID,

% Top Spline Generation
nextpoint = 1;
fprintf(fileID, 'STRING asm_create_grid_xyz_created_ids[VIRTUAL]\n');
fprintf(fileID, 'asm_const_grid_xyz( "%d", "',nextpoint);
fprintf(fileID, s_airfoil_root_top{1});
fprintf(fileID, '" //@\n"');
for i = 2:length(s_airfoil_root_top)
    if mod(i,2)==1
        fprintf(fileID,'"');
    end
    fprintf(fileID,s_airfoil_root_top{i});
    if i == length(s_airfoil_root_top)
        fprintf(fileID,'%s\n','", "Coord 0",  @');
    elseif mod(i,2)==0
        fprintf(fileID,'%s\n','" // @');
    end
end
fprintf(fileID,'asm_create_grid_xyz_created_ids )\n');
fprintf(fileID,'%s\n','$? YESFORALL 1000034');

nextpoint = nextpoint + length(airfoil_root_top);

fprintf(fileID, 'asm_const_grid_xyz( "%d", "',nextpoint);
fprintf(fileID, s_airfoil_tip_top{1});
fprintf(fileID, '" //@\n"');
for i = 2:length(s_airfoil_tip_top)
    if mod(i,2)==1
        fprintf(fileID,'"');
    end
    fprintf(fileID,s_airfoil_tip_top{i});
    if i == length(s_airfoil_tip_top)
        fprintf(fileID,'%s\n','", "Coord 0",  @');
    elseif mod(i,2)==0
        fprintf(fileID,'%s\n','" // @');
    end
end
fprintf(fileID,'asm_create_grid_xyz_created_ids )\n');
fprintf(fileID,'%s\n','$? YESFORALL 1000034');

fprintf(fileID,'STRING sgm_curve_bspline_created_ids[VIRTUAL]\n');
%fprintf(fileID,'sgm_const_curve_bspline( "1", "Point 1:66", 10, TRUE, 1, FALSE,   @\n');
fprintf(fileID,'sgm_const_curve_bspline( "1", "Point 1:%d", 10, TRUE, 1, FALSE,   @\n',length(airfoil_root_top));
fprintf(fileID,'sgm_curve_bspline_created_ids )\n');
%fprintf(fileID,'sgm_const_curve_bspline( "2", "Point 67:132", 10, TRUE, 1, FALSE,   @\n');
fprintf(fileID,'sgm_const_curve_bspline( "2", "Point %d:%d", 10, TRUE, 1, FALSE,   @\n', (length(airfoil_tip_top)+1), length(airfoil_tip_top)*2);
fprintf(fileID,'sgm_curve_bspline_created_ids )\n');


%Bottom Spline Generation
fprintf(fileID,'ga_group_current_set( "Bot_Surf" )\n');
nextpoint = nextpoint + length(s_airfoil_tip_top);
fprintf(fileID, 'STRING asm_create_grid_xyz_created_ids[VIRTUAL]\n');
fprintf(fileID, 'asm_const_grid_xyz( "%d", "',nextpoint);
fprintf(fileID, s_airfoil_root_bot{1});
fprintf(fileID, '" //@\n"'); 
    
for i = 2:length(s_airfoil_root_bot)
    if mod(i,2)==1
        fprintf(fileID,'"');
    end
    fprintf(fileID,s_airfoil_root_bot{i});
    if i == length(s_airfoil_root_bot)
        fprintf(fileID,'%s\n','", "Coord 0",  @');
    elseif mod(i,2)==0
        fprintf(fileID,'%s\n','" // @');
    end
end
fprintf(fileID,'asm_create_grid_xyz_created_ids )\n');
fprintf(fileID,'%s\n','$? YESFORALL 1000034');

nextpoint = nextpoint + length(s_airfoil_root_bot);

fprintf(fileID, 'asm_const_grid_xyz( "%d", "',nextpoint);
fprintf(fileID, s_airfoil_tip_bot{1});
fprintf(fileID, '" //@\n"');
    
for i = 2:length(s_airfoil_tip_bot)
    if mod(i,2)==1
        fprintf(fileID,'"');
    end
    fprintf(fileID,s_airfoil_tip_bot{i});
    if i == length(s_airfoil_tip_bot)
        fprintf(fileID,'%s\n','", "Coord 0",  @');
    elseif mod(i,2)==0
        fprintf(fileID,'%s\n','" // @');
    end
end
fprintf(fileID,'asm_create_grid_xyz_created_ids )\n');
fprintf(fileID,'%s\n','$? YESFORALL 1000034');

fprintf(fileID,'STRING sgm_curve_bspline_created_ids[VIRTUAL]\n');
%fprintf(fileID,'sgm_const_curve_bspline( "3", "Point 1 133:197 66", 10, TRUE, 1, FALSE,   @\n');
fprintf(fileID,'sgm_const_curve_bspline( "3", "Point 1 %d:%d %d", 10, TRUE, 1, FALSE,   @\n',(2*length(airfoil_root_bot)+1),(nextpoint-1),(length(airfoil_root_bot)));
fprintf(fileID,'sgm_curve_bspline_created_ids )\n');
%fprintf(fileID,'sgm_const_curve_bspline( "4", "Point 67 198:262 132", 10, TRUE, 1, FALSE,   @\n');
fprintf(fileID,'sgm_const_curve_bspline( "4", "Point %d %d:%d %d", 10, TRUE, 1, FALSE,   @\n',(length(airfoil_tip_bot)+1),(nextpoint),(length(airfoil_tip_bot)*4 - 2),(length(airfoil_root_bot)*2));
fprintf(fileID,'sgm_curve_bspline_created_ids )\n');

fprintf(fileID,'STRING asm_delete_any_deleted_ids[VIRTUAL]\n');
fprintf(fileID,'asm_delete_point( "Point 1:262", asm_delete_any_deleted_ids )',((length(airfoil_root_bot)*2)));
 

%Generating All Groups
fprintf(fileID,'\n\n$# ================== Generating 3D Geometry ============================\n\n');

fID=fopen('./SessionFiles/3_Generate_Geo_and_Struc.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end  

fID=fopen('./SessionFiles/4_Generate_Mesh.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end  

%might be better practice to define all of this inside SES using dat file. 
fprintf(fileID, '%s\n','$# =================== MPC Nodal Placement ===============');
fprintf(fileID, '%s\n','STRING asm_create_grid_xyz_created_ids[VIRTUAL]');
fprintf(fileID, 'asm_const_grid_xyz( "90001", "[%0.6f 0 0]", "Coord 0",  @\n',qchord_root);
fprintf(fileID, '%s\n', 'asm_create_grid_xyz_created_ids )');
%no need for tip mpc for now
fprintf(fileID, '%s\n', 'STRING fem_create_nodes__nodes_created[VIRTUAL]');
fprintf(fileID, '%s\n', 'fem_create_nodes_1( "Coord 0", "Coord 0", 3, "90001", " Point 90001",  @');
fprintf(fileID, '%s\n', 'fem_create_nodes__nodes_created )');

fID=fopen('./SessionFiles/5_Generate_MPC_and_BC.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end  

fID=fopen('./SessionFiles/6_Aeroelastic.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end

fprintf(fileID, '%s\n', 'fields_create_general( "SUPER_GROUP_AeroSG2D", 2, 5, 2, "Real", "Coord 0", "", @');
fprintf(fileID, '%s\n', '  0, 0, 0, 0 )');
fprintf(fileID, '%s\n','fields_create_general_term( "SUPER_GROUP_AeroSG2D", 0, 0, 0, 128,   @');
%5.993894.*0.754225*4.520746*Coord 0*3*1.226*Full Model" )
fprintf(fileID, '"%0.6f.*%0.6f*%0.6f*Coord 0*3*%0.3f*Half Model" )\n\n',span,mean_chord,area,rho); %set half model or full model here
fprintf(fileID, '%s\n', 'flds_spline_create( "TPSpline_1", "Thin Plate", ["General", "Rigid Attach"], "", [ @');
fprintf(fileID, '%s\n', '"default_group"], "", ["FlatPlate"] )');

%useless 4 lines. Cuz can't do this via pcl for some odd reason.
fprintf(fileID, '%s\n', 'fields_create_general( "STRUCT_MODEL_-1", 2, 5, 2, "Real", "Coord 0", "", 0,  @');
fprintf(fileID, '%s\n', '0, 0, 0 )');
fprintf(fileID, '%s\n', 'fields_create_general_term( "STRUCT_MODEL_-1", 0, 0, 0, 250,  @');
fprintf(fileID, '%s\n', '"TRUE**FALSE*Lumped*1.*TRUE**FALSE*TRUE*20.*TRUE*1.*TRUE" )');


fID=fopen('./SessionFiles/7_Subcase_Create_1.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end  

fprintf(fileID,'mscnastran_subcase.create_real_param( "AERO MACH NUMBER", %0.8f )\n',FlightCondition.Mach);
fprintf(fileID,'mscnastran_subcase.create_real_param( "AERO DYNAMIC PRESSURE", %4.0f. )\n',FlightCondition.dynamic_pressure);

fID=fopen('./SessionFiles/7_Subcase_Create_2.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end  


fprintf(fileID, '$# ========================= Export BDF Fragments =============================\n');

fprintf(fileID, 'flds_partial_bdf_export( "Aero Groups", ["FlatPlate"], [ @\n');
fprintf(fileID, '"Aerodynamic Physical Data", "No Rigid Body Motions", "No General Controllers" @\n');
fprintf(fileID, ', "No Aero Coordinate Frames"], ["C:\\Users\\localbenabbou\\Box Sync\\Research - MDO\\WingOptFramework\\WingModel\\3_PartialBDF-%s.txt", "Overwrite"] )\n',VersionName);

fprintf(fileID, 'flds_partial_bdf_export( "Splines", ["TPSpline_1"], ["Spline ID Based", ""], [ @\n');
fprintf(fileID, '"C:\\Users\\localbenabbou\\Box Sync\\Research - MDO\\WingOptFramework\\WingModel\\3_PartialBDF-%s.txt", "Append"] )\n',VersionName);

fprintf(fileID, 'uil_pref_analysis.set_analysis_preference( "MSC.Nastran", "Structural", ".bdf" @\n');
fprintf(fileID, ', ".op2", "No Mapping" )\n');
fprintf(fileID, 'jobfile.open( "2_MainModelDeck-%s", "ANALYZE NO JOBFILE" )\n',VersionName);
fprintf(fileID, 'msc_delete_old_files( "2_MainModelDeck-%s", ".bdf", ".op2" )\n',VersionName);
fprintf(fileID, 'jobfile.write_spl( "/* Jobfile for PATNAS created %%A%% at %%A%% */", ["18-May-20" @\n');
fprintf(fileID, ', "14:06:33"] )\n');
fprintf(fileID, 'jobfile.writec( "", "TRANSLATOR = pat3nas" )\n');


fprintf(fileID, 'jobfile.writec( "DATABASE", "%s\\%s.db" )\n',DBPath_Long,VersionName);
fprintf(fileID, 'jobfile.writec( "JOBNAME", "2_MainModelDeck-%s" )\n',VersionName );


fID=fopen('./SessionFiles/8_BDF_Export_P2.ses','r');
i=0;  
while ~feof(fID)
    tline = fgetl(fID);
    i = i+1;
    newtline{i} = tline;
end
fclose(fID); 
for j=1:1:i  
    fprintf(fileID,'%s\n',newtline{j});  
end  

fprintf(fileID, 'mscnastran_job.associate_subcases( "101", "2_MainModelDeck-%s", 1, ["Default"] )\n',VersionName );
fprintf(fileID, 'analysis_submit_2( "MSC.Nastran", "2_MainModelDeck-%s" )\n',VersionName ); 


fclose(fileID);
end