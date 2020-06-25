% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

function [total_mass, max_disp , KS_stress, CD, CL, CL_approx, CD_approx] = WingModeler(FlightCondition, airfoil_path, inputs, sampleno)
  
NasPath = 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'; 
VersionName = strcat("newwing-",num2str(sampleno),'-',erase(char(timeofday(datetime)),':'),'-',num2str(inputs(1)),num2str(inputs(2)),num2str(inputs(3)),num2str(inputs(4)),num2str(inputs(5)),num2str(inputs(6)));
mkdircmd = sprintf('mkdir %s',VersionName);
system(mkdircmd);
DBPath_Short =  sprintf('./%s', VersionName);
SessionPath_Short = sprintf('%s/Full_Session-%s.ses', DBPath_Short, VersionName);

InputDataFileGenerator(VersionName, FlightCondition, inputs); %technically a subroutine
SessionGenerator(VersionName, airfoil_path, FlightCondition, inputs);


%
%Running Patran Session
% -b makes it a background process
%ifile init_fld.pcl -skin loads flightloads
%-sfp is... i can't remember. google it.
complete = false;
err_count = 0;
modelbdf_filename = sprintf('2_MainModelDeck-%s.bdf',VersionName);
aerobdf_filename = sprintf('3_PartialBDF-%s.txt',VersionName);
db_filename = sprintf('%s/%s.db',DBPath_Short,VersionName);
while (complete == false)
    delete(modelbdf_filename);
    delete(aerobdf_filename);
    delete(db_filename);
    msg = dos(sprintf('Patran -b -ifile init_fld.pcl -skin -sfp %s',SessionPath_Short));
    if msg ==0
        fprintf('\n%s: Patran - Execution complete: %d\n', VersionName, msg)
    else
        fprintf('\n%s: Patran - Nope, it messed up: %d\n', VersionName, msg)
        err_count = err_count + 1;
    end
    pause(5); %5 seconds should be plenty
    if (exist(modelbdf_filename, 'file') ~= false && exist(aerobdf_filename, 'file') ~= false)
        complete = true;
        break
    else
        fprintf('%s: Cannot find BDF so re-calling Patran\n', VersionName)
    end
    if err_count >= 3
        error('%s: Failed to run PATRAN error-free.', VersionName)
    end
end


%movefile(Target_BDF, VersionName);

% CALL NASTRAN Sol 144

fprintf('%s: Calling Nastran\n', VersionName);
nastran_success = false;
err_count = 0;
license_err_count = 0;
time_lapsed = 0;
%it's ugly because I need to account for crappy internet and access to license servers
while (nastran_success == false)
    dos(sprintf('del %s-fullmodel.*',VersionName));
    BDF_Name = BDF_Assembler(FlightCondition, VersionName); %placing it here so it can be remade after cleaning up temp files unless you can code a delete all but bdf
    dos(sprintf('%s %s',NasPath, BDF_Name));
    %probably not good practice but it works so whatever
    tic;
    pause(30 + 10*err_count*rand());
    fprintf('%s: Attempt #%d - Called Nastran. Attempting data extraction. \n',VersionName, (err_count+1))
    while (true)
        pause(5) %nastran's slow. so this is necessary
        try %in case the file is still empty
        if F06_License_Error_Check(VersionName) == true
            license_err_count = license_err_count + 1;
            fprintf('%s: Attempt #%d - License error #%d. Re-running Nastran! \n',VersionName, (err_count+1),license_err_count);
            if license_err_count >= 10
               error('%s: Nastran fatally failed 8 times.', VersionName)
            else
                break
            end
        end
        catch
           %nothing necessary
        end
        if (exist(sprintf('%s-fullmodel.xdb',VersionName),'file') ~= 0)
            try       
                [total_mass, max_disp, KS_stress, CD, CL, Lift, Drag] = F06DataExtractor_144(VersionName, FlightCondition);
                CL_approx = Lift/(FlightCondition.dynamic_pressure * inputs(1)*((inputs(2)+inputs(3))/2));
                CD_approx = Drag/(FlightCondition.dynamic_pressure * inputs(1)*((inputs(2)*0.12 +inputs(3)*0.12)/2));%0.12 is for NACA XX12
                nastran_success = true;
                fprintf('%s: Attempt #%d - Data Extraction successful! \n',VersionName, (err_count+1))
                break;
            catch
                fprintf('%s: Attempt #%d - Waiting on F06.\n', VersionName, (err_count+1))
            end
        end
        time_lapsed = toc;
        if err_count > 5
            error('%s: Nastran failed 5 times.', VersionName)
        elseif time_lapsed > (200 + 10*err_count) %around 180 seconds per attempt
            err_count = err_count +1;
            license_err_count = 0; %might as well reset the license error count
            fprintf('%s: Attempt #%d - Rerunning Nastran.\n', VersionName, (err_count+1))
            break;
        else
            continue
        end        
    end
end


%clean up
file_set = sprintf('%s*',VersionName);
try
movefile(modelbdf_filename, VersionName);
movefile(aerobdf_filename, VersionName);
movefile(file_set, VersionName);
catch
    fprintf('%s: One or more files still open at the end.\n', VersionName)
end
end


