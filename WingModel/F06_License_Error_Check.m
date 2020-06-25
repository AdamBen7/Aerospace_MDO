% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

function [LICENSE_ERROR_FLAG] = F06_License_Error_Check(VersionName)
f06file_name = sprintf('%s-fullmodel.f06',VersionName);
fid = fopen(f06file_name,'r');
tline = strtrim(fgetl(fid));

lineCounter = 1;
A = {};

marker_ref = [];
j = 1;

LICENSE_ERROR_FLAG = false;
while ischar(tline)
    try
        tline = strtrim(fgetl(fid));
        lineCounter = lineCounter +1;
        if isempty(tline)
            continue
        else
            A{j}= tline;
            if contains(A{j}, 'UNABLE TO OBTAIN LICENSES') == true
				LICENSE_ERROR_FLAG = true;
				break
            end	
            j = j + 1;
        end
    catch
        break;%EOF
    end

end
fclose(fid);

end

