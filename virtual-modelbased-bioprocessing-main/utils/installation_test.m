%% Script to test all prerequisites are installed
%   
%   The following prerequisites are needed for the library to work
%   properly:
%       - Matlab R2022a (or newer)
%       - Matlab Optimization toolbox
%       - Matlab Parallel Computing toolbox
%       - Mosek 10 solver (or newer) (RECOMMENDED otherwise should be
%       disabled, see macroscopic-modelling code)
%   
%   The script will check these dependencies are installed correctly.
%   NOTE: YOU NEED TO ADD MOSEK PATH TO MATLAB USING
%   e.g. addpath('C:\ffff\Mosek\10.2\toolbox\r2017aom')
%   THIS NEEDS TO BE DONE ONCE WHEN USING THE LIBRARY.

disp('=== Checking all prerequisites are correctly installed')

% Checking Matlab version

disp('=== Checking Matlab version...')
vr = version('-release');
if(str2num(vr(1:end-1))>=2022)
    disp('[O] Matlab version is compatible.')
else
    disp(strcat('[X] Detected Matlab release',vr,'. The library requires a release newer than R2022a to function properly.'))
end

% Checking Matlab Toolboxes 

disp('=== Checking Matlab toolboxes installation...')
v = ver(); optimization_toolbox_detected = 0; parallel_computing_toolbox_detected = 0;
for i = 1 : length(v)
    if strcmp(v(i).Name,'Optimization Toolbox')
        disp('[O] Optimization Toolbox detected.')
        optimization_toolbox_detected = 1;
    end
    if strcmp(v(i).Name,'Parallel Computing Toolbox')
        disp('[O] Parallel Computing Toolbox detected.')
        parallel_computing_toolbox_detected = 1;
    end    
end

if optimization_toolbox_detected == 0
    disp('[X] Optimization Toolbox is not installed. The library requires the toolbox to function properly.')
end
if parallel_computing_toolbox_detected == 0
    disp('[X] Parallel Computing Toolbox is not installed. The library requires the toolbox to function properly.')
end

% Checking MOSEK Installation

disp('=== Checking MOSEK installation and version...')

% Try calling mosekdiag and capturing the output
try
    info = evalc('mosekdiag'); % Capture command output as a string
    % Extract the version number using a regular expression
    split_info = splitlines(info);
    mosek_version_whole_string = split_info{4};
    mosek_version_substring = mosek_version_whole_string(strfind(mosek_version_whole_string,':')+2:end);
    parts = split(mosek_version_substring, '.'); % Split by the dot
    versions = str2double(parts); % Convert to numeric array

    if ~isempty(versionMatch)
        majorVersion = versions(1);
        minorVersion = versions(2);
        fprintf('=== MOSEK detected: Version %d.%d\n', majorVersion, minorVersion);
        
        if majorVersion >= 10
            disp('[O] MOSEK version is suitable (>=10).');
        else
            disp('[X] MOSEK version is too old (<10).');
        end
    else
        disp('[X] Could not determine MOSEK version.');
    end

catch
    disp('[X] MOSEK is not installed or not accessible. The library might not work as intended.');
end