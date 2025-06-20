%% Energy Model Simulation (with or w/o Threshold) Calculator
clear; clc; tic

% Access Trimmed Data Folder
TrimFold = 'H:\Flight Data Trim';
TrimList = dir(fullfile(TrimFold, '*.zip'));
allTrimFileNames = {TrimList.name};
nT = numel(allTrimFileNames); % Number of Flights

% Storage for Extracted Data
flightIDs = zeros(nT, 1);
maxSamples = 0;
sampledData = cell(nT, 1);

% Moving Filter Window
k = 101;

% Parallel Pool Start
parpool('local');

parfor i = 1:nT
    try
        % Create Temp Folder for Extraction
        tempFolder = fullfile(tempdir, sprintf('temp_trim_%d', i));
        if ~exist(tempFolder, 'dir')
            mkdir(tempFolder);
        end

        % Unzip CSV
        unzip(fullfile(TrimFold, allTrimFileNames{i}), tempFolder);
        csvList = dir(fullfile(tempFolder, '*.csv'));
        
        % Check If CSV Exists
        if isempty(csvList)
            warning('CSV file not found in folder: %s', tempFolder);
            rmdir(tempFolder, 's');
            continue;
        end

        % Read CSV
        csvFilePath = fullfile(tempFolder, csvList(1).name);
        T = readtable(csvFilePath);
        A = table2array(T); Row = height(A); Col = width(A);
        
        % Identify Parameter Column
        try
            Gs = find(strcmpi(T.Properties.VariableNames, 'gs_mps'));
            Pa = find(strcmpi(T.Properties.VariableNames, 'theta_rad'));
            Ra = find(strcmpi(T.Properties.VariableNames, 'phi_rad'));
            Aa = find(strcmpi(T.Properties.VariableNames, 'aoac_rad'));
            Bh = find(strcmpi(T.Properties.VariableNames, 'hbaro_m'));
            Rh = find(strcmpi(T.Properties.VariableNames, 'hralt_m'));
            Edl = find(strcmpi(T.Properties.VariableNames, 'elv_l_rad'));
            Edr = find(strcmpi(T.Properties.VariableNames, 'elv_r_rad'));
            N11 = find(strcmpi(T.Properties.VariableNames, 'n11_RPM'));
            N12 = find(strcmpi(T.Properties.VariableNames, 'n12_RPM'));
            N13 = find(strcmpi(T.Properties.VariableNames, 'n13_RPM'));
            N14 = find(strcmpi(T.Properties.VariableNames, 'n14_RPM'));
            FQ1 = find(strcmpi(T.Properties.VariableNames, 'fqty_1_kg'));
            FQ2 = find(strcmpi(T.Properties.VariableNames, 'fqty_2_kg'));
            FQ3 = find(strcmpi(T.Properties.VariableNames, 'fqty_3_kg'));
            FQ4 = find(strcmpi(T.Properties.VariableNames, 'fqty_4_kg'));
            Ws = find(strcmpi(T.Properties.VariableNames,'ws_mps'));
            Wd = find(strcmpi(T.Properties.VariableNames,'wdir_rad'));
            Mh = find(strcmpi(T.Properties.VariableNames,'psi_mag_rad'));
        catch
            warning('Missing required columns in file: %s', csvList(1).name);
            continue;
        end

        time = A(:, 1);

        % Thrust Coefficients
        T0 = 450000*4; tau = 2;

        % Weight Estimates 
        OEW = 182480; PW = 120200/3; W = OEW+PW; MTOW = 362870;
        g = 9.81;

        % Aerodynamic and Atmospheric Coefficients
        S = 525; rho = 1.186; CDflp = 0.1; MAC = 8.75; t0 = 288.15; 
        p0 = 101325; R = 287; CLflp = 1.4077; CYb = -0.012;

        % Preliminary Wind Components Calculation
        A(:,Col+1) = A(:,Ws).*sin(A(:,Wd)-A(:,Mh)+pi); % L/R
        A(:,Col+2) = A(:,Ws).*-cos(A(:,Wd)-A(:,Mh)+pi); % H/T
        A(:,Gs) = A(:,Gs) + A(:,Col+2); % Airspeed

        % Set Constant Parameters (for Threshold Setting)
        A(:,Aa) = A(:,Aa)+10*(pi/180;
        A(:,Pa) = A(:,Pa)+10*(pi/180); % Maintain γ
        A(:,Ra) = A(:,Ra)+45*(pi/180);
        A(:,Gs) = A(:,Gs)-30;
        A(:,Col+1) = 15.4333;

        % Preliminary Smoothing
        A(:,Pa) = movmean(A(:,Pa),k,"omitnan");
        A(:,Ra) = movmean(A(:,Ra),k,"omitnan");
        A(:,Aa) = movmean(A(:,Aa),k,"omitnan");
        A(:,Rh) = movmean(A(:,Rh),k,"omitnan");
        A(:,Bh) = movmean(A(:,Bh),k,"omitnan");

        % Sideslip Calculation
        A(:,Col+3) = asin(A(:,Col+1)./A(:,Gs));
        A(:,Col+4) = CYb.*A(:,Col+3)*(180/pi);

        % Intermediate Smoothing I
        A(:,Gs) = movmean(A(:,Gs),k,"omitnan");
        A(:,Col+3) = movmean(A(:,Col+3),k,"omitnan");
        A(:,Col+4) = movmean(A(:,Col+4),k,"omitnan");

        % Air Density, Thrust, Weight, Elevator, Dynamic Pressure
        A(:,Col+5) = (p0.*(1-0.0065.*(A(:,Bh)./t0)).^5.2561)./(R.*(t0-0.0065.*A(:,Bh)))+0.117991451200426;
        A(:,Col+6) = T0.*(((A(:,N11)+A(:,N12)+A(:,N13)+A(:,N14))/400).^...
            (tau));
        A(:,Col+7) = W+((A(:,FQ1)+A(:,FQ2)+A(:,FQ3)+A(:,FQ4))*(10/2.205));
        A(:,Col+8) = ((A(:,Edl)-A(:,Edr)+(pi/2))/2)*180/pi;
        A(:,Col+9) = 0.5.*A(:,Col+5).*(A(:,Gs).^2);

        % Intermediate Smoothing II
        A(:,Col+6) = movmean(A(:,Col+6),k,"omitnan");
        A(:,Col+7) = movmean(A(:,Col+7),k,"omitnan");
        A(:,Col+8) = movmean(A(:,Col+8),k,"omitnan");

        % Aero Coefficients, Force and Velocity Components
        A(:,Col+10) = 0.0813 + 0.0859.*(A(:,Aa).*180/pi) + 0.0066.*A(:,Col+8) + CLflp;
        A(:,Col+11) = ((4.0578e-4).*(A(:,Aa)*180/pi).^2 + (1.36954e-3).*...
            (A(:,Aa)*180/pi) + 0.0148 + CDflp + ((A(:,Col+7)/S)*...
            (3.16e-5)*(MTOW^(-0.215))));

        A(:,Col+12) = -A(:,Col+10).*cos(A(:,Aa)) + A(:,Col+11).*sin(A(:,Aa));
        A(:,Col+13) = -A(:,Col+10).*sin(A(:,Aa)) - A(:,Col+11).*cos(A(:,Aa));

        A(:,Col+14) = A(:,Col+9).*S.*A(:,Col+13) - A(:,Col+7).*g.*sin(A(:,Pa)) + A(:,Col+6);
        A(:,Col+15) = A(:,Col+9).*S.*A(:,Col+4) + A(:,Col+7).*g.*cos(A(:,Pa)).*sin(A(:,Ra));
        A(:,Col+16) = A(:,Col+9).*S.*A(:,Col+12) + A(:,Col+7).*g.*cos(A(:,Pa)).*cos(A(:,Ra));

        % Model Energy Rate and Data Total Energy
        A(:,Col+17) = -sqrt(((A(:,Col+14).*A(:,Gs).*cos(A(:,Aa)).*cos(A(:,Col+3))).^2) +            ((A(:,Col+16).*A(:,Gs).*sin(A(:,Aa)).*cos(A(:,Col+3))).^2) + ...
            ((A(:,Col+15).*A(:,Gs).*sin(A(:,Col+3))).^2));
        A(:,Col+18) = 0.5.*A(:,Col+7).*(A(:,Gs).^2) + A(:,Col+7).*...
            g.*A(:,Rh);

        % Numerical Integration of Energy Rate
        Er = A(:,Col+17)/100;
        Er = flip(Er);
        E = cumtrapz(Er);
        E = E + abs(E(Row)) + 200000000; % Estimated Int Constant

        % Modifiable Output (A(:,Aa) for e.g.)
        para = E;

        % Resample to 1 Hz (due to smoothing increasing it back)
        sampledIndices = Row:-16:1;
        sampledP = para(sampledIndices);
        
        % Store Flight ID and Sampled Data
        flightIDs(i) = str2double(allTrimFileNames{i}(9:14));
        sampledData{i} = sampledP(:)';
        
        % Update Max Sample Count
        maxSamples = max(maxSamples, length(sampledP));

        % Cleanup Temp Folder
        rmdir(tempFolder, 's');

    % Anti-Run Termination Due to Error  
    catch ME
        warning('Skipping file due to error: %s\nFlight ID: %s', ME.message, flightIDs(i));
    end
end

% Create Output Table with Rows of Flights and Cols of Time
outputMatrix = nan(nT, maxSamples + 1);
outputMatrix(:, 1) = flightIDs;

% Fill In Collected Data
for i = 1:nT
    len = length(sampledData{i});
    outputMatrix(i, 2:len+1) = sampledData{i};
end

% Define Max Column Limit (150 + 1 for Flight IDs)
maxColumns = 151; 

% Initialize Column Names
columnNames = cell(1, maxColumns + 1); % +1 for Flight ID
columnNames{1} = 'ID'; % First column is flight ID

for j = 1:maxColumns
    columnNames{j+1} = sprintf('%d', j-1);
end

% Ensure Data Matrix is Limited to maxColumns
numFlights = size(outputMatrix, 1);
outputMatrix = outputMatrix(:, 1:maxColumns + 1); % Include ID column

% Convert the Matrix to Table
outputTable = array2table(outputMatrix, 'VariableNames', columnNames);

% Write to Excel (modifiable filename)
outputFilePath = fullfile(TrimFold, 'Modified_Data_Plot.xlsx');
writetable(outputTable, outputFilePath, 'WriteVariableNames', true);

% Close Parallel Pool
delete(gcp('nocreate'));
toc
