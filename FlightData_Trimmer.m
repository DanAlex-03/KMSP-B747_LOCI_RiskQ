%% Trimming ZIP Files to Approach Only Data
clear; clc; tic

% Access Data Folder
ZipFold = 'H:\Flight Data';
TrimFold = 'H:\Flight Data Trim';
ZipList = dir(fullfile(ZipFold, '*.zip'));
allZipFileNames = {ZipList.name};
nZ = numel(allZipFileNames); % Number of Flights

% Ensure Trim Folder Exists
if ~exist(TrimFold, 'dir')
    mkdir(TrimFold);
end

% Prepare Empty Arrays
flightIDs = zeros(nZ, 1);
timeDifferences = zeros(nZ, 1);

% Parallel Pool Start
parpool('local');

parfor i = 1:nZ
    % Create Unique Temp Subfolder
    tempFolder = fullfile(tempdir, sprintf('temp_%d', i));
    if ~exist(tempFolder, 'dir')
        mkdir(tempFolder);
    end
    
    % Unzip CSV
    unzip(fullfile(ZipFold, allZipFileNames{i}), tempFolder);
    csvList = dir(fullfile(tempFolder, '*.csv'));
    
    % Check If CSV Exists
    if isempty(csvList)
        warning('CSV file not found in folder: %s', tempFolder);
        rmdir(tempFolder, 's');
        continue;
    end

    % Initialize NaN Placeholders
    TP1 = NaN; TP2 = NaN; TP3 = NaN; TP4 = NaN; TP5 = NaN; 
    TP6 = NaN; TP7 = NaN;
    
    % Read CSV
    T = readtable(fullfile(tempFolder, csvList(1).name));
    A = table2array(T); Row = size(A, 1);
    
    % Parameter Column Identification
    try
        Fl = find(strcmpi(T.Properties.VariableNames, 'flap_te_pos'));
        Fp = find(strcmpi(T.Properties.VariableNames, 'fphase_acms'));
        Rh = find(strcmpi(T.Properties.VariableNames, 'hralt_m'));
        Ww = find(strcmpi(T.Properties.VariableNames, 'wow'));
        Gs = find(strcmpi(T.Properties.VariableNames, 'gs_mps'));
    catch
        warning('Missing required columns in file: %s', csvList(1).name);
        continue;
    end
    
    % Extract Time Points (TO to Descent)
    if ~isempty(Fl)
        TP1 = find(A(1:ceil(0.5 * Row), Fl) / 100 > 5, 1, 'first');
    end
    if ~isempty(TP1) && ~isnan(TP1) && ~isempty(Fl)
        TP2 = find(A(TP1:ceil(0.5 * Row), Fl) / 100 < 5, 1, 'first') + TP1 - 1;
    end
    TP1 = (TP1 + TP2) / 2;
    if ~isempty(TP2) && ~isnan(TP2) && ~isempty(Fp)
        TP3 = find(A(TP2:ceil(0.7 * Row), Fp) == 5, 1, 'first') + TP2 - 1;
    end
    if ~isempty(TP3) && ~isnan(TP3) && ~isempty(Fp)
        TP4 = find(A(floor(0.5 * Row):ceil(0.9 * Row), Fp) == 6, 1, 'first') + floor(0.5 * Row) - 1;
    end
    
    % Extract Time Points (Approach)
    if ~isempty(TP4) && ~isnan(TP4) && ~isempty(Rh) && ~isempty(Gs)
        TP5 = find(A(TP4:Row, Rh) < 304.8 & ~isnan(A(TP4:Row, Gs)), 1, 'first') + TP4 - 1;
    end
    if ~isempty(TP5) && ~isnan(TP5) && ~isempty(Rh) && ~isempty(Gs)
        TP6 = find(A(TP5:Row, Rh) < 76.2 & ~isnan(A(TP5:Row, Gs)), 1, 'first') + TP5 - 1;
    end

    % Extract Time Points (Landing to Taxi)
    if ~isempty(TP6) && ~isnan(TP6) && ~isempty(Ww)
       TP7 = find(A(TP6:Row, Ww) == 0, 1, 'first') + TP6 - 1;
    end
    if ~isempty(TP7) && ~isnan(TP7) && ~isempty(Fl)
       TP8 = find(A(TP7:Row, Fl) / 100 < 5, 1, 'first') + TP7 - 1;
    end

try
    % Validate Time Point Indices
    if TP5 >= TP7 || TP5 < 1 || TP7 > Row
        error('Invalid index range for file: %s', csvList(1).name);
    end
    
    % Extract Trimmed Data
    TrimmedData = A(TP5:TP7, :);
    TrimmedTable = array2table(TrimmedData, 'VariableNames', T.Properties.VariableNames);
    
    % Save Trimmed CSV
    trimmedCSVPath = fullfile(tempFolder, csvList(1).name);
    writetable(TrimmedTable, trimmedCSVPath);
    
    % Move and Zip the Trimmed CSV
    trimmedZipPath = fullfile(TrimFold, allZipFileNames{i});
    system(sprintf('powershell -Command "& {Compress-Archive -Path \\"%s\\" -DestinationPath \\"%s\\" -Force}"', trimmedCSVPath, trimmedZipPath));

    % Calculate the Time Difference Between TP5 and TP7
    if ~isempty(TP5) && ~isempty(TP7)
        time = A(:, 1);
        timeDiff = time(TP7) - time(TP5);
        flightIDs(i) = str2double(allZipFileNames{i}(9:14)); 
        timeDifferences(i) = timeDiff;
    end

catch ME
    warning('Skipping file due to error: %s\nFlight ID: %s', ME.message, flightIDs(i));
end
    
    % Cleanup Temp Folder
    rmdir(tempFolder, 's');
end

% Create Table for Time Differences
timeDiffTable = table(flightIDs, timeDifferences, 'VariableNames', {'ID', 'TD57'});

% Write Time Differences to Excel
writetable(timeDiffTable, fullfile(TrimFold, 'Approach_Time.xlsx'));

% Close Parallel Pool
delete(gcp('nocreate'));
toc
