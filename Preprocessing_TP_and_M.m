% Data Preprocessing with M & TP Collection (MATLAB)
clear; clc; tic

% Access Data Folder
ZipFold = 'H:\Flight Data';
ZipList = dir(fullfile(ZipFold, '*.zip'));
allZipFileNames = {ZipList.name};
nZ = numel(allZipFileNames); % Number of Flights

% Clear Main Temp Folder
if exist('temp', 'dir')
    rmdir('temp', 's');
end
mkdir('temp'); % Create New Main Temp Folder

% Initialize Temp Arrays For Parfor Loop
tempTP = NaN(nZ, 9);
tempHD1 = NaN(nZ, 9);
tempHD2 = NaN(nZ, 9);

% Parallel Pool Start
parpool('local');

% Helper Function to Check 9 Elements
padTo9 = @(arr) [arr, NaN(1, 9 – numel(arr))];

% Process Files in Parallel
parfor i = 1:nZ
    % Create Unique Temp Subfolder
    tempFolder = fullfile('temp', sprintf('temp_%d', i));
    mkdir(tempFolder);

    % Unzip CSV
    unzip(fullfile(ZipFold, allZipFileNames{i}), tempFolder);
    csvList = dir(fullfile(tempFolder, '*.csv'));

    % Check If CSV Exists
    if ~isempty(csvList)
        T = readtable(fullfile(tempFolder, csvList(1).name)); 
        A = table2array(T);
        Row = height(A);
        
        % Initialize NaN Placeholders
        ID = NaN; TP1 = NaN; TP2 = NaN; TP3 = NaN; TP4 = NaN; 
        TP5 = NaN; TP6 = NaN; TP7 = NaN; TP8 = NaN;
        M1a = NaN; M2a = NaN; M3a = NaN; M4a = NaN;
        M5a = NaN; M6a = NaN; M7a = NaN; M8a = NaN;
        M1b = NaN; M2b = NaN; M3b = NaN; M4b = NaN;
        M5b = NaN; M6b = NaN; M7b = NaN; M8b = NaN;

        % Extract Flight ID
        try
            flightN = str2double(allZipFileNames{i}(9:14));
            ID = flightN;
        catch
            warning('Flight ID extraction failed for %s', 
            allZipFileNames{i});
        end

        % Parameter Column Identification
        try
            Fl = find(strcmpi(T.Properties.VariableNames, 
                 'flap_te_pos'));
            Fp = find(strcmpi(T.Properties.VariableNames, 
                 'fphase_acms'));
            Rh = find(strcmpi(T.Properties.VariableNames, 
                 'hralt_m'));
            Ww = find(strcmpi(T.Properties.VariableNames, 
                 'wow'));
            Gs = find(strcmpi(T.Properties.VariableNames, 
                 'gs_mps'));
            Pa = find(strcmpi(T.Properties.VariableNames, 
                 'theta_rad'));
            Ra = find(strcmpi(T.Properties.VariableNames, 
                 'phi_rad'));
            Vs = find(strcmpi(T.Properties.VariableNames, 
                 'hdot_2_mps'));
            Ws = find(strcmpi(T.Properties.VariableNames, 
                 'ws_mps'));
            Wd = find(strcmpi(T.Properties.VariableNames, 
                 'wdir_rad'));
            Mh = find(strcmpi(T.Properties.VariableNames, 
                 'psi_mag_rad'));
            Aa = find(strcmpi(T.Properties.VariableNames, 
                 'aoac_rad'));
            Da = find(strcmpi(T.Properties.VariableNames, 
                 'drift_rad'));
        catch
            warning('Missing required columns in file: %s', 
                    csvList(1).name);
            continue;
        end
        
        % Extract Time Points (TO to Descent)
        if ~isempty(Fl)
            TP1 = find(A(1:ceil(0.5 * Row), Fl) / 100 > 5, 1, 
                  'first');
        end
        if ~isempty(TP1) && ~isnan(TP1) && ~isempty(Fl)
            TP2 = find(A(TP1:ceil(0.5 * Row), Fl) / 100 < 5, 1, 
                  'first') + TP1 – 1;
        end
        TP1 = (TP1 + TP2) / 2;
        if ~isempty(TP2) && ~isnan(TP2) && ~isempty(Fp)
            TP3 = find(A(TP2:ceil(0.7 * Row), Fp) == 5, 1, 
                  'first') + TP2 – 1;
        end
        if ~isempty(TP3) && ~isnan(TP3) && ~isempty(Fp)
            TP4 = find(A(floor(0.5 * Row):ceil(0.9 * Row), Fp) == 
                  6, 1, 'first') + floor(0.5 * Row) – 1;
        end
        
        % Extract Time Points (Approach)
        if ~isempty(TP4) && ~isnan(TP4) && ~isempty(Rh) && 
           ~isempty(Gs)
            TP5 = find(A(TP4:Row, Rh) < 151.79 & 
                  ~isnan(A(TP4:Row, Gs)), 1, 'first') + TP4 – 1;
        end
        if ~isempty(TP5) && ~isnan(TP5) && ~isempty(Rh) &&   
           ~isempty(Gs)
        TP6 = find(A(TP5:Row, Rh) < 76.2 & ~isnan(A(TP5:Row, 
              Gs)), 1, 'first') + TP5 – 1;
        end

        % Collect Histogram Data (Approach)
        if ~isempty(TP5) && ~isnan(TP5)
            M1a = A(TP5, Gs); 
            M2a = A(TP5, Pa) * (180 / pi);
            M3a = A(TP5, Ra) * (180 / pi); 
            M4a = A(TP5, Vs);
            M5a = A(TP5, Ws) * sin(A(TP5, Wd) – A(TP5, Mh) + pi);
            M6a = A(TP5, Ws) * -cos(A(TP5, Wd) – A(TP5, Mh)+pi);
            M7a = A(TP5, Aa) * (180 / pi);
            M8a = A(TP5, Da) * (180 / pi);
        end
        if ~isempty(TP6) && ~isnan(TP6)
            M1b = A(TP6, Gs); 
            M2b = A(TP6, Pa) * (180 / pi);
            M3b = A(TP6, Ra) * (180 / pi); 
            M4b = A(TP6, Vs);
            M5b = A(TP6, Ws) * sin(A(TP6, Wd) – A(TP6, Mh) + pi);
            M6b = A(TP6, Ws) * -cos(A(TP6, Wd) – A(TP6, Mh)+pi);
            M7b = A(TP6, Aa) * (180 / pi);
            M8b = A(TP6, Da) * (180 / pi);
        end

        % Extract Time Points (Landing to Taxi)
        if ~isempty(TP6) && ~isnan(TP6) && ~isempty(Ww)
            TP7 = find(A(TP6:Row, Ww) == 0, 1, 'first') + TP6-1;
        end
        if ~isempty(TP7) && ~isnan(TP7) && ~isempty(Fl)
            TP8 = find(A(TP7:Row, Fl) / 100 < 5, 1, 'first') + 
                  TP7 – 1;
        end
        
        % Assign Data to Temp Arrays
        tempTP(i, :) = padTo9([ID, TP1, TP2, TP3, TP4, TP5, TP6, 
                       TP7, TP8]);
        tempHD1(i, :) = padTo9([ID, M1a, M2a, M3a, M4a, M5a, M6a, 
                        M7a, M8a]);
        tempHD2(i, :) = padTo9([ID, M1b, M2b, M3b, M4b, M5b, M6b, 
                        M7b, M8b]);
    else
        warning('CSV file not found in folder: %s', tempFolder);
    end

    % Delete Current Temp Subfolder
    rmdir(tempFolder, 's');
end

% Convert Arrays To Tables
TP = array2table(tempTP, 'VariableNames', ["ID","TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8"]);
HD1 = array2table(tempHD1, 'VariableNames', ["ID","M1","M2","M3","M4","M5","M6","M7","M8"]);
HD2 = array2table(tempHD2, 'VariableNames', ["ID","M1","M2","M3","M4","M5","M6","M7","M8"]);

% Export Results as XLSX
writetable(TP, 'Time_Point_Data.xlsx');
writetable(HD1, 'Histogram_Data_1.xlsx');
writetable(HD2, 'Histogram_Data_2.xlsx');

% Close Parallel Pool
delete(gcp('nocreate'));
toc
