%% Data to Histogram Tool (Para Area per Flight via Integration)
clear; clc;

% Read Excel File and Convert to Table
inputFile = 'AA_Plot.xlsx';  % Replace as needed
dataTable = readtable(inputFile);
flightIDs = dataTable{:, 1};
parameterData = dataTable{:, 2:end};

nFlights = size(parameterData, 1);
integratedValues = zeros(nFlights, 1);

% Integration Loop
for i = 1:nFlights
    dataRow = parameterData(i, :);
    
    % Remove NaNs Just in Case
    validData = dataRow(~isnan(dataRow));
    
    % Trapezoidal Rule with Upper Limit Setup
    if length(validData) >= 95
        validData = validData(1:95);
    end
    integratedValues(i) = trapz(validData);
end

% Write to Excel (modifiable name)
resultTable = table(flightIDs, integratedValues, ...
    'VariableNames', {'ID', 'EA'});
writetable(resultTable, 'AA_Area.xlsx');

%% Area Histogram Plot (used for Model Validation and Parameter-LOC-I Threshold Plotting)
clear; clc;

% Read Excel
T1 = readtable("EM_Area.xlsx");
T2 = readtable("ED_Area.xlsx");

% Filter Outliers
T1f = rmoutliers(T1,'grubbs');
T2f = rmoutliers(T2,'grubbs');

% Extract Parameter Data
M = table2array(T1f(:,2)); D = table2array(T2f(:,2));  

% Error Calculation
E = (abs(table2array(T1(:,2))-table2array(T2(:,2)))./table2array(T2(:,2)))*100;
Ef = rmoutliers(E,'grubbs');

% Mean and St.Dev Comparison
muM = mean(M,'omitnan') ; muD = mean(D,'omitnan'); 
sM = std(M,'omitnan') ; sD = std(D,'omitnan'); 
diffmu = ((muM-muD)/muD)*100; diffs = ((sM-sD)/sD)*100;

% Plotting
figure(1)
h1 = histogram(M);
xlabel('Area (Js)')
ylabel('Frequency')

figure(2)
h2 = histogram(D);
xlabel('Area (Js)')
ylabel('Frequency')

figure(3)
h3a = histogram(M,'BinWidth',3e9);
hold on
h3b = histogram(D,'BinWidth',3e9);
hold off
set(gca,'fontsize',14)
xlabel('Area (Js)','FontSize',18)
ylabel('Frequency','FontSize',18)
legend('Model','Data')

figure(4)
h4 = histogram(Ef);
xlabel('Error (%)')
ylabel('Frequency')

figure(5)
h5 = histogram(D,'BinWidth',3e9);
hold on
xline(32281542989.7463,'--r','LineWidth',2); % LOC-I Line
hold off
xticks(2e10:1e10:10e10);
set(gca,'fontsize',14)
xlabel('Area (Js)','FontSize',18)
ylabel('Frequency','FontSize',18)
legend('Data','LOC-I Limit','Location','northwest')
