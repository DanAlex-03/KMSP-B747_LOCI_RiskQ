%% Parameter Data Plotter
clear; clc;

% Read Excel File
excelFile = 'AA_Plot.xlsx';
rawTable = readtable(excelFile);

% Extract IDs and Data
flightIDs = rawTable{:, 1};
paraData = rawTable{:, 2:end};

% Filter Out Outlier Flights at 120 s (modifiable)
upperLimit = 10e10;
lowerLimit = -10e10;

nFlights = size(paraData, 1);
validMask = true(nFlights, 1);

for i = 1:nFlights
    row = paraData(i, :);
    if ~all(isnan(row)) && numel(row) >= 121
        valueAt120s = row(31);
        if valueAt120s < lowerLimit || valueAt120s > upperLimit
            validMask(i) = false;
        end
    end
end

% Apply Outlier Filter
filteredData = paraData(validMask, :);
filteredIDs = flightIDs(validMask);

% Set Up Data
t = 0:150;
filteredData = filteredData';

% Calculate Trendlines
meanp = nanmean(filteredData, 2);
stdp = nanstd(filteredData, 0, 2);
upperBound = meanp + 3 * stdp;
lowerBound = meanp - 3 * stdp;

% Plotting
figure('Color', 'w', 'Position', [100 100 1000 600]);
hold on;
plot(t, filteredData, 'Color', [0.2 0.2 0.8 0.1], 'LineWidth', 0.5);

hMean = plot(t, meanp, 'k-', 'LineWidth', 4);
hUpper = plot(t, upperBound, 'r--', 'LineWidth', 3);
hLower = plot(t, lowerBound, 'r--', 'LineWidth', 3);

xlabel('Seconds before TP7','FontSize',18);
ylabel('Angle of Attack (°)','FontSize',18);
grid on;
set(gca, 'XDir', 'reverse');
set(gca,'fontsize',14)
xlim([0, 150]);
%ylim([-20, 20]);
xticks(0:10:150);
legend([hMean, hUpper], {'Mean', '± 3σ'}, 'Location', 'northwest');

%% Energy Model-Data Plotter
clear; clc;

% Read Excel Files
dFile = 'ED_Plot.xlsx';
mFile  = 'EM_Plot.xlsx';

% Read Tables
dTable = readtable(dFile);
mTable = readtable(mFile);

% Extract IDs and Data
dData = dTable{:, 2:end};
mData = mTable{:, 2:end};

% Filter Out Outlier Flights at 120 s (modifiable)
upperLimit = 1*(10e8);
lowerLimit = 1*(10e7);

isValid = true(size(dData, 1), 1);
for i = 1:size(dData, 1)
    row = dData(i, :);
    if size(row, 2) >= 121
        val = row(31);
        if val < lowerLimit || val > upperLimit
            isValid(i) = false;
        end
    end
end

dData = dData(isValid, :);

isValid = true(size(mData, 1), 1);
for i = 1:size(mData, 1)
    row = mData(i, :);
    if size(row, 2) >= 121
        val = row(31);
        if val < lowerLimit || val > upperLimit
            isValid(i) = false;
        end
    end
end

mData = mData(isValid, :);

% Set Up Time Vector
t = 0:150;
dData = dData';
mData = mData';

% Compute Trendlines
meanActual = nanmean(dData, 2);
stdActual = nanstd(dData, 0, 2);
ubActual = meanActual + 3 * stdActual;
lbActual = meanActual - 3 * stdActual;

meanModel = nanmean(mData, 2);
stdModel = nanstd(mData, 0, 2);
ubModel = meanModel + 3 * stdModel;
lbModel = meanModel - 3 * stdModel;

% Plotting
figure('Color', 'w', 'Position', [100 100 1000 600]);
hold on;

% Display Flights
plot(t, dData, 'Color', [0.2 0.4 1 0.1]);
plot(t, mData, 'Color', [0 0.6 0 0.1]);

% Plot Trendlines
hActMean = plot(t, meanActual, 'b-', 'LineWidth', 4);
plot(t, ubActual, 'b--', 'LineWidth', 2.4);
plot(t, lbActual, 'b--', 'LineWidth', 2.4);

hModMean = plot(t, meanModel, 'g-', 'LineWidth', 4);
plot(t, ubModel, 'g--', 'LineWidth', 2.4);
plot(t, lbModel, 'g--', 'LineWidth', 2.4);

% Formatting (with Reversal of Time Axis)
set(gca, 'XDir', 'reverse');
set(gca,'fontsize',14)
xlim([0, 150]);
ylim([0, 3*10e8]);
xticks(0:10:150);
xlabel('Seconds before TP7','fontsize',18);
ylabel('Total Energy (J)','fontsize',18);
grid on;
legend([hActMean, hModMean], ...
    {'Data Mean ± 3σ', 'Model Mean ± 3σ'}, ...
    'Location', 'northwest');
