%% Per-Second Subset Simulation
clear; clc;
import tulrfsd.mcmc.*
T = 101; dt = 1; maxRetries = 5;

% Import Excel Files
param_files = {
    'AA_Plot.xlsx',     % alpha
    'MFuel_Plot.xlsx',  % m_fuel
    'N1_Plot.xlsx',     % N1
    'VS_Plot.xlsx',     % V
    'PA_Plot.xlsx',     % theta
    'RA_Plot.xlsx',     % phi
    'EL_Plot.xlsx',     % delta_e
    'LRW_Plot.xlsx',    % crosswind
    'RH_Plot.xlsx'      % altitude
};
n_params = numel(param_files);

thresholds = readmatrix('EThreshold.xlsx');
thresholds = thresholds(2:end,2);

P_F = zeros(T, 1);
energy_dot = zeros(T, 1);
E_vec = zeros(T+1, 1);
rankings = cell(T, 1);

% Load All Parameter Data
param_data = cell(n_params, 1);
for i = 1:n_params
    param_data{i} = readmatrix(param_files{i});
end

% Loop Setup
for t = 0:T-1
    mu = zeros(1, n_params);
    sigma = zeros(1, n_params);
    for i = 1:n_params
        values = param_data{i}(2:end, t+2);  % skip header & flight ID
        mu(i) = mean(values, 'omitnan');
        sigma(i) = std(values, 'omitnan');
    end

    threshold = thresholds(t+1);
    success = false;
    retryCount = 0;

    while ~success && retryCount < maxRetries
        try
            model = @(x) energy_model_dynamic(x, threshold, dt);
        
            % Initialize Algorithm
            algo = tulrfsd.mcmc.SubsetSimulation({ ...
                'ModelHandle', model, ...
                'EvaluationMethod', 'vector', ...
                'SampleDistribution', {mu, sigma}, ...
                'ProposalDistribution', {zeros(1,n_params), 0.5*ones(1,n_params)}, ...
                'NumberOfSamples', 1e5, ...
                'ConditionalProbability', 0.001, ...
                'SamplingMethod', 'adaptive', ...
                'Instrumentation', 'minimal'});
        
            % Suppress Printed Output
            [P_F(t+1), cause_sample, subsets] = algo.incidentprobability();
            
            if isnan(P_F(t+1)) || P_F(t+1) == 0  % You can set this tolerance as needed
                error('Invalid failure probability');
            end

            reportgen = algo.report(subsets);
            rankings{t+1} = reportgen.ranking();
            if ~istable(rankings{t+1})
                warning("At t = %d, reportgen.ranking() is not a table. Type: %s", ...
                t, class(rankings{t+1}));
            end
            success = true;
            fprintf("t = %d | P_F = %.2e\n", t, P_F(t+1));  % optional status

        catch ME
            retryCount = retryCount + 1;
            warning('Retry #%d at t = %d s failed: %s', retryCount, t, ME.message);
        end
    end
    if ~success
        warning('Max retries reached at t = %d s — storing NaN for P_F', t);
        P_F(t+1) = NaN;
    end
end

% Create PF and Ranking Table
writecell({'Time (s)', 'Failure Probability'}, 'PF_Time.xlsx', 'Sheet', 1, 'Range', 'A1');
writematrix((0:100)', 'PF_Time.xlsx', 'Sheet', 1, 'Range', 'A2');
writematrix(P_F,        'PF_Time.xlsx', 'Sheet', 1, 'Range', 'B2');

ranking_table = zeros(9, T);  % 9 rows for x1 to x9, 101 columns for time 0 to 100

for t = 1:T
    this_rank = rankings{t};

if istable(this_rank)
    for i = 1:9
        row_idx = find(strcmp(this_rank.("Distribution Name"), sprintf('x%d', i)));
        if ~isempty(row_idx)
            ranking_table(i, t) = this_rank.("Failure Sensitivity")(row_idx);
        else
            ranking_table(i, t) = NaN;
        end
    end
else
    warning("rankings{%d} is not a table — inserting NaNs", t);
    ranking_table(:, t) = NaN;
end
end

% Write to Excel
headers = arrayfun(@(t) sprintf('T%03d', t-1), 1:T, 'UniformOutput', false);  % T000 to T100
writecell(['Dist', headers], 'Ranking_Heatmap.xlsx', 'Sheet', 1, 'Range', 'A1');

dist_labels = arrayfun(@(i) sprintf('x%d', i), 1:9, 'UniformOutput', false)';
writecell(dist_labels, 'Ranking_Heatmap.xlsx', 'Sheet', 1, 'Range', 'A2');

writematrix(ranking_table, 'Ranking_Heatmap.xlsx', 'Sheet', 1, 'Range', 'B2');

% Energy Model
function margin = energy_model_dynamic(x, threshold, dt)
    alpha   = x(:,1);
    m_fuel  = x(:,2);
    N1      = x(:,3);
    V       = x(:,4);
    theta   = x(:,5);
    phi     = x(:,6);
    delta_e = x(:,7);
    v       = x(:,8);
    h       = x(:,9);

    g = 9.81; S = 525; m_empty = 182480; m_payload = 120200/3;
    T_max = 450000*4; rho = 1.186; tau = 2; mtow = 362870;
       
    m = m_empty + m_payload + m_fuel;
    q = 0.5 .* rho .* V.^2;

    beta = asin(v ./ V);
    u = V .* cos(alpha) .* cos(beta);
    v = V .* sin(beta);
    w = V .* sin(alpha) .* cos(beta);

    C_L = 0.0813 + 0.0859 .* (alpha .* 180/pi) + 0.0066 .* (delta_e .* ...
        180/pi) + 1.4077;
    C_D = ((4.0578e-4).*(alpha*180/pi).^2 + (1.36954e-3).* ...
        (alpha*180/pi) + 0.0148 + 0.1 + ((m/S)*(3.16e-5)* ...
        (mtow^(-0.215))));

    C_Z = -cos(alpha) .* C_L + sin(alpha) .* C_D;
    C_X = -sin(alpha) .* C_L - cos(alpha) .* C_D;
    C_Y = -0.012 .* beta;

    Fx = q .* S .* C_X - m .* g .* sin(theta) + T_max .* ((N1/100) .^ tau);
    Fy = q .* S .* C_Y + m .* g .* cos(theta) .* sin(phi);
    Fz = q .* S .* C_Z + m .* g .* cos(theta) .* cos(phi);

    E_dot = sqrt((Fx .* u).^2 + (Fy .* v).^2 + (Fz .* w).^2)/2;
    E0 = 4e8 + (0.5 .* m .* (v.^2) + m .* g .* h);
    E = E0 + E_dot * dt;

    margin = E - threshold;

    if any(isnan(margin)) || any(isinf(margin))
        margin = 1e6;  % or whatever large penalty you use
        return;
    end
end