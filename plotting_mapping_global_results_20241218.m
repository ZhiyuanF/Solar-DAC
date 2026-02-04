clc
clearvars

% This script plots the global results into global maps
% And additional result analysis for better explantion


%% data and results loading

% reading results, feasibility mask, and lat/lon coordinates 
data_results_global = load('results/Solar_stand_alone_DAC_global_20241219.mat');
data_results_mask = load('results/feasibility_mask.mat');
data_results_latlon = load('results/LatLonData.mat');

raw_results_global = data_results_global.raw_results_global;
feasibility_mask = data_results_mask.feasibility_mask;
Lat = data_results_latlon.Lat; % 529*1440 size
Lon = data_results_latlon.Lon; % 529*1440 size

feasibility_mask_reshaped = reshape(feasibility_mask, [], 1);
feasible_indices = find(feasibility_mask_reshaped == 1); % Indices where mask == 1

% read the raw solar input data
data_full_solar = load('B:\Ph.D_research\DAC\DAC_optimization\DAC_optimization\MATLAB_SolarHeat\data\reduced_solar_radiation.mat');
full_solar = data_full_solar.feasible_data;

%% LCCO2 vs Solar capacity factor plot

MarkerFaceAlpha = 0.3;
% Data for scatter plot
solar_CF = raw_results_global(:, 2);
LCCO2_uncorrected_raw = raw_results_global(:, 9);

% Create scatter plot
figure(1);
scatter(solar_CF, LCCO2_uncorrected_raw, 'b', 'filled', 'MarkerFaceAlpha', MarkerFaceAlpha);
hold on;

% Find Pareto Frontier
% Combine data into one matrix for sorting
data = [solar_CF, LCCO2_uncorrected_raw];

% Sort data by the first column (solar_CF)
sorted_data = sortrows(data, 1);

% Initialize Pareto frontier
pareto_frontier = sorted_data(1, :); % Start with the first point

% Loop through sorted points to find Pareto frontier
for i = 2:size(sorted_data, 1)
    if sorted_data(i, 2) < pareto_frontier(end, 2)
        pareto_frontier = [pareto_frontier; sorted_data(i, :)];
    end
end

% Extract x and y values of the Pareto frontier
pareto_x = pareto_frontier(:, 1);
pareto_y = pareto_frontier(:, 2);

% Fit a quadratic curve (y = ax^2 + bx + c)
coefficients_1 = polyfit(pareto_x, pareto_y, 2);

% Define the x range for the fitted curve
fitted_x = linspace(0.145, 0.26, 100);

% Evaluate the quadratic curve
fitted_y = polyval(coefficients_1, fitted_x);

% Plot the quadratic curve
plot(fitted_x, fitted_y, 'r-', 'LineWidth', 2);

% Highlight the Pareto points
scatter(pareto_x, pareto_y, 'r', 'filled', 'MarkerFaceAlpha', MarkerFaceAlpha);

% Add labels and title
xlabel('Solar Capacity Factor');
ylabel('Levelized Cost of CO2 ($/ton-CO2) - LCCO2');
title('Quadratic Fit of the Lower-bound for LCCO2 Given Solar Capacity Factor');
legend('Global Results', 'Quadratic Fit', 'LCCO2 lower bound', 'Location', 'best');


% Extract the coefficients
a_1 = coefficients_1(1);
b_1 = coefficients_1(2);
c_1 = coefficients_1(3);
% Add quadratic equation to the plot
equation_text = sprintf('y = %.0fx^2 %.0fx + %.0f', a_1, b_1, c_1);
text_position_x = min(pareto_x); % Slightly offset from the minimum x
text_position_y = max(pareto_y) - 75;  % Slightly offset from the maximum y
text(text_position_x, text_position_y, equation_text, 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w');

hold off;

%% LCCO2 offset vs Daily solar capacity factor SD across the year

% Step 1: Compute daily averages
% Reshape to 3D: [points, hours_per_day, days_per_year]
full_solar_reshaped = reshape(full_solar, size(full_solar, 1), 24, 365)./(1000*3600);
% Take mean across the 24-hour dimension (dim 2)
full_solar_averages = mean(full_solar_reshaped, 2);
% Squeeze to remove singleton dimension
full_solar_averages = squeeze(full_solar_averages); % Result: 107713x365

% Step 2: Compute standard deviation of daily averages
% Compute std across the days (dim 2)
std_full_solar_averages = std(full_solar_averages, 0, 2); % Result: 107713x1

% calculate the lower-bound theoretical optimal for each point
% Apply the formula f(x) = ax^2 + bx + c
lower_bound_LCCO2 = a_1 .* (solar_CF.^2) + b_1 .* solar_CF + c_1;
LCCO2_increase_coefficient = LCCO2_uncorrected_raw./lower_bound_LCCO2 - 1;

% Fit a quadratic curve (y = ax^2 + bx + c)
coefficients_2 = polyfit(std_full_solar_averages, LCCO2_increase_coefficient, 2);

% scatter plot to show correlation
figure(2);
scatter(std_full_solar_averages, LCCO2_increase_coefficient, 'b', 'filled', 'MarkerFaceAlpha', MarkerFaceAlpha);
% Add labels and title
xlabel('Solar Capacity Factor Standard Deviation Across 365 Days');
ylabel('LCCO2 increase (%) Relative to Theoretical Lower-bound');
title('Impact of Solar Capacity Variability on LCCO2 Increase');

hold on
% Define the x range for the fitted curve
fitted_x_2 = linspace(0.01, 0.12, 100);

% Evaluate the quadratic curve
fitted_y_2 = polyval(coefficients_2, fitted_x_2);

% Plot the quadratic curve
plot(fitted_x_2, fitted_y_2, 'r-', 'LineWidth', 2);

legend('Global Results', 'Quadratic Fit', 'Location', 'best');

% Extract the coefficients
a_2 = coefficients_2(1);
b_2 = coefficients_2(2);
c_2 = coefficients_2(3);
% Add quadratic equation to the plot
equation_text_2 = sprintf('y = %.2fx^2 %.2fx + %.2f', a_2, b_2, c_2);
text_position_x_2 = min(fitted_x_2); % Slightly offset from the minimum x
text_position_y_2 = max(fitted_y_2) - 0.1;  % Slightly offset from the maximum y
text(text_position_x_2, text_position_y_2, equation_text_2, 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w');

hold off


%% Global LCCO2 mapping plot

figure_size = [600 100 1000 600];

% Results reconstruction by filling NaN
LCCO2_uncorrected_reconstruct = NaN(size(feasibility_mask_reshaped)); % Result: 761760 x 1
% levelized cost of CO2
LCCO2_uncorrected_column = raw_results_global(:, 9); % Extract the 9th column of size 107731 x 1,
% Assign these values back to the original array
LCCO2_uncorrected_reconstruct(feasible_indices) = LCCO2_uncorrected_column;

LCCO2_uncorrected = reshape(LCCO2_uncorrected_reconstruct, 1440, 529);

% Set up the map
figure(3);
set(gcf,'Position',figure_size); 
worldmap('world'); % Create a world map
load coastlines; % Load coastlines for reference
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black'); % Plot coastlines

% Plot the average solar radiation data
pcolorm(Lat, Lon, LCCO2_uncorrected'); % Transpose data to match dimensions

% Define a yellow-to-red colormap
% custom_colormap = [linspace(1, 1, 256)', linspace(1, 0, 256)', linspace(0, 0, 256)']; % Yellow to red


%custom_colormap = [linspace(1, 1, 256)', linspace(0, 1, 256)', linspace(0, 0, 256)'; % Red to Yellow
%                   linspace(1, 0, 256)', linspace(1, 1, 256)', linspace(0, 0, 256)'; % Yellow to Green
%                   linspace(0, 0, 256)', linspace(1, 0, 256)', linspace(0, 1, 256)']; % Green to Blue

%{
% --- Sequential colormap: low cost = blue, mid = cyan, high cost = off-white ---

n = 256;

color_1  = [0.00, 0.00, 1.00];    % low cost
color_2  = [0.20, 0.70, 0.45];    % mid cost
offwhite = [0.98, 0.97, 0.94];    % high cost

n1 = floor(n/2);
n2 = n - n1;

cmap1 = [ ...
    linspace(color_1(1), cyan(1), n1)', ...
    linspace(color_1(2), cyan(2), n1)', ...
    linspace(color_1(3), cyan(3), n1)' ...
];

cmap2 = [ ...
    linspace(color_2(1), offwhite(1), n2)', ...
    linspace(color_2(2), offwhite(2), n2)', ...
    linspace(color_2(3), offwhite(3), n2)' ...
];

custom_colormap = [cmap1; cmap2];
%colormap(custom_colormap);
%}
% --- Modify parula so high end fades to white instead of yellow ---

n = 256;
cmap = parula(n);

% Define how much of the top end to replace (e.g., top 20%)
frac = 0.35;
k = round(frac * n);

% Choose off-white (better than pure white)
offwhite = [0.98, 0.97, 0.94];

% Starting color for the transition (parula at cutoff)
startColor = cmap(n-k, :);

% Replace the top k colors with a smooth fade to off-white
fade = [ ...
    linspace(startColor(1), offwhite(1), k)', ...
    linspace(startColor(2), offwhite(2), k)', ...
    linspace(startColor(3), offwhite(3), k)' ...
];

cmap(n-k+1:n, :) = fade;

colormap(cmap);
%colormap(parula);
%colormap(flipud(jet));

% Add color bar and labels
colorbar; % Display color scale
clim([min(LCCO2_uncorrected(:)) max(LCCO2_uncorrected(:))]); % Adjust color range to data
title('Levelized Cost of CO2 capture ($/ton-CO2) - Benchmark');
xlabel('Longitude');
ylabel('Latitude');

%% Global LCCO2 sensitivity to ambient correction mapping plot

% Results reconstruction by filling NaN
LCCO2_corrected_reconstruct = NaN(size(feasibility_mask_reshaped)); % Result: 761760 x 1
% levelized cost of CO2
LCCO2_corrected_column = raw_results_global(:, 10); % Extract the 9th column of size 107731 x 1,
% Assign these values back to the original array
LCCO2_corrected_reconstruct(feasible_indices) = LCCO2_corrected_column;

LCCO2_corrected = reshape(LCCO2_corrected_reconstruct, 1440, 529);

% element-wise LCCO2 difference
LCCO2_ratio_correction = (LCCO2_corrected./LCCO2_uncorrected - 1)*100;

% Set up the map
figure(4);
set(gcf,'Position',figure_size); 
worldmap('world'); % Create a world map
load coastlines; % Load coastlines for reference
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black'); % Plot coastlines

% Plot the average solar radiation data
pcolorm(Lat, Lon, LCCO2_ratio_correction'); % Transpose data to match dimensions

% Define a yellow-to-red colormap
% custom_colormap = [linspace(1, 1, 256)', linspace(1, 0, 256)', linspace(0, 0, 256)']; % Yellow to red

% custom_colormap_correction = [linspace(0, 1, 256)', linspace(1, 0, 256)', linspace(0, 0, 256)']; % Green to Red
colormap(flipud(hot)); % Apply the custom colormap
% Add color bar and labels
colorbar; % Display color scale
clim([min([LCCO2_ratio_correction(:);0]) max(LCCO2_ratio_correction(:))]); % Adjust color range to data
title('LCCO2 Increase (%) Subjected to Ambient Temperature and Relative Humidity Correction');
xlabel('Longitude');
ylabel('Latitude');

