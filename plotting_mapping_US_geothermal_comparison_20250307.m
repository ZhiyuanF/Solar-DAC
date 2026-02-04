clc
clearvars

% This script plots the US specific results into US maps
% It compares the LCCO2 between geothermal and solar results


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

%% Solar DAC global data reconstruct

% Results reconstruction by filling NaN
LCCO2_uncorrected_reconstruct = NaN(size(feasibility_mask_reshaped)); % Result: 761760 x 1
% levelized cost of CO2
LCCO2_uncorrected_column = raw_results_global(:, 9); % Extract the 9th column of size 107731 x 1,
% Assign these values back to the original array
LCCO2_uncorrected_reconstruct(feasible_indices) = LCCO2_uncorrected_column;

LCCO2_uncorrected = reshape(LCCO2_uncorrected_reconstruct, 1440, 529);

%% Additional data loading and processing

% load the geothermal raw data
data_geothemal_raw = readtable('data/Geothermal_US_baseload.csv');
lat_geo = data_geothemal_raw.project_lat;
lon_geo = data_geothemal_raw.project_long;
LCOE_geo = data_geothemal_raw.LCOE;

%% data resolution reconstruction

% Flatten the grid coordinates into a list
lat_transpose = Lat';
lon_transpose = Lon';
lat_grid = lat_transpose(:);
lon_grid = lon_transpose(:);
gridPoints = [lat_grid, lon_grid];

% Set the latitude and longitude limits for the continental US.
latlim = [24 50];   % approximate latitudinal bounds
lonlim = [-125 -66]; % approximate longitudinal bounds

% Create a mask to select only grid points within the US bounds.
us_grid_mask = (lat_grid >= latlim(1)) & (lat_grid <= latlim(2)) & ...
               (lon_grid >= lonlim(1)) & (lon_grid <= lonlim(2));

% keep both feasible and in US rough regions
new_mask = us_grid_mask.*feasibility_mask_reshaped;
us_indices = find(new_mask);

gridPoints_us = [us_indices, gridPoints(us_indices, :)];

% Create copies for the updated feasible mask and a new LCOE grid
updatedFeasMask = zeros(size(feasibility_mask_reshaped)); 
LCOE_updated = nan(size(updatedFeasMask)); % initialize with NaNs

% Loop over each geothermal data point
for i = 1:length(lat_geo)

    if mod(i, 1000) == 0
        fprintf('Processed %d geothermal data points.\n', i);
    end
    % current geothermal point coordinates
    pt = [lat_geo(i), lon_geo(i)];
    
    % Compute the Euclidean distance from this point to every grid point
    % (This is a simple approximation assuming degrees are roughly uniform.)
    dists = sqrt((gridPoints_us(:,2) - pt(1)).^2 + (gridPoints_us(:,3) - pt(2)).^2);
    
    % Find the index of the grid point with the smallest distance
    [minDist, idx] = min(dists);

    original_idx = gridPoints_us(idx,1);
    
    % If the nearest grid point is within 0.25 degrees, update the arrays.
    if minDist <= 0.25
        updatedFeasMask(original_idx) = 1;      % Mark as feasible (if not already)
        LCOE_updated(original_idx) = LCOE_geo(i); % Attach the LCOE value.
        % If more than one geo point is attached to the same grid point,
        % this simply overwrites the previous one (arbitrary choice).
    end
end

%% Global LCCO2 mapping plot

figure_size = [600 100 1000 600];

updatedFeasMask_reshape = reshape(updatedFeasMask, 1440, 529);
LCCO2_uncorrected = LCCO2_uncorrected.*updatedFeasMask_reshape;
LCCO2_uncorrected(LCCO2_uncorrected == 0) = NaN;

% Set up the map
figure(3);
set(gcf,'Position',figure_size); 

% Set the latitude and longitude limits for the continental US.
latlim = [24 50];   % approximate latitudinal bounds
lonlim = [-125 -66]; % approximate longitudinal bounds

worldmap(latlim, lonlim); % Create a world map
load coastlines; % Load coastlines for reference
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black'); % Plot coastlines


% Plot the LCCO2 results
pcolorm(Lat, Lon, LCCO2_uncorrected'); % Transpose data to match dimensions

% Define a yellow-to-red colormap
% custom_colormap = [linspace(1, 1, 256)', linspace(1, 0, 256)', linspace(0, 0, 256)']; % Yellow to red

custom_colormap = [linspace(1, 1, 256)', linspace(0, 1, 256)', linspace(0, 0, 256)'; % Red to Yellow
                   linspace(1, 0, 256)', linspace(1, 1, 256)', linspace(0, 0, 256)'; % Yellow to Green
                   linspace(0, 0, 256)', linspace(1, 0, 256)', linspace(0, 1, 256)']; % Green to Blue
colormap(custom_colormap); % Apply the custom colormap
% Add color bar and labels
colorbar; % Display color scale
clim([170 270]); % Adjust color range to data


% --- Add state boundary lines ---
% Read the state boundaries shapefile (usastatelo.shp is typically provided with MATLAB)
% Note: Provide the bounding box as a double matrix.
bbox = [lonlim(1) latlim(1); lonlim(2) latlim(2)];
states = shaperead('usastatelo.shp', 'UseGeoCoords', true, 'BoundingBox', bbox);

% Plot the state boundaries
geoshow(states, 'DisplayType', 'polygon', 'FaceColor', 'none', ...
    'EdgeColor', 'black', 'LineWidth', 1, 'LineStyle', '-');


title('LCCO2 ($/ton-CO2) of Solar-driven DAC - Benchmark');
xlabel('Longitude');
ylabel('Latitude');



%% Geothermal mapping
% Set up the map
figure(4);
set(gcf,'Position',figure_size); 

LCOE_updated_reshape = reshape(LCOE_updated, 1440, 529);

geo_DAC = LCOE_updated_reshape*(0.739+0.10/0.8*3.15)+92.11+40;

% Set the latitude and longitude limits for the continental US.
latlim = [24 50];   % approximate latitudinal bounds
lonlim = [-125 -66]; % approximate longitudinal bounds

worldmap(latlim, lonlim); % Create a world map
load coastlines; % Load coastlines for reference
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black'); % Plot coastlines


% Plot the LCCO2 results
pcolorm(Lat, Lon, geo_DAC'); % Transpose data to match dimensions

% Define a yellow-to-red colormap
% custom_colormap = [linspace(1, 1, 256)', linspace(1, 0, 256)', linspace(0, 0, 256)']; % Yellow to red

custom_colormap = [linspace(1, 1, 256)', linspace(0, 1, 256)', linspace(0, 0, 256)'; % Red to Yellow
                   linspace(1, 0, 256)', linspace(1, 1, 256)', linspace(0, 0, 256)'; % Yellow to Green
                   linspace(0, 0, 256)', linspace(1, 0, 256)', linspace(0, 1, 256)']; % Green to Blue
colormap(custom_colormap); % Apply the custom colormap
% Add color bar and labels
colorbar; % Display color scale
clim([170 270]); % Adjust color range to data


% --- Add state boundary lines ---
% Read the state boundaries shapefile (usastatelo.shp is typically provided with MATLAB)
% Note: Provide the bounding box as a double matrix.
bbox = [lonlim(1) latlim(1); lonlim(2) latlim(2)];
states = shaperead('usastatelo.shp', 'UseGeoCoords', true, 'BoundingBox', bbox);

% Plot the state boundaries
geoshow(states, 'DisplayType', 'polygon', 'FaceColor', 'none', ...
    'EdgeColor', 'black', 'LineWidth', 1, 'LineStyle', '-');


title('LCCO2 ($/ton-CO2) of Geothermal DAC - Benchmark');
xlabel('Longitude');
ylabel('Latitude');

%% LCCO2 differences mapping
% Set up the map
figure(5);
set(gcf,'Position',figure_size); 

diff_lcco2 = LCCO2_uncorrected - geo_DAC;

% Set the latitude and longitude limits for the continental US.
latlim = [24 50];   % approximate latitudinal bounds
lonlim = [-125 -66]; % approximate longitudinal bounds

worldmap(latlim, lonlim); % Create a world map
load coastlines; % Load coastlines for reference
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black'); % Plot coastlines


% Plot the LCCO2 results
pcolorm(Lat, Lon, diff_lcco2'); % Transpose data to match dimensions


% Define a yellow-to-red colormap
% custom_colormap = [linspace(1, 1, 256)', linspace(1, 0, 256)', linspace(0, 0, 256)']; % Yellow to red

%custom_colormap = [linspace(1, 1, 256)', linspace(0, 1, 256)', linspace(0, 0, 256)'; % Red to Yellow
%                   linspace(1, 0, 256)', linspace(1, 1, 256)', linspace(0, 0, 256)'; % Yellow to Green
%                   linspace(0, 0, 256)', linspace(1, 0, 256)', linspace(0, 1, 256)']; % Green to Blue
%colormap(flipud(jet)); % Apply the custom colormap

% --- Diverging colormap: Blue -> White -> Red (centered at zero) ---
n = 256;                        % total colormap resolution
half = floor(n/2);

blue  = [0.00, 0.45, 0.74];     % “cool” end (negative)
white = [1.00, 1.00, 1.00];     % neutral (near zero)
red   = [0.85, 0.33, 0.10];     % “warm” end (positive)

cmap_neg = [linspace(blue(1),  white(1), half)', ...
            linspace(blue(2),  white(2), half)', ...
            linspace(blue(3),  white(3), half)'];

cmap_pos = [linspace(white(1), red(1), n-half)', ...
            linspace(white(2), red(2), n-half)', ...
            linspace(white(3), red(3), n-half)'];

cmap = [cmap_neg; cmap_pos];
colormap(flipud(cmap));

% Add color bar and labels
colorbar; % Display color scale
clim([min(diff_lcco2(:)) -min(diff_lcco2(:))]); % Adjust color range to data


% --- Add state boundary lines ---
% Read the state boundaries shapefile (usastatelo.shp is typically provided with MATLAB)
% Note: Provide the bounding box as a double matrix.
bbox = [lonlim(1) latlim(1); lonlim(2) latlim(2)];
states = shaperead('usastatelo.shp', 'UseGeoCoords', true, 'BoundingBox', bbox);

% Plot the state boundaries
geoshow(states, 'DisplayType', 'polygon', 'FaceColor', 'none', ...
    'EdgeColor', 'black', 'LineWidth', 1, 'LineStyle', '-');

% changing the background gray to clearly separate "0" with "no-data"
setm(gca, 'FFaceColor', [0.85 0.85 0.85])

title('Benchmark LCCO2 Differences ($/ton-CO2): Solor - Geothermal');
xlabel('Longitude');
ylabel('Latitude');



