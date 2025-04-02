clc
clearvars

%% Common inputs

% maximum profit limit for plot
profit_max = 50000;

% maximum CO2 limit for plot
CO2_max = 500;

% initializing figure index
figure_index = 1;


%% Section 1: Optimal storage capacity and temperature plot

% for particular slice of cp_s and cr (solar heating collector capactiy and 
% concentation ratio), what is the optimal storage size (investment cost)
% and the optimal operation temperature

slice_data_raw = readtable('TX_results_summary_stand_alone_sand_solar.csv');

% x-value is the cp: capacity of solar heating collectors
x1 = slice_data_raw{:,"h_hat"};
% y-value is the cr: concentration ratio
y1 = slice_data_raw{:,"T_design"};

% getting the profit per unit investment
slice_data_raw.unit_profit = slice_data_raw{:,"Total_profit"}./slice_data_raw{:,"CAPEX_total"};
% getting the CO2 per unit investment
slice_data_raw.unit_CO2 = slice_data_raw{:,"Total_CO2"}./slice_data_raw{:,"CAPEX_total"};

% z-value is the unit profit
z_profit = slice_data_raw{:,"unit_profit"};
z_CO2 = slice_data_raw{:,"unit_CO2"};

% Identify unique x and y values
uniqueX1 = unique(x1);
uniqueY1 = unique(y1);

% Initialize a matrix for z-values: unit profit
unit_profit_table = NaN(length(uniqueY1), length(uniqueX1));
% Initialize a matrix for z-values: unit profit
unit_CO2_table = NaN(length(uniqueY1), length(uniqueX1));

% Populate the unit profit Matrix with values
for i = 1:length(x1)
    % Find the indices for the current x and y
    xIndex = find(uniqueX1 == x1(i));
    yIndex = find(uniqueY1 == y1(i));
    
    % Assign the unit profit to the correct position in the matrix
    unit_profit_table(yIndex, xIndex) = z_profit(i);
    unit_CO2_table(yIndex, xIndex) = z_CO2(i);
end

% reverse the order of uniqueY to match the flipped zMatrix
uniqueY1 = flipud(uniqueY1(:));
% Flip the zMatrix upside down
unit_profit_table = flipud(unit_profit_table);
unit_CO2_table = flipud(unit_CO2_table);

% the data is well suited for plotting
figure(figure_index);
% plot the heat map
heatmap(uniqueX1, uniqueY1, unit_profit_table, 'Colormap', jet);
clim([0 profit_max]); % Sets the color limits from 0 to 20,000
% Add labels and title if needed
title('Profit per Unit Investment $/(yr-M$-CAPEX): cp_s=3, cr=1');
xlabel('Sand thermal storage volume');
ylabel('Sand thermal storage operation Temperature');
% xlim([0, 5]); % defaul
% ylim([0, 0]); % default
% grid on;  % Optionally, add gridlines

%% Section 2: Optimal PV and battery storage

% for particular slice of cp_s and cr (solar heating collector capactiy and 
% concentation ratio), what is the optimal battery size (investment cost)
% and the optimal PV size

slice_data_raw = readtable('TX_results_summary_stand_alone_sand_PV.csv');

% x-value is the cp: capacity of solar heating collectors
x1 = slice_data_raw{:,"PV"};
% y-value is the cr: concentration ratio
y1 = slice_data_raw{:,"Battery"};

% getting the profit per unit investment
slice_data_raw.unit_profit = slice_data_raw{:,"Total_profit"}./slice_data_raw{:,"CAPEX_total"};
% getting the CO2 per unit investment
slice_data_raw.unit_CO2 = slice_data_raw{:,"Total_CO2"}./slice_data_raw{:,"CAPEX_total"};

% z-value is the unit profit
z_profit = slice_data_raw{:,"unit_profit"};
z_CO2 = slice_data_raw{:,"unit_CO2"};

% Identify unique x and y values
uniqueX1 = unique(x1);
uniqueY1 = unique(y1);

% Initialize a matrix for z-values: unit profit
unit_profit_table = NaN(length(uniqueY1), length(uniqueX1));
% Initialize a matrix for z-values: unit profit
unit_CO2_table = NaN(length(uniqueY1), length(uniqueX1));

% Populate the unit profit Matrix with values
for i = 1:length(x1)
    % Find the indices for the current x and y
    xIndex = find(uniqueX1 == x1(i));
    yIndex = find(uniqueY1 == y1(i));
    
    % Assign the unit profit to the correct position in the matrix
    unit_profit_table(yIndex, xIndex) = z_profit(i);
    unit_CO2_table(yIndex, xIndex) = z_CO2(i);
end

% reverse the order of uniqueY to match the flipped zMatrix
uniqueY1 = flipud(uniqueY1(:));
% Flip the zMatrix upside down
unit_profit_table = flipud(unit_profit_table);
unit_CO2_table = flipud(unit_CO2_table);

% the data is well suited for plotting
figure(figure_index);
% plot the heat map
heatmap(uniqueX1, uniqueY1, unit_profit_table, 'Colormap', jet);
clim([0 profit_max]); % Sets the color limits from 0 to 20,000
% Add labels and title if needed
title('Profit per Unit Investment $/(yr-M$-CAPEX): cp_s=3, cr=1');
xlabel('PV capacity');
ylabel('Battery capacity');
% xlim([0, 5]); % defaul
% ylim([0, 0]); % default
% grid on;  % Optionally, add gridlines



