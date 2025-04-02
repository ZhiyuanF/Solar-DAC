clc
clearvars

%% initiating the calculation, 5-min temporal resolution
%Set up length of simulation
steps_per_day = 288; % 5-min temporal resolution
days = 365; % 1 year annual data/results

duration = steps_per_day * days;
startingday = 1;
start_step = steps_per_day * (startingday - 1) + 1;

%% input the controling parameters for the simulation

pi_co2 = 200; % selling co2 at pi_co2 per ton

% adjust the optmization and look-ahead horizon based on cycle length
increment = 288;
look_ahead = 288;

% cycle capacity definition
X_hat = 1;

% cycle length definition
beta_a_1 = 0.2*X_hat;
beta_a_2 = -0.2;
beta_d_1 = 0*X_hat;
beta_d_2 = 0.4;

% DAC tech parameters
S = 40*X_hat; % material cost (cycle cost)
P_a_unit = 0.080250; % adsorption energy consumption
P_d_unit = 0.016167; % desorption energy consumption
P_a = P_a_unit*X_hat;
P_d = P_d_unit*X_hat;

% DAC normalized cycle depth
depth_lower = 0.1;
depth_upper = 0.75;



%% Initial conditions
X = 0;
k = 0;
h = 0;

%% Geography selection, data loading and processing

disp('Note: Double-check input files and paraters before running');
user_input1 = input(['Enter "NY" for New York, "CA" for California,' ...
                    ' or "TX" for Texas to load emission data: '], 's');

% New York, California, or Texas User Input for Emission/Price/Solar Data
if strcmpi(user_input1, 'NY')
    
    % Load emissions data from CSV for New York
    myISO = readtable('data\NYISO_full_filtered_emission.csv');
    input_mat_file = 'data\expanded_RAW_NYISO.mat';
    input_solar_file = 'data\TX_solar_testing_input.xlsx';
    
elseif strcmpi(user_input1, 'CA')
    
    % Load emissions data from CSV for California
    myISO = readtable('data\CAISO_full_filtered_emission.csv');
    input_mat_file = 'data\expanded_RAW_CAISO.mat';
    input_solar_file = 'data\CA_solar_testing_input.xlsx';

elseif strcmpi(user_input1, 'TX')
    
    % Load emissions data from CSV for Texas
    myISO = readtable('data\TX_full_filtered_emission.csv');
    input_mat_file = 'data\expanded_RAW_TX.mat';
    input_solar_file = 'data\TX_solar_testing_input.xlsx';
      

else
    error('Invalid input. Please enter "NY", "CA", or "TX".');
end

% CO2 data extration and processing
co2_data = transpose(myISO{1:end, 11}); % Extract CO2 emissions from the 11th column

% pricing data extration and processing
fileln = load(input_mat_file);
myRTP = fileln.RTP;
price_data = myRTP(start_step:start_step+duration-1); % Load prices for the specified time chunk
myRTP = price_data; % Assign the transposed price data to myRTP

% solar data extraction and processing
dataRange = 'A2:C105121'; % skip the first row where header is read
SH_raw = readmatrix(input_solar_file,'Range', dataRange); % read the numeric data only
mySH = SH_raw(:,2); % the 3rd column is the "Solar Heating Binary" input
% renaming
SH_data = transpose(mySH);

%% Solar heating/storage system configuration

% Solar heating and thermal storage parameters
cp_s_norm = 3;
h_hat_norm = 70;

T_min = 300; % (do NOT change) minimum temperature for initiating storage (degree C)
T_design = 400; % Input parameter for designed sand storage target temperature (degree C)

cp_s = cp_s_norm*X_hat; % solar heating charging capacity, the maximum solar charging in MWh per time period
h_hat = h_hat_norm*X_hat; % solar thermal energy storage capacity, the maximum solar storage in MWh
h_effect = h_hat_norm*X_hat*(T_design-T_min)/100; % effective thermal energy storage capacity in MWh

eta_heat_transfer = 0.9; % heat transfer efficiency for heating DAC
Ht = 3.15/eta_heat_transfer*X_hat/4; % thermal energy consumption per time step, in MWh per time period
% default: total 3.15 MWh/ton-CO2 energy consumption, in 4 time periods (# of steps in desorption phase)
eta_thermal_loss = 0.998; % solar thermal energy storage efficiency, 
% e.g., 0.995 meaning it losses 0.005 of total currently stored energy each time step
% eta = 0.995 means it losses 70% energy in 20 hours

% solar system design parameters
if h_hat_norm > 0
    Temp = T_design; % target system temperature (degree C), 500 for heating with molten salt storage
    eta_storage = 0.95; % when thermal storage is included, extra loss will occur
elseif h_hat_norm == 0
    Temp = 150; % target system temperature (degree C), 150 for heating with storage
    eta_storage = 1; % without thermal storage, by-pass the storage loss with efficiency 1
end
cr = 1; % solar concentration ratio (>=1)

% calculate the solar heating collector efficiency
eta_c_t = 0.78 - 8.8e-7 * Temp^(2) * (1.1./(SH_data + 0.1)) ./cr;
eta_c_t(eta_c_t < 0) = 0; % correct any negative efficiency to 0
% calculate the available heat
s_t = SH_data .* eta_c_t * eta_storage; % solar capacity factor * solar heating collector efficiency

%!!! The s_t time series is the new solar heating value input for execution
% of the code !!!

%% Whole system compiling

% parameter compiling
parameters = [X_hat, S, pi_co2, P_a_unit, P_d_unit, P_a, P_d,...
              beta_a_1, beta_a_2, beta_d_1, beta_d_2,...
              cp_s, h_effect, Ht, eta_thermal_loss, depth_lower, depth_upper];


%% MAIN: code and function execution

tic

% the "execution" function runs the single case with given data, settings,
% initial conditions and parameters
final_results = execution(price_data, s_t ,increment, look_ahead,...
                          X, k, h, parameters);

disp('Primary code check went through successfully')
toc

%% Results writing to .csv file

full_results_filename = sprintf('%s-X%.2f-s%.2f-h%.2f-T%.2f_full_sand.csv',...
                        user_input1, X_hat, cp_s, h_hat, T_design);
% filename for sequence printing
fileID = fopen(full_results_filename, 'w');
% column names to the sequence CSV file
fprintf(fileID, 'step,chunk,optimal lambda,u,v,SHB,k,z,L,SolarInput,h,X,a,d,profit\n');

% loop through the cell array and write each row as a CSV line
for row = 1:size(final_results, 1)
    fprintf(fileID, '%f,%d,%f,%.0f,%.0f,%.0f,%.0f,%.0f,%f,%f,%f,%f,%f,%f,%f\n', final_results{row,:});
end

% close the file
fclose(fileID);

disp('Primary code full results written')
toc

%{

%% impact of pi_co2 selling price

% initial CAPEX assumptions
CAPEX_DAC = 10.5; % $10.5 million for X_hat=1; MOF technology
CAPEX_Solar_heating = 2; % $2 million for cp_s_norm = 1
CAPEX_thermal_storage = 0.01; % $0.01 million
% calculate the collector capacity
collect_cap = X_hat*cp_s_norm*cr;
% total CAPEX
CAPEX_total = CAPEX_DAC*X_hat +...
              CAPEX_Solar_heating*collect_cap*(1-0.15)^(log2(collect_cap)) + ...
              CAPEX_thermal_storage*X_hat*h_hat_norm;

case_pi_co2 = 20:20:500;
case_cp_s = 3;
case_h_hat = 10:5:100;
case_cr = 1;
case_T = 400;

case_number = numel(case_pi_co2)*numel(case_cp_s)*numel(case_h_hat)*numel(case_cr)*numel(case_T);
case_counter = 0;

% create space for storing result summary
summary_columns = 10;
% the 10 columns are: (1) pi_co2 (2) cp_s (3) h_hat (4) cr (5) T_target 
% (6) CF (7) SHB_CF (8) CAPEX (9) Total_profit (10) Total_CO2
summary_results = cell(numel(case_number), summary_columns);

% construct each case in for loops
for pi_co2_input = case_pi_co2
    for cp_s_input = case_cp_s
        for h_hat_input = case_h_hat
            for cr_input = case_cr
                for T_input = case_T
                    % update the parameters using inputs
                    pi_co2 = pi_co2_input; % selling co2 at pi_co2 per ton
                    X_hat = 1;
                    % Solar heating and thermal storage parameters
                    cp_s_norm = cp_s_input;
                    h_hat_norm = h_hat_input;
                    cr = cr_input;
                    % cycle length definition
                    beta_a_1 = 0.2*X_hat;
                    beta_a_2 = -0.2;
                    beta_d_1 = 0*X_hat;
                    beta_d_2 = 0.4;
                    
                    % DAC tech parameters
                    S = 40*X_hat; % material cost (cycle cost)
                    P_a_unit = 0.080250; % adsorption energy consumption
                    P_d_unit = 0.016167; % desorption energy consumption
                    P_a = P_a_unit*X_hat;
                    P_d = P_d_unit*X_hat;

                    T_min = 300; % (do NOT change) minimum temperature for initiating storage (degree C)
                    T_design = T_input; % Input parameter for designed sand storage target temperature (degree C)
                    
                    % solar thermal inputs
                    cp_s = cp_s_norm*X_hat; % solar heating charging capacity, the maximum solar charging in MWh per time period
                    h_hat = h_hat_norm*X_hat; % solar thermal energy storage capacity, the maximum solar storage in MWh
                    h_effect = h_hat_norm*X_hat*(T_design-T_min)/100; % effective thermal energy storage capacity
                    Ht = 3.15*X_hat/4; % thermal energy consumption per time step, in MWh per time period
                    % default: total 3.15 MWh/ton-CO2 energy consumption, in 4 time periods (# of steps in desorption phase)
                    eta = 0.995; % solar thermal energy storage efficiency, 
                    % e.g., 0.995 meaning it losses 0.005 of total currently stored energy each time step
                    % eta = 0.995 means it losses 70% energy in 20 hours
        
                    depth_lower = 0.1;
                    depth_upper = 0.75;
                    % solar system design parameters
                    if h_hat_norm > 0
                        Temp = T_design; % target system temperature (degree C), 500 for heating with storage
                        eta_storage = 0.95; % when thermal storage is included, extra loss will occur
                    elseif h_hat_norm == 0
                        Temp = 150; % target system temperature (degree C), 150 for heating with storage
                        eta_storage = 1; % without thermal storage, by-pass the storage loss with efficiency 1
                    end
    
                    % calculate the solar heating collector efficiency
                    eta_c_t = 0.78 - 8.8e-7 * Temp^(2) * (1.1./(SH_data + 0.1)) ./cr;
                    eta_c_t(eta_c_t < 0) = 0; % correct any negative efficiency to 0
                    % calculate the available heat
                    s_t = SH_data .* eta_c_t * eta_storage; % solar capacity factor * solar heating collector efficiency
        
                    % parameter compiling
                    parameters = [X_hat, S, pi_co2, P_a_unit, P_d_unit, P_a, P_d,...
                                beta_a_1, beta_a_2, beta_d_1, beta_d_2,...
                                cp_s, h_effect, Ht, eta_thermal_loss, depth_lower, depth_upper];
        
                    % redefine initial conditions
                    X = 0;
                    k = 0;
                    h = 0;
                    
                    % simulation execution
                    final_results = execution(price_data, s_t ,increment, look_ahead,...
                                  X, k, h, parameters);
        
                    % counter
                    case_counter = case_counter + 1;
                    % calculate the collector capacity
                    collect_cap = X_hat*cp_s_norm*cr;
                    % write results
                    CF = sum(sum(cell2mat(final_results(:, 4:5))))/(2*duration); % total capacity factor
                    SHB_CF = sum(cell2mat(final_results(:, 6)))/duration; % total SHB_CF
                    CAPEX_total = CAPEX_DAC*X_hat + ...
                                  CAPEX_Solar_heating*collect_cap*(1-0.15)^(log2(collect_cap)) + ...
                                  CAPEX_thermal_storage*h_hat; % total CAPEX
                    Total_profit = sum(cell2mat(final_results(:, 15)));
                    Total_CO2 = sum(cell2mat(final_results(:, 14)));
        
                    % results storing
                    summary_results(case_counter, :) = {pi_co2, cp_s, h_hat, cr, T_design, CF,...
                                            SHB_CF, CAPEX_total, Total_profit, Total_CO2};
        
                    fprintf('Finished %.0f out of %.0f.\n', case_counter, case_number);
                end
            end
        end
    end
end

% Summary Results writing to .csv file

summary_results_filename = sprintf('%s_pi_co2_summary_sand.csv',...
                        user_input1);
% filename for sequence printing
fileID = fopen(summary_results_filename, 'w');
% column names to the sequence CSV file
fprintf(fileID, 'pi_co2, cp_s, h_hat, cr, T_design, CF,SHB_CF, CAPEX_total, Total_profit, Total_CO2\n');

% loop through the cell array and write each row as a CSV line
for row = 1:size(summary_results, 1)
    fprintf(fileID, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', summary_results{row,:});
end

% close the file
fclose(fileID);

disp('Primary code full results written')
toc


%}



%{

%% Get full summary of the results

% initial CAPEX assumptions
CAPEX_DAC = 10.5; % $10.5 million for X_hat=1; MOF technology
CAPEX_Solar_heating = 2; % $2 million for cp_s_norm = 1
CAPEX_thermal_storage = 0.01; % $0.01 million
% calculate the collector capacity
collect_cap = X_hat*cp_s_norm*cr;
% total CAPEX
CAPEX_total = CAPEX_DAC*X_hat +...
              CAPEX_Solar_heating*collect_cap*(1-0.15)^(log2(collect_cap)) + ...
              CAPEX_thermal_storage*X_hat*h_hat_norm;

% run the system simulation with multiple cases and store results

case_X = 0.2:0.2:3;
case_cp_s = 3;
case_h_hat = 70;
case_cr = 1;
case_T = 400;

case_number = numel(case_X)*numel(case_cp_s)*numel(case_h_hat)*numel(case_cr)*numel(case_T);
case_counter = 0;

% create space for storing result summary
summary_columns = 10;
% the 10 columns are: (1) X_hat (2) cp_s (3) h_hat (4) cr (5) T_target 
% (6) CF (7) SHB_CF (8) CAPEX (9) Total_profit (10) Total_CO2
summary_results = cell(numel(case_number), summary_columns);

% construct each case in for loops
for X_input = case_X
    for cp_s_input = case_cp_s
        for h_hat_input = case_h_hat
            for cr_input = case_cr
                for T_input = case_T
                    % update the parameters using inputs
                    X_hat = X_input;
                    % Solar heating and thermal storage parameters
                    cp_s_norm = cp_s_input;
                    h_hat_norm = h_hat_input;
                    cr = cr_input;
                    % cycle length definition
                    beta_a_1 = 0.2*X_hat;
                    beta_a_2 = -0.2;
                    beta_d_1 = 0*X_hat;
                    beta_d_2 = 0.4;
                    
                    % DAC tech parameters
                    S = 40*X_hat; % material cost (cycle cost)
                    P_a_unit = 0.080250; % adsorption energy consumption
                    P_d_unit = 0.016167; % desorption energy consumption
                    P_a = P_a_unit*X_hat;
                    P_d = P_d_unit*X_hat;

                    T_min = 300; % (do NOT change) minimum temperature for initiating storage (degree C)
                    T_design = T_input; % Input parameter for designed sand storage target temperature (degree C)
                    
                    % solar thermal inputs
                    cp_s = cp_s_norm*X_hat; % solar heating charging capacity, the maximum solar charging in MWh per time period
                    h_hat = h_hat_norm*X_hat; % solar thermal energy storage capacity, the maximum solar storage in MWh
                    h_effect = h_hat_norm*X_hat*(T_design-T_min)/100; % effective thermal energy storage capacity
                    Ht = 3.15*X_hat/4; % thermal energy consumption per time step, in MWh per time period
                    % default: total 3.15 MWh/ton-CO2 energy consumption, in 4 time periods (# of steps in desorption phase)
                    eta = 0.995; % solar thermal energy storage efficiency, 
                    % e.g., 0.995 meaning it losses 0.005 of total currently stored energy each time step
                    % eta = 0.995 means it losses 70% energy in 20 hours
        
                    depth_lower = 0.1;
                    depth_upper = 0.75;
                    % solar system design parameters
                    if h_hat_norm > 0
                        Temp = T_design; % target system temperature (degree C), 500 for heating with storage
                        eta_storage = 0.95; % when thermal storage is included, extra loss will occur
                    elseif h_hat_norm == 0
                        Temp = 150; % target system temperature (degree C), 150 for heating with storage
                        eta_storage = 1; % without thermal storage, by-pass the storage loss with efficiency 1
                    end
    
                    % calculate the solar heating collector efficiency
                    eta_c_t = 0.78 - 8.8e-7 * Temp^(2) * (1.1./(SH_data + 0.1)) ./cr;
                    eta_c_t(eta_c_t < 0) = 0; % correct any negative efficiency to 0
                    % calculate the available heat
                    s_t = SH_data .* eta_c_t * eta_storage; % solar capacity factor * solar heating collector efficiency
        
                    % parameter compiling
                    parameters = [X_hat, S, pi_co2, P_a_unit, P_d_unit, P_a, P_d,...
                                beta_a_1, beta_a_2, beta_d_1, beta_d_2,...
                                cp_s, h_effect, Ht, eta_thermal_loss, depth_lower, depth_upper];
        
                    % redefine initial conditions
                    X = 0;
                    k = 0;
                    h = 0;
                    
                    % simulation execution
                    final_results = execution(price_data, s_t ,increment, look_ahead,...
                                  X, k, h, parameters);
        
                    % counter
                    case_counter = case_counter + 1;
                    % calculate the collector capacity
                    collect_cap = X_hat*cp_s_norm*cr;
                    % write results
                    CF = sum(sum(cell2mat(final_results(:, 4:5))))/(2*duration); % total capacity factor
                    SHB_CF = sum(cell2mat(final_results(:, 6)))/duration; % total SHB_CF
                    CAPEX_total = CAPEX_DAC*X_hat + ...
                                  CAPEX_Solar_heating*collect_cap*(1-0.15)^(log2(collect_cap)) + ...
                                  CAPEX_thermal_storage*h_hat; % total CAPEX
                    Total_profit = sum(cell2mat(final_results(:, 15)));
                    Total_CO2 = sum(cell2mat(final_results(:, 14)));
        
                    % results storing
                    summary_results(case_counter, :) = {X_hat, cp_s, h_hat, cr, T_design, CF,...
                                            SHB_CF, CAPEX_total, Total_profit, Total_CO2};
        
                    fprintf('Finished %.0f out of %.0f.\n', case_counter, case_number);
                end
            end
        end
    end
end

%% Summary Results writing to .csv file

summary_results_filename = sprintf('%s_results_summary_sand.csv',...
                        user_input1);
% filename for sequence printing
fileID = fopen(summary_results_filename, 'w');
% column names to the sequence CSV file
fprintf(fileID, 'X_hat, cp_s, h_hat, cr, T_design, CF,SHB_CF, CAPEX_total, Total_profit, Total_CO2\n');

% loop through the cell array and write each row as a CSV line
for row = 1:size(summary_results, 1)
    fprintf(fileID, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', summary_results{row,:});
end

% close the file
fclose(fileID);

disp('Primary code full results written')
toc


%}
