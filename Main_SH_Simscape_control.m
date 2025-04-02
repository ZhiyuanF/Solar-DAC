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
    input_solar_file = 'data\TX_solar_testing_input.xlsx';

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
h_hat_norm = 30;
cr = 1; % solar concentration ratio (>=1)

T_min = 300; % (do NOT change) minimum temperature for initiating storage (degree C)
T_design = 400; % Input parameter for designed sand storage target temperature (degree C)

cp_s = cp_s_norm*X_hat; % solar heating charging capacity, the maximum solar charging in MWh per time period
h_hat = h_hat_norm*X_hat; % solar thermal energy storage capacity, the maximum solar storage in MWh
h_effect = h_hat_norm*X_hat*(T_design-T_min)/100; % effective thermal energy storage capacity

eta_heat_transfer = 0.9; % heat transfer efficiency for heating DAC
Ht = 3.15/eta_heat_transfer*X_hat/4; % thermal energy consumption per time step, in MWh per time period
% default: total 3.15 MWh/ton-CO2 energy consumption, in 4 time periods (# of steps in desorption phase)
eta = 0.998; % solar thermal energy storage efficiency, 
% e.g., 0.995 meaning it losses 0.005 of total currently stored energy each time step
% eta = 0.995 means it losses 70% energy in 20 hours

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

%!!! The s_t time series is the new solar heating value input for execution
% of the code !!!

%% Whole system compiling

% parameter compiling
parameters = [X_hat, S, pi_co2, P_a_unit, P_d_unit, P_a, P_d,...
              beta_a_1, beta_a_2, beta_d_1, beta_d_2,...
              cp_s, h_effect, Ht, eta, depth_lower, depth_upper];


%% MAIN: code and function execution

tic

% the "execution" function runs the single case with given data, settings,
% initial conditions and parameters
final_results = execution(price_data, s_t ,increment, look_ahead,...
                          X, k, h, parameters);

disp('Primary code check went through successfully')
toc

%% Results writing to .csv file

full_results_filename = sprintf('results\\%s-X%.2f-s%.2f-h%.2f-T%.2f_full_sand.csv',...
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

%% Get full summary of the results



% initial CAPEX assumptions
CAPEX_DAC = 10.5; % $10.5 million for X_hat=1; MOF technology
CAPEX_Solar_heating = 2; % $2 million for cp_s_norm = 1
CAPEX_thermal_storage = 0.07/3; % $0.03 million for h_hat_norm = 1
% calculate the collector capacity
collect_cap = X_hat*cp_s_norm*cr;
% total CAPEX
CAPEX_total = CAPEX_DAC*X_hat +...
              CAPEX_Solar_heating*collect_cap*(1-0.15)^(log2(collect_cap)) + ...
              CAPEX_thermal_storage*X_hat*h_hat_norm;





%% Simscape model connecting

% Define the Simulink model name
modelName = 'solar_thermal_parameterized';
% Load the Simulink model
load_system(modelName);

%% generating input parameters

sand_specific_heat = 840; % specific heat of sand is 840 J/(K*kg)
sand_storage_tonnage = h_effect*3600*1000000/840/(T_design-T_min)/1000; % mass of sand thermal storage



%% generating input time series signals

% Getting the solar input signal
% Create the timeseries object
Solar_heat_sig_5min = s_t*cp_s; % solar charging at 5 min resolution, MWh/5min
% Expand data from 5-min resolution to 1-second resolution
Solar_heat_sig_1sec_uncorrected = repelem(Solar_heat_sig_5min, 300);
% divide the 5-min energy into 300 second
Solar_heat_sig_1sec_corrected = Solar_heat_sig_1sec_uncorrected/300;
% perferm MWh to Joule unit conversion
Solar_heat_sig_double = Solar_heat_sig_1sec_corrected*3600*1000000;
Solar_heat_sig_double = Solar_heat_sig_double'; % transpose

% Getting the desorb initiation signal
full_result_raw = readmatrix(full_results_filename);
DAC_status_raw = full_result_raw(:,7);
% Expand data from 5-min resolution to 1-second resolution
DAC_status_1sec_uncorrected = repelem(DAC_status_raw, 300);
% DAC_status_1sec_uncorrected = cellfun(@double, DAC_status_1sec_uncorrected);
% Convert cell array to numeric array
% DAC_status_1sec_uncorrected = cell2mat(DAC_status_1sec_uncorrected);
% Initialize the output signal with zeros
Desorb_start_sig_double = false(size(DAC_status_1sec_uncorrected));
% Identify transitions from 1 to 0
Desorb_start_sig_double(2:end) = (DAC_status_1sec_uncorrected(2:end) == false) & (DAC_status_1sec_uncorrected(1:end-1) == true);


% Time vector
total_time = 3600*24*365;  % Total simulation time in seconds, 365 days 
dt = 1;  % Time step
time = 1:dt:total_time;
% reset signal at very first time step shall be >600
reset_signal_raw = zeros(size(time));
reset_signal_raw(1:300) = 640; % random number >600
reset_sig_double = reset_signal_raw;

% get the maximum of solar signal to calculate capacity factor
solar_sig_max = max(Solar_heat_sig_double);

%% Control signal generation



%{   
% Plot the solar profile
figure;
plot(time(1:3600*24), Desorb_start_sig_double(1:3600*24));
xlabel('Time');
ylabel('Desorbtion initiation signal');
title('Desorb Signal Check');
grid on;

% Extract the first 3600*24 data points
num_points = 3600 * 24; % Define the number of points to plot
time_to_plot = Desorb_start_sig.Time(1:num_points);
data_to_plot = Desorb_start_sig.Data(1:num_points);
% Plot the timeseries
figure;
plot(time_to_plot, data_to_plot, '-o'); % '-o' for line with markers (optional)
xlabel('Time');
ylabel('Signal');
title('First 3600*24 Data Points of Timeseries');
grid on;
%}

Desorb_start_sig = timeseries(boolean(Desorb_start_sig_double), time);
Solar_heat_sig = timeseries(Solar_heat_sig_double, time);
Reset_sig = timeseries(reset_sig_double, time);




%% Set Up Simulation Inputs in MATLAB:

% Simulate the model and provide the time series to specific inports
disp('Runnig MATLAB SimScape model')
simOut = sim(modelName, ...
    'ExternalInput', '[Desorb_start_sig, Solar_heat_sig]');
toc