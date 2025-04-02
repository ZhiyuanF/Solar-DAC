%% This script test the control of the "solar_thermal_control.slx" simulink model

clc
clearvars

% Define the Simulink model name
modelName = 'solar_thermal_control';
% Load the Simulink model
load_system(modelName);

%% Control signal generation

% Time vector
total_time = 3600*24*3;  % Total simulation time in seconds
dt = 1;  % Time step
time = 0:dt:total_time;

%% Cycle reset signal

signal3 = false(size(time));

% Set reset signal to 1 every 3600 seconds
signal3(mod(time, 3600*1.4) == false) = true;


Desorb_start_sig = timeseries(boolean(signal3), time);

% Desorb_start_sig = timeseries(boolean(true(size(time))), time);

%% Solar signal generation

% Time vector for 3 days with 5-min resolution
dt_solar = minutes(5);
t = datetime(2024, 11, 16, 0, 0, 0):dt_solar:datetime(2024, 11, 18, 23, 55, 0);

% Initialize solar profile
solar_profile = zeros(size(t));

% Day 1: Good Solar
day1 = (t >= datetime(2024, 11, 16, 6, 0, 0) & t <= datetime(2024, 11, 16, 18, 0, 0));
solar_profile(day1) = 1000 * sin(pi * hours(t(day1) - datetime(2024, 11, 16, 6, 0, 0)) / 12);

% Day 2: Little Solar (Cloudy)
day2 = (t >= datetime(2024, 11, 17, 6, 0, 0) & t <= datetime(2024, 11, 17, 18, 0, 0));
solar_profile(day2) = 300 * sin(pi * hours(t(day2) - datetime(2024, 11, 17, 6, 0, 0)) / 12);

% Day 3: Normal Solar
day3 = (t >= datetime(2024, 11, 18, 6, 0, 0) & t <= datetime(2024, 11, 18, 18, 0, 0));
solar_profile(day3) = 700 * sin(pi * hours(t(day3) - datetime(2024, 11, 18, 6, 0, 0)) / 12);

% Clip negative values (night time)
solar_profile(solar_profile < 0) = 0;

% Expand each time step to 300 seconds
solar_profile_1s = repelem(solar_profile, 300);

% Create a new time vector with 1-second resolution
t_1s = datetime(2024, 11, 16, 0, 0, 0):seconds(1):datetime(2024, 11, 18, 23, 59, 59);

% Ensure the new time vector length matches the expanded solar profile
t_1s = t_1s(1:length(solar_profile_1s));

% Plot the solar profile
figure;
plot(t_1s, solar_profile_1s);
xlabel('Time');
ylabel('Solar Irradiance (W/m^2)');
title('3-Day Solar Profile (5-min Resolution)');
grid on;

% Assuming t_1s and solar_profile_1s from earlier steps
t_1s_seconds = seconds(t_1s - t_1s(1)); % Convert datetime to seconds relative to start time

% solar thermal correction factor
scale_solar = 0.2;

% Create the timeseries object
Solar_heat_sig = timeseries(solar_profile_1s*scale_solar, t_1s_seconds);

%% Set Up Simulation Inputs in MATLAB:

% Simulate the model and provide the time series to specific inports
simOut = sim(modelName, ...
    'ExternalInput', '[Desorb_start_sig, Solar_heat_sig]');