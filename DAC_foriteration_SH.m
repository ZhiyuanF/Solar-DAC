function [X, k, h, PROFIT] = DAC_foriteration_SH(price_cap, price_data, SH_data, X, k, h, parameters)


    % Fixed parameters
    X_hat = parameters(1,1); % Paying S dollar per each cycle for material/operation cost
    S = parameters(1,2); % Paying S dollar per each cycle for material/operation cost
    pi_co2 = parameters(1,3); % Selling co2 at pi_co2 per ton, incentive
    P_a = parameters(1,6); % Correction of absorption power consumption by the installed capacity
    P_d = parameters(1,7); % Correction of desorption power consumption by the installed capacity
    beta_a_1 = parameters(1,8); % Piecewise linear approximation for absorption, first-order term
    beta_a_2 = parameters(1,9); % Piecewise linear approximation for absorption, second-order term
    beta_d_1 = parameters(1,10); % Piecewise linear approximation for desorption, first-order term
    beta_d_2 = parameters(1,11); % Piecewise linear approximation for desorption, second-order term

    % Solar heating and thermal storage fixed parameter
    cp_s = parameters(1,12); % solar heating charging capacity, the maximum solar charging in MWh per time period
    h_hat = parameters(1,13); % solar thermal energy storage capacity, the maximum solar storage in MWh
    Ht = parameters(1,14); % thermal energy consumption per time step, in MWh per time period
    eta = parameters(1,15); % solar thermal energy storage efficiency

    % cycling depth parameters
    depth_lower = parameters(1,16);
    depth_upper = parameters(1,17);

    % Initialize variables
    PROFIT = 0;
    %price_cap = 100000000000;

    for t = 1:numel(price_data)

        h = h*eta + SH_data(t)*cp_s; % calculate the solar storage after charging
        SHB = (h > Ht); % check whether solar storage is enough for desorption

        if price_data(t) > price_cap  %if price is higher than lambda, idle
            u = false;
            v = false;
            z = false;

        elseif X > depth_upper*X_hat %if saturation is high, desorb
            u = false;
            v = SHB; % whether switch to adsorb or not depend on SHB result
            z = false;
            k = ~v; % if switch, change k=0 (v=1), if not, keep k=1 (v=0)

        elseif X < depth_lower*X_hat %if saturation is low, adsorb

            % check if most recently desorbing,
            % and if so, indicate a new cycle has begun
            u = true;
            v = false;
            z = ~k;
            k = true;
            % switch to adsorption does NOT depend on SH input

        else
            % continue adsorbing or desorbing as previously
            u = k;
            v = logical((~k)*(SHB)); % desorbing or not depending on SHB result
            z = false;

        end
        
        % calculating adsorption and desorption rates, state of saturation,
        % tons of CO2 sold, and profit
        a = (beta_a_1 + beta_a_2 * X) * u;
        d = (beta_d_1 + beta_d_2 * X) * v;
        X = X + a - d;
        % calculate the thermal storage capacity
        h = min(h - v*Ht, h_hat); % thermal storage does not exceed capacity

        profit = (pi_co2 * d) - (price_data(t) * (P_a * u + P_d * v) + S * z);
        PROFIT = PROFIT+profit;
        
        %no boost, no optimization happening
    end
end