function [best_lambda_opt, best_prof] = SingleLambdaOptimizer_SH(price_data, SH_data, X, k, h, parameters)

    % define the function that wraps the DAC_Cycle function
    objectiveFunc = @(lambda)(-DAC_foropt_SH(lambda, price_data, SH_data, X, k, h, parameters));
    
    % initialize variables
    best_prof = 0; % initialize with zero
    best_lambda_opt = 0; % initialize with zero
    counter = 0; % initialize with zero
    
    series = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16,...
        18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50,...
        50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,...
        250, 300, 350, 400, 450, 500];
    
    for lambda_guess = series
                
        % perform optimization using fminsearch
        [lambda_opt, fake_profit] = fminsearch(objectiveFunc, lambda_guess);
    
        % update the best optimized output and corresponding lambda_opt
        if -fake_profit > best_prof
            best_prof = -fake_profit;
            best_lambda_opt = lambda_opt;
            counter = counter + 1;
        end
    end

    %double check if the best lambda opt has ever been changed, if not, use
    %lowest price value for the year to ensure idle
    
    if counter == 0
            best_lambda_opt = min(price_data)-1;
    end
        
    % once lambda_opt has been found, take the unadjusted profit from that run
    [fake_prof, boost] = DAC_foropt_SH(best_lambda_opt, price_data, SH_data, X, k, h, parameters);
    best_prof = fake_prof-boost;
end