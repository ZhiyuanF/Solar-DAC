function [lambda_opts,profits] = LambdaOptimizer_SH(price_data, SH_data, increment, look_ahead, X, k, h, parameters)

    %find how many chunks there are,ceil rounds up last chunk
    chunks = ceil(numel(price_data)/increment);

    %each chunk yields a lambda value and a profit
    lambda_opts = cell(chunks,1);
    profits = cell(chunks,1);

    for chunk = 1:chunks
        
        %find start and end indices for each chunk in price_data
        %including look ahead
        small_price_data_to_opt_start = (chunk-1) * increment + 1;
        small_price_data_to_opt_end = min([chunk * increment + look_ahead , numel(price_data)]); 
        %min solves indexing issues for going beyond the end of price_data
        
        %pull out the range of price_data that we want to optimize over
        small_price_data = price_data(1, small_price_data_to_opt_start : small_price_data_to_opt_end);
        small_SH_data = SH_data(1, small_price_data_to_opt_start : small_price_data_to_opt_end);

        %optimize using SingleLambdaOptimizer function
        [lambda_opt, ~] = SingleLambdaOptimizer_SH(small_price_data, small_SH_data, X, k, h, parameters);
        
        %store lambda value for given chunk + look_ahead
        lambda_opts{chunk} = lambda_opt;
        
        %find start and end indices for each chunk in price_data
        %excluding look ahead
        small_price_data_to_run_start = (chunk-1) * increment + 1;
        small_price_data_to_run_end = min([chunk * increment , numel(price_data)]);
        
        %now that optimization is complete, we use the optimal lambda value to run
        %the simulation and find the intial X and k values for the next chunk
        small_price_data = price_data(1, small_price_data_to_run_start : small_price_data_to_run_end);
        small_SH_data = SH_data(1, small_price_data_to_run_start : small_price_data_to_run_end);
        [X, k, h, profit] = DAC_foriteration_SH(lambda_opt, small_price_data, small_SH_data, X, k, h, parameters);
        
        %store the profit for given chunk
        profits{chunk} = profit;
        
    end
end