function [final_results] = execution(price_data, SH_data ,increment, look_ahead, X, k, h, parameters)
    % Run LambdaOptimizer, get the optimal lambdas and profits for each chunk for
    % the whole data set
    [lambda_opts] = LambdaOptimizer_SH(price_data, SH_data ,increment, look_ahead, X, k, h, parameters);
    
    % create a cell with
    % # of rows = # of price data points
    % # of columns = 15
    total_columns = 15;
    
    % the 13 rows are (1) time steps (2) chunks (3) optimal lambda for chunk 
    % (4) u (5) v (6) SHB (7) k (8) z (9) L (10) SolarInput (11) h (12) X 
    % (13) a (14) d (15) profit
    final_results = cell(numel(price_data), total_columns);
    final_results(:,1) = num2cell(1:numel(price_data));
    
    column_chunk = 2;
    column_op_lambda = 3;
    
    % loop through each chunk
    for chunk = 1:numel(lambda_opts)
        
        % length of the chunk is increment, except for at the very end of the
        % data set, where min rounds down
        chunk_len = min([increment, numel(price_data) - (chunk-1) * increment]);
        
        %   every item of chunk column is = chunk number
        %   every item of optimal lambda column is the optimal lambda for that chunk
        c1 = cell(chunk_len,1);
        c1(:) = {chunk};
        c2 = cell(chunk_len,1);
        c2(:) = {lambda_opts{chunk}};
            
        % take the # of chunks that you're on, multiply by increment and then
        % add 1 to get the start point
        % take the # of chunks that you're on, multiply by increment and then
        % add the chunk length to get the end point
        final_results((chunk-1) * increment +1 :...
            (chunk-1) * increment + chunk_len, column_chunk) = c1(:);
        final_results((chunk-1) * increment +1 :...
            (chunk-1) * increment + chunk_len, column_op_lambda) = c2(:);
    
        % retrieve right start and end point for the chunk in price_data and
        % SH_data
        small_price_data = price_data(1, (chunk-1) * increment +1 : (chunk-1) * increment + chunk_len);
        small_SH_data = SH_data(1, (chunk-1) * increment +1 : (chunk-1) * increment + chunk_len);
        
        % the rest of columns from DAC_fordata function that runs on whatever the
        % current optimal lambda value is for that chunk
        
        % note: this is what takes the longest, it's not the optimization
        % it's the data writing
        final_results((chunk-1) * increment +1 :...
            (chunk-1) * increment + chunk_len, column_op_lambda+1:total_columns) = ...
            DAC_fordata_SH(lambda_opts{chunk}, small_price_data, small_SH_data, X, k, h, parameters);
        
        % update X, k, and h at the beginning of the next chunk
        [X, k, h] = DAC_foriteration_SH(lambda_opts{chunk}, small_price_data, small_SH_data, X, k, h, parameters);
    end

end

