function [S, Recon, S_iters] = SolidOptimization(Recon, NUM_ITERATIONS, start, useWeights)

    S = start;
    
    % If the user has asked to return a copy of each iteration's structure
    % then we need to allocate it.
    if(nargout > 2)
        S_iters = zeros([size(S) NUM_ITERATIONS+1]);
        S_iters(:, :, :, 1) = S;
    end
    
    elapsed_search = 0;
    elapsed_opt = 0;

    for ii=1:NUM_ITERATIONS
        fprintf(1, 'Iteration: %d\n', ii);

        tic;
        fprintf(1, '    Finding nearest neighbors ... ');
        Recon = SearchNNB(S, Recon, useWeights);
        elapsed_search = elapsed_search + toc;
        fprintf(1, 'Done: %f seconds\n', toc);

        tic;
        fprintf(1, '    Optimizing solid ... ');
        [Snew, Recon.TexelWeights, Recon.SourceTable] = OptimizeSolidMEX(S, Recon);

        % For each exemplar, lets calculate a histogram that keeps track
        % of how often that neighborhood has been matched. This will give
        % us and idea of how much of the exemplar is being used in each
        % iteration
        Recon.NBHoodHist = cell(length(Recon.EXEMPLARS),1);
        for pp=1:size(Recon.EXEMPLARS, 1)
            % Calculate a histogram that counts many times each
            % neighborhood has been selected as the best matching
            % neighborhood.
            H = zeros(size(Recon.EXEMPLARS{pp})); 
            NBs = Recon.NNB_Table(:, :, :, pp); 
            NBs = NBs(:);
            NBs = NBs(NBs > 0);
            for kk=1:length(NBs) 
                xy = Recon.NB_ExemplarLookup{pp}(NBs(kk), :); 
                H(xy(1), xy(2)) = H(xy(1), xy(2)) + 1; 
            end
            H = H ./ sum(H(:));
            
            Recon.NBHoodHist{pp} = H;
        end


        PlotIteration(Recon, S, Snew, useWeights);
       
        % Store a copy of the current structure every iteration if the user
        % wants it.
        if(nargout > 2)
            S_iters(:, :, :, ii+1) = Snew;
        end

        % Calculate the percentage of the image that has changed
        perVoxelPercentChange = abs(S(:) - Snew(:)) ./ S(:); 
        
        % If we have a NaN, that means 0 / 0, which we will say is no
        % relative change.
        perVoxelPercentChange(isnan(perVoxelPercentChange)) = 0;
        
        % If we have Inf, then we have something like x/0. Which we will
        % just say is 1.
        perVoxelPercentChange(isinf(perVoxelPercentChange)) = 1;
        
        % Calculate the total percentage change over the entire image.
        % This is the relative percent change divided by the total possible
        % change.
        percentChange = 100*(sum(perVoxelPercentChange) / numel(S));
          
        % It converged!
        if( ii == NUM_ITERATIONS )
            fprintf(1, 'Done.\nIt converged!\n');
            break;
        else
            S = Snew;
        end

        fprintf(1, 'Done: %f, vf=%f seconds, percentChange=%f\n', toc, mean(S(:)), percentChange);
        elapsed_opt = elapsed_opt + toc;

    end
    fprintf(1, '   Iterations: %d   SearchTime: %f secs   OptimTime: %f secs\n', ii, elapsed_search, elapsed_opt);

    % Resize the iterations matrix if they want it. We might not need it 
    % all if the reconstruciton converged earlier that NUM_ITERATIONS
    if(nargout > 2)
        S_iters = S_iters(:, :, :, 1:ii);
    end

end
