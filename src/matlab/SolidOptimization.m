function [S, Recon] = SolidOptimization(Recon, NUM_ITERATIONS, start, useWeights)

    S = start;
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

        Recon.NBHoodHist = cell(length(Recon.EXEMPLARS),1);
        for pp=1:size(Recon.EXEMPLARS, 1)
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

        perVoxelPercentChange = abs(S(:) - Snew(:)) ./ S(:); 
        perVoxelPercentChange(isnan(perVoxelPercentChange)) = 0;
        perVoxelPercentChange(isinf(perVoxelPercentChange)) = 1;
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

end
