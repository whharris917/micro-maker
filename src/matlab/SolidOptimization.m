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
        
        PlotIteration(Recon, S, Snew, useWeights);
          
        % It converged!
        if( ii == NUM_ITERATIONS )
            fprintf(1, 'Done.\nIt converged!\n');
            break;
        else
            S = Snew;
        end

        fprintf(1, 'Done: %f, vf=%f seconds', toc);
        elapsed_opt = elapsed_opt + toc;

    end
    
    fprintf(1, '   Iterations: %d   SearchTime: %f secs   OptimTime: %f secs\n', ii, elapsed_search, elapsed_opt);

end
