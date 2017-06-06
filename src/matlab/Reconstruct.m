function [S_star, Recon] = Reconstruct(Recon, max_iterations, startm)

    FULL_RECON_SIZE = Recon.ReconObjects{1}.RECON_SIZE;
    START_RECON_SIZE = Recon.ReconObjects{Recon.NUM_LEVELS}.RECON_SIZE;
    startm = double(rand(START_RECON_SIZE) > 0.5);
    useWeights = Recon.UseNeighborhoodSearchWeights;
    S_star = startm;

    % We will be plotting a lot of figures during the optimization
    % process, close anything that is open.
    close all;
    
    ll = 1;

    for ii=Recon.NUM_LEVELS:-1:1

        [S_star Recon.ReconObjects{ii}] = SolidOptimization(Recon.ReconObjects{ii}, max_iterations(ii), S_star, useWeights);

        if(ii == 1)
            break;
        end

        % I commented out a large section of code which we will come back to later to perform the multi-resolution reconstruction.

        ll = ll + 1;

    end % For Each Level Of Reconstruction

end
