function [S_star, Recon] = Reconstruct(Recon, max_iterations)

    FULL_RECON_SIZE = Recon.ReconObjects{1}.RECON_SIZE;
    START_RECON_SIZE = Recon.ReconObjects{Recon.NUM_LEVELS}.RECON_SIZE;
    startm = double(rand(START_RECON_SIZE) > 0.5);
    useWeights = Recon.UseNeighborhoodSearchWeights;

    S_star = startm;
    
    % I commented out a large section of code here which implements neighborhood weighting. I will come 
    % back to this later since the math is complex and I don't fully understand it yet. At this point,
    % we will only do unweighted reconstructions.

    [S_star, Recon] = SolidOptimization(Recon, max_iterations, S_star, useWeights);

end 
