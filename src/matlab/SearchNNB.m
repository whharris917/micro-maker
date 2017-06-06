function [Recon] = SearchNNB(S, Recon, useWeights)

    for ExIndex=1:size(Recon.EXEMPLARS, 1)
    
        nbOffsets = squeeze(Recon.nbOffsets(ExIndex, :, :));
        [NB_queries, isPeriodic] = Extract3DNeighborhoods(S, nbOffsets, Recon.NUM_CORES);
        
        [Xnnidx] = flann_search(Recon.Exemplar_Index{ExIndex}, NB_queries', 1, Recon.Exemplar_Params{ExIndex});

        % Now, mark any neighborhoods that are on the edges as a -1, this
        % will tell the code later to ignore these guys in our calculations
        Xnnidx(isPeriodic) = -1;
        
        % We need to reshape the search results into the 4D array that the
        % rest of our code expects. We need to do the permute because of
        % how we arranged the neighborhoods when extracting them from the
        % 3D image. We did it in the opposite of the MATLAB convention for
        % some reason, maybe because they are extracted in C++ code.
        Recon.NNB_Table(:,:,:,ExIndex) = permute(reshape(Xnnidx, size(S)), [3 2 1]);
        
    end

end
