function [S_star, Recon] = Reconstruct(Recon, max_iterations)

    FULL_RECON_SIZE = Recon.ReconObjects{1}.RECON_SIZE;
    START_RECON_SIZE = Recon.ReconObjects{Recon.NUM_LEVELS}.RECON_SIZE;
    startm = double(rand(START_RECON_SIZE) > 0.5);
    useWeights = Recon.UseNeighborhoodSearchWeights;

    S_star = startm;

    if(useWeights)
    
        for ii=1:Recon.NUM_LEVELS
            if(isfield(Recon.ReconObjects{ii}, 'NBWeights'))
                Recon.ReconObjects{ii} = rmfield(Recon.ReconObjects{ii}, 'NBWeights');
            end
        end
        
        ll = 1;
        
        % We will be plotting a lot of figures during the optimization
        % process, close anything that is open.
        close all;

        % Lets start the reconstruction, we loop through the reconstruction
        % objects from last to first because this is from lowest resolution
        % to highest.
        for ii=Recon.NUM_LEVELS:-1:1

            if(ii == 1)
                useWeights = 0;
            end

            % Make sure to assign the return reconstruction object because
            % this functions is expected to produce side effects on the
            % object.
            [S_star Recon.ReconObjects{ii}] = SolidOptimization(Recon.ReconObjects{ii}, max_iterations(ii), S_star, useWeights);
   
            
            % If we are not at the final resolution, we need to do some
            % things before going into the next resolution. Basically, we
            % need upsample three things, the current reconstruction
            % results, the weight table for texels, and the weight table
            % for neighborhoods. Otherwise, just go to the next iteration
            % (break)
            if(ii == 1)
                break;
            end

            % Upsample the volume. Get the next iterations size. 
            nxSize = Recon.ReconObjects{ii-1}.RECON_SIZE;

            % Upsample the current results, we will use simple nearest
            % neighbor interpolation to avoid blurring. Hopefully, it
            % produces good results.
            S_star = resize(S_star, nxSize, 'nearest');

            % If we are going to the final level next then we need to threshold
            % the input. Lets choose a threshold such that the volume fraction
            % matches the exemplar.
            if(ii == 2)

                % Compute the average volume fraction accross all
                % exemplars.
                vf_to_match = 0;
                for zz=1:length(Recon.EXEMPLARS)
                    T = Recon.EXEMPLARS{zz};
                    vf_to_match = vf_to_match + mean(T(:));
                end
                vf_to_match = vf_to_match / length(Recon.EXEMPLARS);

                % Threshold the image so the volume fraction is as
                % close as we can get.
                S_star = ThresholdToVf(S_star, vf_to_match);

            end

            % Resize the texel weights so we can use them to initialize
            % the next level.
            for pp=1:length(Recon.ReconObjects{ii}.EXEMPLARS)
                Recon.ReconObjects{ii-1}.TexelWeights{pp} = ...
                    resize(Recon.ReconObjects{ii}.TexelWeights{pp}, ... 
                           size(Recon.ReconObjects{ii-1}.EXEMPLARS{pp}), 'nearest');
            end

            % Ok, if we have a neighborhood weighting table, resize it 
            % as well and use it as the input into the next level. Only 
            % do this if we are not going into the final resolution 
            % because we will not use neighborhood weighting at the 
            % final resolution.
            if(useWeights && ii > 2)

                % The neighborhood weight table is not stored like an image but instead
                % a list. We need to make it an image so we can resize\interpolate spatially
                for pp=1:length(Recon.ReconObjects{ii}.EXEMPLARS)
                    weight_table = zeros( size(Recon.ReconObjects{ii}.EXEMPLARS{pp})-Recon.ReconObjects{ii}.NB_SIZE+1 );

                    % Reformat NBWeights{pp} into an image.
                    for kk=1:length(Recon.ReconObjects{ii}.NBWeights{pp})
                        % Get the x,y pixel location of this neighborhood in the exemplar
                        xy = Recon.ReconObjects{ii}.NB_ExemplarLookup{pp}(kk, :);

                        weight_table(xy(1), xy(2)) = Recon.ReconObjects{ii}.NBWeights{pp}(kk);
                    end

                    % Upsample the weight_table with simple nearest neighbor interpolation
                    nx_Ex_Size = size(Recon.ReconObjects{ii-1}.EXEMPLARS{pp})-Recon.ReconObjects{ii}.NB_SIZE+1;
                    weight_table = resize(weight_table, nx_Ex_Size, 'nearest');

                    % Reformat it as a list.
                    Recon.ReconObjects{ii-1}.NBWeights{pp} = zeros(1, size(Recon.ReconObjects{ii-1}.Exemplar_NBs{pp},1));
                    for kk=1:length(Recon.ReconObjects{ii-1}.NBWeights{pp})
                        xy = Recon.ReconObjects{ii-1}.NB_ExemplarLookup{pp}(kk, :);
                        Recon.ReconObjects{ii-1}.NBWeights{pp}(1, kk) = weight_table(xy(1), xy(2));   
                    end

                end

            end
            
            ll = ll + 1;
        
        end % For Each Level Of Reconstruction

    else
        % We have just a single level. We can just call SolidOptimization on
        % it.
        [S_star, Recon] = SolidOptimization(Recon, max_iterations, S_star, useWeights);
    end

end
