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

    for ii=Recon.NUM_LEVELS:-1:1

        if(ii == 1)
            useWeights = 0;
        end

        [S_star Recon.ReconObjects{ii}] = SolidOptimization(Recon.ReconObjects{ii}, max_iterations(ii), S_star, useWeights);

        if(ii == 1)
            break;
        end

        nxSize = Recon.ReconObjects{ii-1}.RECON_SIZE;

        S_star = resize(S_star, nxSize, 'nearest');

        if(ii == 2)

            vf_to_match = 0;
            for zz=1:length(Recon.EXEMPLARS)
                T = Recon.EXEMPLARS{zz};
                vf_to_match = vf_to_match + mean(T(:));
            end
            vf_to_match = vf_to_match / length(Recon.EXEMPLARS);

            S_star = ThresholdToVf(S_star, vf_to_match);

        end

        for pp=1:length(Recon.ReconObjects{ii}.EXEMPLARS)
            Recon.ReconObjects{ii-1}.TexelWeights{pp} = ...
                resize(Recon.ReconObjects{ii}.TexelWeights{pp}, ... 
                       size(Recon.ReconObjects{ii-1}.EXEMPLARS{pp}), 'nearest');
        end

        if(useWeights && ii > 2)

            for pp=1:length(Recon.ReconObjects{ii}.EXEMPLARS)
                weight_table = zeros( size(Recon.ReconObjects{ii}.EXEMPLARS{pp})-Recon.ReconObjects{ii}.NB_SIZE+1 );

                for kk=1:length(Recon.ReconObjects{ii}.NBWeights{pp})
                    xy = Recon.ReconObjects{ii}.NB_ExemplarLookup{pp}(kk, :);

                    weight_table(xy(1), xy(2)) = Recon.ReconObjects{ii}.NBWeights{pp}(kk);
                end

                nx_Ex_Size = size(Recon.ReconObjects{ii-1}.EXEMPLARS{pp})-Recon.ReconObjects{ii}.NB_SIZE+1;
                weight_table = resize(weight_table, nx_Ex_Size, 'nearest');

                Recon.ReconObjects{ii-1}.NBWeights{pp} = zeros(1, size(Recon.ReconObjects{ii-1}.Exemplar_NBs{pp},1));
                for kk=1:length(Recon.ReconObjects{ii-1}.NBWeights{pp})
                    xy = Recon.ReconObjects{ii-1}.NB_ExemplarLookup{pp}(kk, :);
                    Recon.ReconObjects{ii-1}.NBWeights{pp}(1, kk) = weight_table(xy(1), xy(2));   
                end
            end
        end
        ll = ll + 1;
    end
else
    [S_star, Recon] = SolidOptimization(Recon, max_iterations, S_star, useWeights);
end
end
