function RS = SetupRecon(params)

    add_paths;

    RS = params;
    HALF_NB_SIZE = floor(RS.NB_SIZE/2);
    NUM_EXEMPLARS = length(RS.EXEMPLARS);
    RS.Exemplar_NBs = cell(NUM_EXEMPLARS, 1);
    RS.Exemplar_Index = cell(NUM_EXEMPLARS, 1);
    RS.Exemplar_Params = cell(NUM_EXEMPLARS, 1);
    RS.Exemplar_NBLookup = cell(NUM_EXEMPLARS,1);
    RS.NB_ExemplarLookup = cell(NUM_EXEMPLARS,1);

    nbOffsets = MakeNBOffsets(RS.NB_SIZE);
    nbOffsets = nbOffsets(RS.NB_INDICES);
    
    RS.nbOffsets = zeros([length(RS.NB_INDICES), size(nbOffsets{1})]);
    for ii=1:length(RS.NB_INDICES)
        RS.nbOffsets(ii, :) = nbOffsets{ii}(:);
    end

    if(strcmp(RS.ANN_ALGO,'FLANN'))
        build_params.cores = RS.NUM_CORES;
        build_params.algorithm = 'kdtree';
        build_params.trees = 8;
        build_params.checks = 131;
    elseif(strcmp(RS.ANN_ALGO,'PatchTable'))
       build_params = [];
    else
        error('Invalid ANN_ALGO parameter. Must be either FLANN or PatchTable');
    end

    for ii=1:NUM_EXEMPLARS
        fprintf(1, '\tGetting neighborhoods for exemplar %d of %d ...', ii, NUM_EXEMPLARS);

        if(isfield(RS, 'NB_GrainList'))
            [NB_Hoods RS.NB_ExemplarLookup{ii} RS.NB_GrainList{ii}] = GetAllNHoodsMEX(RS.EXEMPLARS{ii}, RS.NB_SIZE);
            RS.NB_NumGrains{ii} = sum(RS.NB_GrainList{ii} > 0, 2);
            RS.NB_GrainList{ii} = uint16(RS.NB_GrainList{ii});
            RS.NB_GrainList{ii} = RS.NB_GrainList{ii}(:, 1:max(RS.NB_NumGrains{ii}));
        else
            [NB_Hoods RS.NB_ExemplarLookup{ii}] = GetAllNHoodsMEX(RS.EXEMPLARS{ii}, RS.NB_SIZE);
        end
        RS.Exemplar_NBs{ii} = NB_Hoods; 

        % If there is a grain boundary map. Then get the neighborhoods for this as well
        if(isfield(RS, 'E_Edges'))
            [RS.Exemplar_NBs_Edges{ii} RS.NB_ExemplarLookup{ii}] = GetAllNHoodsMEX(RS.E_Edges{ii}, RS.NB_SIZE);
        end

        % If we have an angle map. Then the images are indexed and we need to convert
        % the neighborhoods to their actual angles. We will use the 6,7,8 GSH coeffcients
        % because they have indpendent information.
        if(isAngMap)
            tmp = zeros(size(NB_Hoods, 1), 3*size(NB_Hoods, 2));
            for ll=1:size(ANG_MAP, 2)
                idx_start = (ll-1)*size(NB_Hoods, 2)+1;
                tmp(:, idx_start:(idx_start+size(NB_Hoods,2)-1)) = reshape(ANG_MAP(NB_Hoods, ll), size(NB_Hoods));
            end
            RS.Exemplar_NBs{ii} = tmp;
        end

        % Make a map for the original exemplar image that tells where each
        % neighborhood is stored in the the neighborhood table
        RS.Exemplar_NBLookup{ii} = zeros(size(RS.EXEMPLARS{ii}));
        for ll=1:size(RS.NB_ExemplarLookup{ii}, 1)
           RS.Exemplar_NBLookup{ii}(RS.NB_ExemplarLookup{ii}(ll, 1), RS.NB_ExemplarLookup{ii}(ll, 2)) = ll;  
        end

        fprintf(1, 'Done\n');

        % Build nearest neighbor query indices for these points.
        fprintf(1, '\tBuilding ANN index for exemplar %d of %d ...', ii, NUM_EXEMPLARS);
        RS.Exemplar_Params{ii} = build_params;
        tic;
        if(isfield(RS, 'Exemplar_NBs_Edges') && ~isAngMap)
            if(isfield(RS, 'NB_NumGrains'))
                if(strcmp(RS.ANN_ALGO, 'FLANN'))
                    [RS.Exemplar_Index{ii} RS.Exemplar_Params{ii}] = flann_build_index([RS.Exemplar_NBs_Edges{ii} RS.NB_NumGrains{ii}*1000]', build_params);
                elseif(strcmp(RS.ANN_ALGO, 'PatchTable'))
                    
                end
            else
                if(strcmp(RS.ANN_ALGO, 'FLANN'))
                    [RS.Exemplar_Index{ii} RS.Exemplar_Params{ii}] = flann_build_index(RS.Exemplar_NBs_Edges{ii}', build_params);
                elseif(strcmp(RS.ANN_ALGO, 'PatchTable'))
                    
                end
            end
        else
            if(strcmp(RS.ANN_ALGO, 'FLANN'))
                [RS.Exemplar_Index{ii} RS.Exemplar_Params{ii}] = flann_build_index([RS.Exemplar_NBs{ii}]' , build_params);
            elseif(strcmp(RS.ANN_ALGO, 'PatchTable'))
                
            end
        end

        elapsed_time = toc;
        fprintf(1, 'Done. Time to Build = %f seconds\n', elapsed_time);
    end

    % Setup the reconstruction change table. This tells us which
    % neighborhoods have changed in the reconstruction over the last
    % iteration. This allows us to only search for new matching
    % neighborhoods on those that have changed.
    RS.change_table = ones(RS.RECON_SIZE);

    % At each iteration we need to build a table that says for each voxel
    % in the reconstruction where the best matching neighborhoods is in each
    % exemplar.
    RS.NNB_Table = zeros([RS.RECON_SIZE NUM_EXEMPLARS]);

end
