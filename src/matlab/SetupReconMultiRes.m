function ReconH = SetupReconMultiRes(params)

ReconH = params;
ReconH.ReconObjects = cell(params.NUM_LEVELS, 1);
ReconH.NUM_EXEMPLARS = length(ReconH.EXEMPLARS);

if(iscell(ReconH.EXEMPLARS{1}))
    P = ReconH.EXEMPLARS;
    NUM_LEVELS = length(P{1});
    ReconH.ReconObjects = cell(NUM_LEVELS, 1);
else
    P = cell(ReconH.NUM_EXEMPLARS, 1);
    for ii=1:ReconH.NUM_EXEMPLARS
        P{ii} = BuildExemplarPyramid(ReconH.EXEMPLARS{ii});
    end
end

RECON_SIZE = params.FULL_RECON_SIZE;

for level=1:ReconH.NUM_LEVELS
    fprintf(1, 'Level %d\n', level);
    
    Exyz = cell(ReconH.NUM_EXEMPLARS,1);
    for ex=1:ReconH.NUM_EXEMPLARS
        Exyz{ex} = P{ex}{level};
    end
    
    if(~isfield(params, 'NUM_CORES'))
        rparams.NUM_CORES = 1;
    else
        rparams.NUM_CORES = params.NUM_CORES;
    end
    
    rparams.NB_SIZE = params.NB_SIZES(level);
    rparams.NB_INDICES = params.NB_INDICES;
    rparams.RECON_SIZE = RECON_SIZE;
    rparams.EXEMPLARS = Exyz;
    
    if(isfield(params, 'ANN_ALGO'))
        rparams.ANN_ALGO = params.ANN_ALGO;
    end
    
    if(isfield(params, 'UseNeighborhoodSearchWeights'))
        rparams.UseNeighborhoodSearchWeights = params.UseNeighborhoodSearchWeights;
    end
    
    ReconH.ReconObjects{level} = SetupRecon(rparams);
    RECON_SIZE = ceil(RECON_SIZE / 2);
end
