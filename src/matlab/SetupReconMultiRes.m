function ReconH = SetupReconMultiRes(params)

ReconH = params;
ReconH.ReconObjects = cell(params.NUM_LEVELS, 1);
ReconH.NUM_EXEMPLARS = length(ReconH.EXEMPLARS);
P = ReconH.EXEMPLARS;
NUM_LEVELS = length(P{1});
ReconH.ReconObjects = cell(NUM_LEVELS, 1);
RECON_SIZE = params.FULL_RECON_SIZE;

for level=1:ReconH.NUM_LEVELS

    fprintf(1, 'Level %d\n', level);
    
    Exyz = cell(ReconH.NUM_EXEMPLARS,1);
    for ex=1:ReconH.NUM_EXEMPLARS
        Exyz{ex} = P{ex}{level};
    end
    
    rparams.NUM_CORES = params.NUM_CORES;
    rparams.NB_SIZE = params.NB_SIZES(level);
    rparams.NB_INDICES = params.NB_INDICES;
    rparams.RECON_SIZE = RECON_SIZE;
    rparams.EXEMPLARS = Exyz;
    rparams.ANN_ALGO = params.ANN_ALGO;
    rparams.UseNeighborhoodSearchWeights = params.UseNeighborhoodSearchWeights;
    
    ReconH.ReconObjects{level} = SetupRecon(rparams);
    
    RECON_SIZE = ceil(RECON_SIZE / 2);
    
end
