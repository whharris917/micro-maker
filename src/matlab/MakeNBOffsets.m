function [nbOffsets] = MakeNBOffsets(NB_SIZE)

    HS = floor(NB_SIZE/2);
    
    [x y z] = ndgrid((0-HS):(0+HS), (0-HS):(0+HS), 0); 
    nbOffsets{1} = [x(:) y(:) z(:)]; % (0,0,1)
    nbOffsets{2} = [z(:) x(:) y(:)]; % (0,1,0)
    nbOffsets{3} = [x(:) z(:) y(:)]; % (1,0,0)

end
