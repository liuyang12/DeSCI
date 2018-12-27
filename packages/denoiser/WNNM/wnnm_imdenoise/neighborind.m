function [neighborindarr,neighbornumarr,selfindarr] = neighborind(dim,opt)
%NEIGHBORIND Pre-calculation of the indexes of the neighbors within the
%search window.
%   [neighborindarr,neighbornumarr,selfindarr]=NEIGHBORIND(dim,opt) returns
%    neighborindarr   - the array of all neighbor patch indexes within
%                       search window for each key patch, which is a
%                       two-dimentional array with each column denotes all
%                       neighbor patch indexes for each key patch,
%    neighbornumarr   - number of neighbors in gridind (grid index) within
%                       search window, which is a row vector and each
%                       column corresponds to that of neighborindarr,
%    selfindarr       - index of the key patch itself, which is a row
%                       vector and each column corresponds to that of
%                       neighborindarr,
%   where the inputs are
%    dim              - dimension of the image (transfer of the image, nor
%                       the image sequence is not necessary and could be 
%                       memory-consuming),
%    opt              - options for that function, and those options are
%                         $option$                 $[default] value$   
%                         opt.patchsize            
%                         opt.windowsize
%                         opt.patchstep
%   See also WNNM_IMDENOISE.
w = opt.windowsize; % size of the search window for neighbors
rownum = dim(1)-opt.patchsize+1; % number of patch rows
colnum = dim(2)-opt.patchsize+1; % number of patch columns
% indexes of all the patches in the image
ind = reshape([1:rownum*colnum],[rownum,colnum]);

% consider the adjacent patches along the patch step, which forms a grid
rowgridind = [1:opt.patchstep:rownum]; % indexes of row grid
  rowgridind = [rowgridind,rowgridind(end)+1:rownum]; % residual steps
colgridind = [1:opt.patchstep:colnum]; % indexes of column grid
  colgridind = [colgridind,colgridind(end)+1:colnum]; % residual steps

rowgridnum = length(rowgridind);
colgridnum = length(colgridind);

% initialization of the output arrays
neighborindarr = zeros([(2*w+1)*(2*w+1),rowgridnum*colgridnum],'uint32');
neighbornumarr = zeros([1,rowgridnum*colgridnum],'uint32');
selfindarr     = zeros([1,rowgridnum*colgridnum],'uint32');

% assignment of the neighbors within search window for each key patch in
% the grid index
for rgrid = 1:rowgridnum
    r = rowgridind(rgrid); % index of row grid
    for cgrid = 1:colgridnum
        c = colgridind(cgrid); % index of column grid
        patchind = (c-1)*rownum+r; % index of key patch
        patchgridind = (cgrid-1)*rowgridnum+rgrid; % grid index of key patch
        
        top    = max(r-w,1); % index of the top of search window
        button = min(r+w,rownum); % button
        left   = max(c-w,1); % left
        right  = min(c+w,colnum); % right
        
        neighborind = ind(top:button,left:right); % indexes of the neighbors of the key patch
        
        neighbornumarr(patchgridind) = numel(neighborind); % number of neighbors in grid index (within search window)
        neighborindarr(1:neighbornumarr(patchgridind),patchgridind) ... 
            = neighborind(:); % each column is the indexes of the neighbors in grid index
        selfindarr(patchgridind) = patchind; % index of key patch itself
    end
end

end

