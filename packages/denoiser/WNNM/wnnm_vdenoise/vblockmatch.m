function [nonlocalarr,vindarr] = vblockmatch(cframe,patchmat,neighborindarr,neighbornumarr,selfindarr,para)
%VBLOCKMATCH Block matching of patches with non-local similarity for
%video/volume denoising.
%   nonlocalarr=VBLOCKMATCH(cframe,patchmat,neighborindarr,neighbornumarr,
%   selfindarr,para) returns the indexes of the patches with non-local
%   similarity.
%   See slao WNNM_VDENOISE.
gridsize = length(neighbornumarr);

nonlocalarr = zeros([para.patchnum gridsize],'uint32');
vindarr = zeros([para.patchnum gridsize],'uint8');
% todo: adapt the function to parfor for improved performance

% backup of version without parfor
for i = 1:gridsize
    patch = patchmat(:,selfindarr(i),cframe); % key patch
    neighbors = patchmat(:,neighborindarr(1:neighbornumarr(i),i),:); % all neighbors
    [R,C] = ndgrid(1:size(neighbors,2),1:size(neighbors,3));
    distance = sum(bsxfun(@minus,neighbors,patch).^2,1); % ell_2 distance
    [~,ind] = sort(distance(:)); % sort distance in an ascending order
    % indexes of the most similar (shortest distance) para.patchnum of neighbors
    Rind = R(ind);
    Cind = C(ind);
    nonlocalarr(:,i) = neighborindarr(Rind(1:para.patchnum),i);
    vindarr(:,i) = Cind(1:para.patchnum);
end

% % untested version with parfor
% patchnum = para.patchnum;
% parfor (i = 1:gridsize, gridsize)
%     selfind = selfindarr(i);
%     neighbornum = neighbornumarr(i);
%     [nonlocalvec,vindvec] = matchpatneighbor(patchmat,neighborindarr,...
%         i,cframe,selfind,neighbornum,patchnum);
%     nonlocalarr(:,i) = nonlocalvec;
%     vindarr(:,i) = vindvec;
% end

% for i = 1:gridsize
%     patch = patchmat(:,selfindarr(i),cframe); % key patch
%     neighbors = patchmat(:,neighborindarr(1:neighbornumarr(i),i),:); % all neighbors
%     [R,C] = ndgrid(1:size(neighbors,2),1:size(neighbors,3));
%     distance = sum(bsxfun(@minus,neighbors,patch).^2,1); % ell_2 distance
%     [~,ind] = sort(distance(:)); % sort distance in an ascending order
%     % indexes of the most similar (shortest distance) para.patchnum of neighbors
%     Rind = R(ind);
%     Cind = C(ind);
%     nonlocalarr(:,i) = neighborindarr(Rind(1:para.patchnum),i);
%     vindarr(:,i) = Cind(1:para.patchnum);
% end

end

    