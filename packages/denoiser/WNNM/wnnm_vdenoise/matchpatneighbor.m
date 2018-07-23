function [ nonlocalvec, vindvec ] = matchpatneighbor( ...
    patchmat, neighborindarr, ipat, cframe, selfind, neighbornum, patchnum )
%MATCHPATNEIGHBOR Match all the neighbors of a key patch.

patch = patchmat(:,selfind,cframe); % key patch
neighbors = patchmat(:,neighborindarr(1:neighbornum,ipat),:); % all neighbors
[R,C] = ndgrid(1:size(neighbors,2),1:size(neighbors,3));
distance = sum(bsxfun(@minus,neighbors,patch).^2,1); % ell_2 distance
[~,ind] = sort(distance(:)); % sort distance in an ascending order
% indexes of the most similar (shortest distance) para.patchnum of neighbors
Rind = R(ind);
Cind = C(ind);
nonlocalvec = neighborindarr(Rind(1:patchnum),ipat);
vindvec = Cind(1:patchnum);

end

