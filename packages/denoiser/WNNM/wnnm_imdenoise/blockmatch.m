function nonlocalarr = blockmatch(patchmat,neighborindarr,neighbornumarr,selfindarr,para)
%BLOCKMATCH Block matching of patches with non-local similarity.
%   nonlocalarr=BLOCKMATCH(patchmat,neighborindarr,neighbornumarr,
%   selfindarr,para) returns the indexes of the patches with non-local
%   similarity.
%   See slao WNNM_IMDENOISE.
gridsize = length(neighbornumarr);

nonlocalarr = zeros([para.patchnum,gridsize],'uint32');
for i = 1:gridsize
    patch = patchmat(:,selfindarr(i)); % key patch
    neighbors = patchmat(:,neighborindarr(1:neighbornumarr(i),i)); % all neighbors
    distance = sum(bsxfun(@minus,neighbors,patch).^2); % ell_2 distance
    [~,ascendind] = sort(distance); % sort distance in an ascending order
    % indexes of the most similar (shortest distance) para.patchnum of neighbors
    nonlocalarr(:,i) = neighborindarr(ascendind(1:para.patchnum),i);
end

end

    