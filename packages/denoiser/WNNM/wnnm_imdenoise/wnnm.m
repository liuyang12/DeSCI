function X = wnnm(Y,c,nsigma,iternum)
%WNNM Weighted nuclear norm minimization (WNNM).
%   X=WNNM(Y,c,nsigma) returns the low-rank approximation under weighted
%   nuclear norm minimization, where Y is the patch group with each column
%   a single patch (vectorized), c is the parameter for the determinzation
%   of the weight, and nsigma is the estimated noise level of the patch.
%   Note that weights here are in a non-descending order w.r.t. the
%   non-ascending order of the singular value, and the WNNM result is the
%   guaranteed convergence of the fix point.
%   See also PATCHESTIMATE, WNNM_IMDENOISE.
[U,SigmaY,V] = svd(full(Y),'econ'); % singular value decomposition (economic)
[~,patchnum] = size(Y);
% % estimation according to WNNM paper in CVPR'14
% SigmaY = diag(SigmaY);
% SigmaX = sqrt(max(bsxfun(@minus,SigmaY.^2,patchnum*nsigma^2),0));
% w = c*sqrt(patchnum)./(SigmaX+eps); % weight vector
% SigmaX = max(SigmaY-w,0); % estimated singular values of X
% X = U*diag(SigmaX)*V';
% closed-form estimation according to WNNM application paper in IJCV'17
SigmaY = diag(SigmaY);
C0 = c*sqrt(patchnum)*2*nsigma^2;
Delta = (SigmaY+eps).^2-4*C0;
effind = Delta>0;
SigmaX = max((SigmaY(effind)-eps+sqrt(Delta(effind)))./2,0);
X = U(:,effind)*diag(SigmaX)*V(:,effind)';

end

