function X = nnm(Y,lambda)
%NNM Nuclear norm minimization (NNM).
%   X=NNM(Y,lambda) returns the low-rank approximation under nuclear norm 
%   minimization.
%   Problem formulation
%    X_=argmin_X ||Y-X||_F^2+lambda||X||_*.
%   See also PATCHESTIMATE, WNNM_IMDENOISE.
[U,SigmaY,V] = svd(full(Y),'econ'); % singular value decomposition (economic)
X = U*max(SigmaY-lambda/2,0)*V';

end

