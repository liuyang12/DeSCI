function [L_hat,S_hat] = rpca(patchgroup,Omega,sigma)
%RPCA Robust principal component analysis (PCA).
%   See also PATCHESTIMATE, WNNM_IMDENOISE.
[nr,nc] = size(patchgroup);
p = nnz(Omega)/(nr*nc);
lambda = 0.8*1/sqrt(max(nr,nc));
mu = 2*(sqrt(nr)+sqrt(nc))*sigma;
iternum = 30;
tol = 1e-2;
[L_hat,S_hat,iter] = partial_proximal_gradient_rpca2(patchgroup,Omega,lambda,mu,iternum,tol);

end

