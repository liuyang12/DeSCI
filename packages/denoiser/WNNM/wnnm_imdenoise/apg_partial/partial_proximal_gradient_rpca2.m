function [X_k_A,X_k_E,numIter] = partial_proximal_gradient_rpca2(D, omega, lambda,  mu, maxIter, tol)
%the following code aims to solve
%       \min_{L,S} mu*\|A\|_*+lambda*mu*\|E\|_1+0.5*\|P_omega(D-A-E)\|_F^2
%       ---- where omega is 1 for existing entry and 0 for missing entry
%modified by Yuhong Xu, Temasek Laboratory, National University of Singapore
%while the original code comes from ----- 

% D - m x n matrix of observations/data (required input)
% lambda - weight on sparse error term in the cost function (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maximum number of iterations
%         - DEFAULT 10000, if omitted or -1.
% [A_hat, E_hat] - estimates for the low-rank part and error part, respectively
% numIter - number of iterations until convergence

% References
% "Robust PCA: Exact Recovery of Corrupted Low-Rank Matrices via Convex Optimization", J. Wright et al., preprint 2009.
% "An Accelerated Proximal Gradient Algorithm for Nuclear Norm Regularized Least Squares problems", K.-C. Toh and S. Yun, preprint 2009.
%
% 
% Arvind Ganesh, Minming Chen, Summer 2009. Questions? abalasu2@illinois.edu, v-minmch@microsoft.com.
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

%addpath PROPACK;

DISPLAY_EVERY = 20 ;
maxLineSearchIter = 200 ; % maximum number of iterations in line search per outer iteration

%% Initializing optimization variables

[m,n] = size(D) ;

t_k = 1 ; % t^k
t_km1 = 1 ; % t^{k-1}

tau_0 = 2 ; % square of Lipschitz constant for the RPCA problem

X_km1_A = zeros(m,n) ; X_km1_E = zeros(m,n) ; % X^{k-1} = (A^{k-1},E^{k-1})
X_k_A = zeros(m,n) ; X_k_E = zeros(m,n) ; % X^{k} = (A^{k},E^{k})

mu_k = mu ;
tau_k = tau_0;

converged = 0 ;
numIter = 0 ;
% addpath PROPACK;
sv = 5;
%% Start main loop
while ~converged
    
    Y_k_A = X_k_A + ((t_km1 - 1)/t_k)*(X_k_A-X_km1_A) ;
    Y_k_E = X_k_E + ((t_km1 - 1)/t_k)*(X_k_E-X_km1_E) ;
    
        %projection onto the set of Omega, i.e., setting other entries zeros
        X_km1_A = Y_k_A - (1/tau_k)*(Y_k_A+Y_k_E-D).*omega ;
        X_km1_E = Y_k_E - (1/tau_k)*(Y_k_A+Y_k_E-D).*omega ;
        
        if choosvd(n, sv) == 1
            [U S V] = lansvd(X_km1_A, sv, 'L');
        else
            [U S V] = svd(X_km1_A, 'econ');
        end
        diagS = diag(S);
        svp = length(find(diagS > mu_k/tau_k));
        if svp < sv
            sv = min(svp + 1, n);
        else
            sv = min(svp + round(0.05*n), n);
        end
        
        X_kp1_A = U(:, 1:svp) * diag(diagS(1:svp) - mu_k/tau_k) * V(:, 1:svp)';    
        X_kp1_E = sign(X_km1_E) .* pos( abs(X_km1_E) - lambda*mu_k/tau_k );
        
        rankA  = sum(diagS>mu_k/tau_k);
        cardE = sum(sum(double(abs(X_kp1_E)>0)));
        
    
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    %projection onto the set of Omega, checking stopping criteria
    X_km1_A = (X_kp1_A + X_kp1_E - Y_k_A - Y_k_E).*omega ;
    Y_k_A = tau_k*(Y_k_A-X_kp1_A) + X_km1_A ;
    Y_k_E = tau_k*(Y_k_E-X_kp1_E) + X_km1_A ;
    
    s1 = sqrt(Y_k_A(:)'*Y_k_A(:)+Y_k_E(:)'*Y_k_E(:));
    s2 = sqrt(X_kp1_A(:)'*X_kp1_A(:)+X_kp1_E(:)'*X_kp1_E(:));
    stoppingCriterion = s1 / (tau_k*max(1,s2));
    if stoppingCriterion <= tol && numIter > 5
        converged = 1 ;
    end    
   
    t_km1 = t_k ;
    t_k = t_kp1 ;
    X_km1_A = X_k_A ; X_km1_E = X_k_E ;
    X_k_A = X_kp1_A ; X_k_E = X_kp1_E ;
    
    numIter = numIter + 1 ;
    
    if mod(numIter,DISPLAY_EVERY) == 0
        disp(['Iteration ' num2str(numIter) '  rank(A) ' num2str(rankA) ...
            ' ||E||_0 ' num2str(cardE) '  Stopping Criterion ' ...
            num2str(stoppingCriterion)]) ;
%         mu_k
    end
    
    if nargin > 8
        fprintf(fid, '%s\n', ['Iteration ' num2str(numIter) '  rank(A)  ' num2str(rankA) ...
            '  ||E||_0  ' num2str(cardE) '  Stopping Criterion   ' ...
            num2str(stoppingCriterion)]) ;
    end
        
    
    if ~converged && numIter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
    
end
