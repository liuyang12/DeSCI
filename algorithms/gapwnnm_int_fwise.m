function [v,psnrall] = gapwnnm_int_fwise( y, opt )
%GAPWNNM_INT_FWISE Integrated GAP and WNNM for snapshot compressive 
%imaging.
%   Note integrated version of GAP-WNNM instead of embedded for
%   acceleration of WNNM with less block-match process for multiple
%   iterations.
%   See also GAPDENOISE.
if nargin<2
    opt = [];
end
% [0] default parameter configuration, to be specified
A  = @(x) M_func(x);
At = @(z) Mt_func(z);
if isfield(opt,'Mfunc'),  A  = @(x) opt.Mfunc(x);  end
if isfield(opt,'Mtfunc'), At = @(z) opt.Mtfunc(z); end
% Phisum   = ;
% v0       = At(y); % start point (initialization of iteration)
lambda   = 1;     % correction coefficiency
maxiter  = 300;     % maximum number of iteration
acc      = 1;       % enable acceleration
sigma    = 10/255;  % noise deviation 
flag_iqa = true;    % flag of showing image quality assessments
save_iter_image = false; % flag of saving all the images of the iteration 
                         %  process (true only necessary)
                         %  require assigned directory opt.iter_image_dir

if isfield(opt,'Phisum'),     Phisum = opt.Phisum;   end
if isfield(opt,'v0'),             v0 = opt.v0;       end
if isfield(opt,'lambda'),   lambda   = opt.lambda;   end
if isfield(opt,'maxiter'),   maxiter = opt.maxiter;  end
if isfield(opt,'acc'),           acc = opt.acc;      end
if isfield(opt,'sigma'),       sigma = opt.sigma;    end
if isfield(opt,'flag_iqa'), flag_iqa = opt.flag_iqa; end
if isfield(opt,'save_iter_image'), save_iter_image = opt.save_iter_image; end

if ~exist('v0','var') || isempty(v0)
    v0 = At(y); % start point (initialization of iteration)
end
y1 = zeros(size(y),'like',y);
psnrall = []; % return empty with no ground truth
% ssimall = []; % return empty with no ground truth
% [1] start iteration
v = v0; % initialization
% opt.sigma = sigma;
k = 1; % current number of iteration
nonlocalarr = uint32.empty;
vindarr     = uint8.empty;
for isig = 1:length(sigma)
    opt.sigma = sigma(isig);
    for iter = 1:maxiter(isig)
        for ii = 1:1
            % [1.1] Euclidean projection
            yb = A(v);
            if acc % enable acceleration
                y1 = y1+(y-yb);
                v = v+lambda*(At((y1-yb)./Phisum)); % v=v+lambda*(At*A)^-1*At*dy
            else
                v = v+lambda*(At((y-yb)./Phisum));
            end
        end
        % [1.2] Denoising to match the video prior
        % WNNM video denoising (MATLAB-style matrix version)
        % v = wnnm_vdenoise(v,[],opt); % opt.sigma
        noisyv = v;
        para = vdefparaconf(opt); % default parameter configuration
        para.lambda = para.vlambda; % resolve the conflict of the parameter name

        [nrow,ncol,nframe] = size(noisyv); % grayscale video (color todo)
        % [1] pre-calculation of the indexes of the neighbors within the search
        %     window
        [neighborindarr,neighbornumarr,selfindarr] = neighborind([nrow ncol],para);

        % [2] WNNM denoisng for several iterations
        estv = noisyv;
        for iit = 1:para.iternum
            % correction between adjacent iterations
            if ~para.adaptboost % 
                estv = estv + para.delta*(noisyv-estv); 
            else % adaptive boosting for WNNM-based denoising
                Eestv = mean(abs(estv(:)).^2);
                Enoise  = para.abeta*abs(para.nsigma^2-var(noisyv(:)-estv(:)));
                rho = sqrt(Eestv)/(sqrt(Eestv)+sqrt(max(Eestv-Enoise,0)));
                fprintf('    Iteration % 2d, rho = %.3f.\n',iit,rho);
                estv = estv + (1-rho)*(noisyv-estv); 
            end
            
            blockmatch_period = para.blockmatch_period;
            % inner loop to reuse the block-matching results
                if mod(k,blockmatch_period) == 1 % block matching periodically
                    nonlocalarr = uint32.empty;
                    parfor (iframe = 1:nframe,nframe) % maximum nframe parpool for optimal
                        estim = estv(:,:,iframe);
                        noisyim = noisyv(:,:,iframe);
                        % [2.1] video to patches
                        [rawpatchmat,~] = im2patch(estim,noisyim,para);
                        % [2.2] calculate the patches with non-local similarity for each 
                        %       key patch
                        curnonlocalarr = blockmatch(rawpatchmat,...
                            neighborindarr,neighbornumarr,selfindarr,para);
                        nonlocalarr(:,:,iframe) = curnonlocalarr;
                    end
                end
            % frame-wise denosing using parfor instead of for to get accelerated
            % performance, requiring Parallel Computing Toolbox.
            
            parfor (iframe = 1:nframe,nframe) % maximum nframe parpool for optimal
            % for iframe = 1:nframe % maximum nframe parpool for optimal
                estim = estv(:,:,iframe);
                noisyim = noisyv(:,:,iframe);
                % [2.1] video to patches
                [rawpatchmat,nsigmamat] = im2patch(estim,noisyim,para);
                if iit==1 % initial noise level of each patch
                    nsigmamat = para.nsigma*ones(size(nsigmamat));
                end
                curnonlocalarr = nonlocalarr(:,:,iframe);
                % [2.3] patch estimation by means of WNNM
                [estpatchmat,frqpatchmat] = patchestimate(curnonlocalarr,...
                    rawpatchmat,nsigmamat,selfindarr,para);
                % [2.4] aggregate overlapped patches to the whole image
                curim = patch2im(estpatchmat,frqpatchmat,size(noisyim),...
                    para.patchsize);
                
                estv(:,:,iframe) = curim;
            end

            if mod(iit-1,para.innerloop)==0 
                % % [2.2] calculate the patches with non-local similarity for 
                % %       each key patch
                % [nonlocalarr,vindarr] = vblockmatch(cframe,rawpatchmat,...
                %     neighborindarr,neighbornumarr,selfindarr,para);
                % less non-local patches with lower noise level
                para.patchnum = para.patchnum-10; 
            end
            % [1.4] save and show intermediate results of psnr and ssim
            if flag_iqa && isfield(opt,'orig') && (~isempty(opt.orig))
                psnrall(k) = psnr(double(estv),double(opt.orig)); % record all psnr
                % ssimall(k) = ssim(double(estv),opt.orig); % record all ssim
                if (mod(k,5)==0) 
                    fprintf('  GAP-%s iteration % 4d, sigma %.1f, PSNR %2.2f dB.\n',...
                        upper(opt.denoiser),k,opt.sigma*255,psnrall(k));
                end
                if save_iter_image % save all the images of the iteration process
                    MAXB = opt.MAXB;
                    % save each frame in the assigned directory
                    nim = size(estv,ndims(estv)); % number of frames in the video
                    for iim = 1:nim
                        im = estv(:,:,iim);
                        impsnr = psnr(im,opt.orig(:,:,iim));
                        imssim = ssim(im,opt.orig(:,:,iim));
                        if k == 1 % mkdir for the first iteration
                            subdir = sprintf('%s/frame%02d',opt.iter_image_dir,iim); % frame-wise
                            if ~exist(subdir,'dir')
                                mkdir(subdir);
                            end
                        end
                        imwrite(uint8(im*MAXB),sprintf('%s/frame%02d/frame%02d_iter%03d_sigma%.1f_psnr%2.2f_ssim%.4f.png',...
                            opt.iter_image_dir,iim,iim,k,opt.sigma*MAXB,impsnr,imssim));
                    end
                end
            end
            k = k+1;
        end % WNNM loop [iternum]
        v = estv;
        % % [1.3] update noise standard deviation
        % nsigma = eta*sqrt(abs(sigma^2-var(v(:)-v0(:))));
        % opt.sigma = nsigma;
        
    end % GAP loop [maxiter]
end % sigma loop [length(sigma)]

end            
            
            
