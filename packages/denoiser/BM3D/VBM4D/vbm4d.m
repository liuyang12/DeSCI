function [ y_est ] = vbm4d( z, sigma, profile, do_wiener, sharpen, deflicker, verbose )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  V-BM4D is an algorithm for attenuation of additive white Gaussian noise  
%  in videos. This algorithm reproduces the results from the articles:
%
%  [1] M. Maggioni, G. Boracchi, A. Foi, K. Egiazarian, "Video Denoising 
%      Using Separable 4D Nonlocal Spatiotemporal Transforms", 
%      Proc. SPIE Electronic Imaging 2011, San Francisco, CA, USA.
%
%  [2] M. Maggioni, G. Boracchi, A. Foi, K. Egiazarian, "Video Denoising, 
%      Deblocking and Enhancement Through Separable 4-D Nonlocal 
%      Spatiotemporal Transforms", IEEE Trans. on Image Proc., 
%      Vol. 21, No. 9, Sep. 2012. doi:10.1109/TIP.2012.2199324
%
%
%  FUNCTION INTERFACE:
% 
%  INPUTS:
%     1) z         (3D or 4D array) : noisy video
%     2) sigma             (double) : noise standard deviation, a value
%                                     equal to -1 enables noise estimation
%                                   : (default is -1)
%     3) profile             (char) : 'lc' --> low complexity profile 
%                                   : 'np' --> normal profile 
%                                   : 'mp' --> modified profile (high complexity)
%                                     (default is 'np')
%     4) do_wiener        (logical) : perform collaborative Wiener filtering
%                                     (default is 1)
%     5) sharpen           (double) : sharpening parameter (should be >0):
%                                      if < 1 --> smoothing
%                                      if = 1 --> sharpening disabled
%                                      if > 1 --> sharpening
%                                     (default is 1)
%     6) deflicker         (double) : deflickering parameter (should be <0):
%                                      if < 1 --> deflickering
%                                      if = 1 --> deflickering disabled
%                                      if > 1 --> sharpening
%                                     (default is 1)
%     7) verbose          (logical) : 0 --> do not print output information
%                                     1 --> print information to screen
%                                     (default is 0)
%
%   Only input z is required, all other inputs are optional. Optional
%   inputs can be omitted, and assume their default value when set to 'empty' [] .
%   Input z and sigma must be given with respect to the same intensity
%   range. Any intensity range is valid.
%
%  OUTPUTS:
%     1) y_est      (3D or 4D array) : Final estimate
% 
%  TYPICAL USAGE EXAMPLES:
%
%  Case: known noise standard deviation and default parameters
%   y = read_video('path/to/video.avi');
%   sigma = 25;
%   z = y + sigma*randn(size(y));
%   y_est = vbm4d(z, sigma);
%
%  Case: unknown noise standard deviation and default parameters
%   y_est = vbm4d(z);
%
%  Case: unknown noise standard deviation and low complexity parameters
%   y_est = vbm4d(z, [], 'lc');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2010-2014 Tampere University of Technology.
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% AUTHOR:
%     Matteo Maggioni, email: matteo.maggioni _at_ tut.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Arguments
if ~exist('z','var')
    error('Input video "z" is required.')
end
if ~exist('sigma','var') || isempty(sigma)
    sigma = -1;
end
if ~exist('profile','var') || isempty(profile)
    profile = 'np';
end
if ~exist('do_wiener','var') || isempty(do_wiener)
    do_wiener = 1; 
end
if ~exist('sharpen','var') || isempty(do_wiener)
    sharpen = 1; 
end
if ~exist('deflicker','var') || isempty(do_wiener)
    deflicker = 1; 
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = 1;
end
data_type = 'single';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input check
%%%%
if ~strcmpi(profile,'lc') && ~strcmpi(profile,'np') && ~strcmpi(profile,'mp')
    warning(['Invalid profile argument "',profile,'". Assuming normal profile.']);
    profile = 'np';
end
if ~exist('z','var') || ndims(z)>4
    error('Invalid input image dimension. Only 3-D grayscale or 4-D color videos are accepted.');
end
z      = cast(z, data_type);
is_rgb = ndims(z)==4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Normalization to [0,1] range
%%%%
maxz = max(z(:));
minz = min(z(:));
scale = 0.7;
shift = (1-scale)/2;
z = (z-minz)/(maxz-minz);
z = z*scale+shift;
if sigma~=-1
    sigma = sigma/(maxz-minz)*scale;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Algorithm parameters
%%%%

% Transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
transform_2D_HT_name = 'bior1.5';    % 2-D spatial transform (Hard thresholding)
transform_3rd_dim_HT_name = 'dct';   % 1-D temporal transform
transform_4th_dim_HT_name = 'haar';  % 1-D nonlocal transform
transform_2D_Wie_name = 'dct';       % 2-D spatial transform (Wiener filtering)
transform_3rd_dim_Wie_name = 'dct';  % 1-D temporal transform
transform_4th_dim_Wie_name = 'haar'; % 1-D nonlocal transform

% Motion estimation parameters:
motion_est_type = 1;    % motion estimation strategy
                        %  0 --> full search
                        %  1 --> fast search

% Hard-thresholding (HT) parameters:
N = 8;                  % block size
h_minus = 4;            % backward temporal extent
h_plus = 4;             % forward temporal extent
Nstep = 4;              % step between each processed volume
Nme = 11;               % motion estimation search window
tau_traj = -1;          % similarity threshold for block matching in motion estimation
                        % if -1, then tau_traj will be a function of sigma
M = 16;                 % maximum size of the 4-D group
Nnl = 13;               % nonlocal volume search window
tau_match = 0.5;        % similarity threshold for volume matching in nonlocal grouping
lambda_thr4D = 2.7;     % hard-threshold lambda parameter
alphaDC = sharpen;      % sharpening parameter for 4-D DC hyper-plane
alphaAC = deflicker;    % sharpening parameter for 4-D AC hyper-plane
beta = 2;               % Kaiser window parameter

% Wiener filtering parameters:
N_wiener = 8;           % block size
h_minus_wiener = 4;     % backward temporal extent
h_plus_wiener = 4;      % forward temporal extent
Nstep_wiener = 4;       % step between each processed volume
Nme_wiener = 11;        % motion estimation search window
tau_traj_wiener = 0.05; % similarity threshold for block matching in motion estimation
M_wiener = 16;          % maximum size of the 4-D group
Nnl_wiener = 11;        % nonlocal volume search window
tau_match_wiener = 0.5; % similarity threshold for volume matching in nonlocal grouping
beta_wiener = 1.3;      % Kaiser window parameter

if strcmpi(profile,'lc')
    N = 8;
    h_minus = 3;
    h_plus = 3;
    Nstep = 6;
    M = 1;
    N_wiener = 8;
    h_minus_wiener = h_minus;
    h_plus_wiener = h_plus;
    Nstep_wiener = Nstep;
    M_wiener = M;
    motion_est_type = 1;
elseif strcmpi(profile,'mp')
    h_minus = 7;
    h_plus = 7;
    Nstep = 4;
    M = 32;
    Nnl = 15;
    h_minus_wiener = h_minus;
    h_plus_wiener = h_plus;
    Nstep_wiener = Nstep;
    M_wiener = N;
    Nnl_wiener = Nnl;
    motion_est_type = 0;
end

% Transform matrices
H = h_plus + h_minus + 1;
[Tfor, Tinv]  = get_transform_matrix( N, transform_2D_HT_name, 0 );
Tfor = cast(Tfor, data_type);
Tinv = cast(Tinv, data_type);
Tfor3 = cell(1,H);
Tinv3 = cell(1,H);
for i=1:H
    [Tfor3{i}, Tinv3{i}] = get_transform_matrix(i, transform_3rd_dim_HT_name);
    Tfor3{i} = cast(Tfor3{i}, data_type);
    Tinv3{i} = cast(Tinv3{i}, data_type);
end
Tfor4 = cell(1,M);
Tinv4 = cell(1,M);
for i=2.^(0:log2(M))
    [Tfor4{i}, Tinv4{i}] = get_transform_matrix(i, transform_4th_dim_HT_name, 0);
    Tfor4{i} = cast(Tfor4{i}, data_type);
    Tinv4{i} = cast(Tinv4{i}, data_type);
end

H_wiener = h_plus_wiener + h_minus_wiener + 1;
[TforW, TinvW] = get_transform_matrix( N_wiener, transform_2D_Wie_name, 0 );
TforW = cast(TforW, data_type);
TinvW = cast(TinvW, data_type);
Tfor3W = cell(1,H_wiener);
Tinv3W = cell(1,H_wiener);
for i=1:H_wiener
    [Tfor3W{i}, Tinv3W{i}] = get_transform_matrix(i, transform_3rd_dim_Wie_name);
    Tfor3W{i} = cast(Tfor3W{i}, data_type);
    Tinv3W{i} = cast(Tinv3W{i}, data_type);
end
Tfor4W = cell(1,M_wiener);
Tinv4W = cell(1,M_wiener);
for i=2.^(0:log2(M_wiener))
    [Tfor4W{i}, Tinv4W{i}] = get_transform_matrix(i, transform_4th_dim_Wie_name, 0);
    Tfor4W{i} = cast(Tfor4W{i}, data_type);
    Tinv4W{i} = cast(Tinv4W{i}, data_type);
end

% Kaiser windows
Kwin = cast(kaiser(N,beta) * kaiser(N,beta)', data_type);
Kwin_wiener = cast(kaiser(N_wiener,beta_wiener) * kaiser(N_wiener,beta_wiener)', data_type);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Filtering
%%%%

if verbose
    fprintf('Filtering video of resolution %dx%dpx and %d frames \n', size(z,1), size(z,2), size(z,3+is_rgb));
    if sigma~=-1
        fprintf('\tGaussian noise standard deviation %.2f \n', sigma*(maxz-minz)/scale);
    else
        fprintf('\tNoise estimation enabled \n');
    end
    fprintf('\tIntensity range of noisy data is [%.2f, %.2f] \n', minz, maxz);
    if is_rgb
        fprintf('\tColor filtering \n');
    else
        fprintf('\tGrayscale filtering \n');
    end
    if do_wiener && alphaDC==1 && alphaAC==1
        fprintf('\tWiener filtering enabled \n');
    else
        fprintf('\tWiener filtering disabled \n');
    end
    if alphaDC~=1
        fprintf('\tSharpening enabled \n');
    end
    if alphaAC~=1
        fprintf('\tDeflickering enabled \n');
    end
    fprintf(['\tParameter profile "',profile,'" \n']);
end

if is_rgb
    [ z, l2 ] = rgb2yuv(z);
    if sigma~=-1
        sigma = transform_sigma( sigma, l2 );
    end
end

% Hard thresholding
basic = tic;
y_hat = vbm4d_thr_mex(z, sigma, N, h_minus, h_plus, Nstep, Nme, tau_traj, ...
    M, Nnl, tau_match, lambda_thr4D, alphaDC, alphaAC, ...
    Tfor, Tinv, Tfor3, Tinv3, Tfor4, Tinv4, ...
    transform_2D_HT_name, transform_3rd_dim_HT_name, transform_4th_dim_HT_name, ...
    Kwin, motion_est_type);
if verbose
    fprintf('Basic estimate completed (%.1fs) \n', toc(basic));
end

% Wiener filtering
if ~do_wiener || alphaDC~=1 || alphaAC~=1
    y_est = y_hat;
else
	wiener = tic;
    y_est = vbm4d_wie_mex(z, y_hat, sigma, N_wiener, h_minus_wiener, h_plus_wiener, ...
        Nstep_wiener, Nme_wiener, tau_traj_wiener, ...
        M_wiener, Nnl_wiener, tau_match_wiener, ...
        TforW, TinvW, Tfor3W, Tinv3W, Tfor4W, Tinv4W, ...
        transform_2D_Wie_name, transform_3rd_dim_HT_name, transform_4th_dim_HT_name, ...
        Kwin_wiener, motion_est_type);
    if verbose
        fprintf('Final estimate completed (%.1fs) \n', toc(wiener))
    end
end

if is_rgb
    y_est = yuv2rgb(y_est);
end

if verbose
    fprintf('Total execution time: %.1f sec \n', toc(basic));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input de-normalization
%%%%
y_est = (y_est-shift)/scale;
y_est = y_est*(maxz-minz)+minz;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          AUXILIARY FUNCTIONS                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ TM, ITM ] = get_transform_matrix( N, type, Nden )
    % Create forward and inverse transform matrices, which allow for 
    % perfect reconstruction. The forward transform matrix is normalized so
    % that the l2-norm of each basis element is 1.
    %
    % N     size of the transform (for wavelets, must be 2^K)
    % type  type of the transform, 'dct', 'dst', 'hadamard', or anything 
    %       that is listed by 'help wfilters' (bi-orthogonal wavelets)
    %       'DCrand' -- an orthonormal transform with a DC and all the 
    %       other basis elements of random nature
    % Nden  If a wavelet transform is generated, this is the desired 
    %       decomposition level. Must be in the range [0, log2(N)-1], where
    %       '0' implies full decomposition.
    %
    % TM    NxN forward transform matrix
    % ITM   NxN inverse transform matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist('Nden','var')
        Nden = 0;
    end
    
    if N==1,
        TM = 1;
    elseif strcmp(type, 'eye')
        TM = eye(N);
    elseif strcmp(type, 'dct')
        TM = dct(eye(N));
    elseif strcmp(type, 'dst')
        TM = dst(eye(N));
    elseif strcmp(type, 'DCrand')
        x = randn(N); x(1:end,1) = 1; [Q,R] = qr(x);
        if (Q(1) < 0),
            Q = -Q;
        end;
        TM = Q';
    elseif strcmp(type, 'hadamard')
        TM = hadamard(N);
    else
        % wavelet transform
        % set periodic boundary conditions, to preserve bi-orthogonality
        dwtmode('per','nodisp');

        [LO_D, HI_D, LO_R, HI_R] = wfilters(type);
        for i = 1:N
              % construct transform matrix
            TM(i,:) = waverec(circshift([1 zeros(1,N-1)],[Nden i-1]), 2.^[Nden Nden:log2(N)], LO_D, -HI_D);
        end
    end

    % normalize the basis elements
    TM = single((TM' * diag(sqrt( 1./sum(TM.^2,2) )) )');
    % compute the inverse transform matrix
    ITM = inv(TM);
end



function [ yuv, l2 ] = rgb2yuv( rgb, type )
	if  ~exist('type','var') || strcmp(type,'opp')
        A=[1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
	elseif strcmp(type,'yuv')
        A = [0.29900000000000   0.58700000000000   0.11400000000000; 
            -0.16873660714285  -0.33126339285715   0.50000000000000; 
             0.50000000000000  -0.41868750000000  -0.08131250000000];
	elseif strcmp(type,'dct')
        A=[1/3 1/3 1/3;1/sqrt(6) 0 -1/sqrt(6);1/sqrt(18) -sqrt(2)/3 1/sqrt(18)];
	end
    yuv = zeros(size(rgb), 'single');
    for T=1:size(rgb,4)
        yuv(:,:,1,T) = rgb(:,:,1,T)*A(1,1) + rgb(:,:,2,T)*A(1,2) + rgb(:,:,3,T)*A(1,3) + ...
                (1-sum(A(1,:).*(A(1,:)>0),2)-sum(A(1,:).*(A(1,:)<0),2))/2;
        yuv(:,:,2,T) = rgb(:,:,1,T)*A(2,1) + rgb(:,:,2,T)*A(2,2) + rgb(:,:,3,T)*A(2,3) + ...
                (1-sum(A(2,:).*(A(2,:)>0),2)-sum(A(2,:).*(A(2,:)<0),2))/2;
        yuv(:,:,3,T) = rgb(:,:,1,T)*A(3,1) + rgb(:,:,2,T)*A(3,2) + rgb(:,:,3,T)*A(3,3) + ...
                (1-sum(A(3,:).*(A(3,:)>0),2)-sum(A(3,:).*(A(3,:)<0),2))/2;
    end
    l2 = sqrt(sum(A.^2,2));
end
function [ rgb ] = yuv2rgb( yuv, type )
    if ~exist('type','var') || strcmp(type,'opp')
        A = [1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
        B = [1   1   2/3; 1    0  -4/3; 1     -1    2/3];
    elseif strcmp(type,'yuv')
        A = [0.29900000000000   0.58700000000000   0.11400000000000; 
            -0.16873660714285  -0.33126339285715   0.50000000000000;
             0.50000000000000  -0.41868750000000  -0.08131250000000];
        B = inv(A);
    elseif strcmp(type,'dct')
        A=[1/3 1/3 1/3;1/sqrt(6) 0 -1/sqrt(6);1/sqrt(18) -sqrt(2)/3 1/sqrt(18)];
        B=[1 sqrt(3/2) 1/sqrt(2);1 0 -sqrt(2);1 -sqrt(3/2) 1/sqrt(2)];
    end
    rgb = zeros(size(yuv), 'single');
    for T=1:size(rgb,4)
        yuv1 = yuv(:,:,1,T) - (1-sum(A(1,:).*(A(1,:)>0),2)-sum(A(1,:).*(A(1,:)<0),2))/2;
        yuv2 = yuv(:,:,2,T) - (1-sum(A(2,:).*(A(2,:)>0),2)-sum(A(2,:).*(A(2,:)<0),2))/2;
        yuv3 = yuv(:,:,3,T) - (1-sum(A(3,:).*(A(3,:)>0),2)-sum(A(3,:).*(A(3,:)<0),2))/2;
        
        rgb(:,:,1,T) = yuv1*B(1,1) + yuv2*B(1,2) + yuv3*B(1,3);
        rgb(:,:,2,T) = yuv1*B(2,1) + yuv2*B(2,2) + yuv3*B(2,3);
        rgb(:,:,3,T) = yuv1*B(3,1) + yuv2*B(3,2) + yuv3*B(3,3);
    end
end
function [ sigma_yuv ] = transform_sigma( sigma, l2 )
    sigma_yuv    = zeros(1,3);
    sigma_yuv(1) = sigma * l2(1);
    sigma_yuv(2) = sigma_yuv(1) * l2(2) / l2(1);
    sigma_yuv(3) = sigma_yuv(1) * l2(3) / l2(1);
end

