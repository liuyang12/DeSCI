% -------------------------------------------------------------------------
%
%     Demo software for V-BM4D volumetric data denoising 
%               Release ver. 1.0  (11 December 2014)
%
% -------------------------------------------------------------------------
%
% The software implements the V-BM4D denoising algorithm described in:
%
%    M. Maggioni, G. Boracchi, A. Foi, K. Egiazarian, "Video Denoising 
%      Using Separable 4D Nonlocal Spatiotemporal Transforms", 
%      Proc. SPIE Electronic Imaging 2011, San Francisco, CA, USA.
%
%    M. Maggioni, G. Boracchi, A. Foi, K. Egiazarian, "Video Denoising, 
%      Deblocking and Enhancement Through Separable 4-D Nonlocal 
%      Spatiotemporal Transforms", IEEE Trans. on Image Proc., 
%      Vol. 21, No. 9, Sep. 2012. doi:10.1109/TIP.2012.2199324
%
% -------------------------------------------------------------------------
%
% authors:               Matteo Maggioni
%                        Alessandro Foi
%
% web page:              http://www.cs.tut.fi/~foi/GCF-BM3D
%
% contact:               firstname.lastname@tut.fi
%
% -------------------------------------------------------------------------
% Copyright (c) 2010-2014 Tampere University of Technology.
% All rights reserved.
% This work should be used for nonprofit purposes only.
% -------------------------------------------------------------------------
%
% Disclaimer
% ----------
%
% Any unauthorized use of these routines for industrial or profit-oriented 
% activities is expressively prohibited. By downloading and/or using any of 
% these files, you implicitly agree to all the terms of the TUT limited 
% license (included in the file Legal_Notice.txt).
% -------------------------------------------------------------------------

clear all;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modifiable parameters
file_name = '';      % Video file name. If empty, a prompt dialog will appear
sigma     = 25;      % Noise standard deviation. it should be in the same 
                     % intensity range of the video
profile = 'lc';      % V-BM4D parameter profile
                     %  'lc' --> low complexity
                     %  'np' --> normal profile
do_wiener = 1;       % Wiener filtering
                     %   1 --> enable Wiener filtering
                     %   0 --> disable Wiener filtering
sharpen = 1;         % Sharpening
                     %   1 --> disable sharpening
                     %  >1 --> enable sharpening
deflicker = 1;       % Deflickering
                     %   1 --> disable deflickering
                     %  <1 --> enable deflickering
verbose = 1;         % Verbose mode

est_noise = 0;       % Enable noise estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MODIFY BELOW THIS POINT ONLY IF YOU KNOW WHAT YOU ARE DOING       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If no file_name is given, let the user select a file from local disk
if isempty(file_name)
    [ file_name, folder ] = uigetfile({'*.avi'}, 'Select Video');
    if isequal(file_name,0)
        error('No file selected.')
    end
    file_name = [folder,file_name];
end
if strcmpi(file_name(end-2:end),'mat')
    % if a matlab file is given, the noise-free video is assumed to be
    % saved as a 3-D matrix called "y"
    load(file_name)
else
    % video formats mast be readable by VideoReader ('help VideoReader')
    y = read_video(file_name);
end

% Scaling data 
S = 255;
I_MAX = 255/S;
% Create synthetic Gaussian noise
randn('seed',0);
y = cast(y, 'single')/S;
sigma = sigma/S;
z = y + sigma*randn(size(y));
% Set sigma to -1 to enable noise estimation
if est_noise
    sigma = -1;
end

% V-BM4D filtering
disp('Denoising started')
tic;
y_est = vbm4d( z, sigma, profile, do_wiener, sharpen, deflicker, verbose );
time = toc;
PSNR = 10*log10(I_MAX^2/mean((y(:)-y_est(:)).^2));
fprintf('Denoising completed (%.1fs / %.1ffps): PSNR %.2fdB \n', time, size(z,ndims(z))/time, PSNR);

% Show results
f = round(size(y,3)/2);
figure,imshow([y(:,:,f) z(:,:,f) y_est(:,:,f)],[])



