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
% if isempty(file_name)
%     [ file_name, folder ] = uigetfile({'*.avi'}, 'Select Video');
%     if isequal(file_name,0)
%         error('No file selected.')
%     end
%     file_name = [folder,file_name];
% end
% if strcmpi(file_name(end-2:end),'mat')
%     % if a matlab file is given, the noise-free video is assumed to be
%     % saved as a 3-D matrix called "y"
%     load(file_name)
% else
%     % video formats mast be readable by VideoReader ('help VideoReader')
%     y = read_video(file_name);
% end
load('Traffic.mat');
y = X(:,:,1:12);
% Scaling data 
S = 1;
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
figure,
for f=1:size(y,3)
imshow([y(:,:,f) z(:,:,f) y_est(:,:,f)],[]); pause(0.1);
end