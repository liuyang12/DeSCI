function [ssimmean,ssimall] = vssim(v,ref)
%VSSIM Structure similarity index (SSIM) for a video or volume.
%   ssimmean=VSSIM(v,ref) returns the mean SSIM of the image sequence, 
%   where v is the video or volume and ref is the reference video or volume 
%   with the same dimension as v.
%   See also SSIM.
nframe = size(v,ndims(v));
ssimall = zeros([1 nframe]);
for iframe = 1:nframe
    ssimall(iframe) = ssim(v(:,:,iframe),ref(:,:,iframe));
end
ssimmean = mean(ssimall);

