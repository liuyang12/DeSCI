function [psnrmean,psnrall] = vpsnr(v,ref,maxval)
%VPSNR Peak signal-to-noise ratio (PSNR) for a video or volume.
%   psnrmean=VPSNR(v,ref,maxval) returns the mean PSNR of the image
%   sequence, where v is the video or volume, ref is the reference 
%   video or volume with the same dimension as v, and maxval is the 
%   maximum value of the video or volume, say 255 for uint8.
%   See also: PSNR.
if nargin<3 % no maximum value
    maxval = [];
end
nframe = size(v,ndims(v));
psnrall = zeros([1 nframe]);
for iframe = 1:nframe
    if isempty(maxval)
        psnrall(iframe) = psnr(v(:,:,iframe),ref(:,:,iframe));
    else
        psnrall(iframe) = psnr(v(:,:,iframe),ref(:,:,iframe),maxval);
    end
end
psnrmean = mean(psnrall);

    