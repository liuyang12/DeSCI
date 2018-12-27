function [ im_norm ] = imnorm( im_orig )
%IMNORM Maximum-minimum normalization of the 2-d image to the range of [0,1].

im_double = double(im_orig);
im_norm = (im_double - min(min(min(im_double))))...
          / (max(max(max(im_double))) - min(min(min(im_double))));

end

