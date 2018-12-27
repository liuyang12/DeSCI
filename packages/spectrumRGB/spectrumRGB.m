function sRGB = spectrumRGB(lambdaIn, varargin)
%spectrumRGB   Converts a spectral wavelength to RGB.
%
%    sRGB = spectrumRGB(LAMBDA) converts the spectral color with wavelength
%    LAMBDA in the visible light spectrum to its RGB values in the sRGB
%    color space.
%
%    sRGB = spectrumRGB(LAMBDA, MATCH) converts LAMBDA to sRGB using the
%    color matching function MATCH, a string.  See colorMatchFcn for a list
%    of supported matching functions.
%
%    Note: LAMBDA must be a scalar value or a vector of wavelengths.
%
%    See also colorMatchFcn, createSpectrum, spectrumLabel.

%    Copyright 1993-2005 The MathWorks, Inc.

if (numel(lambdaIn) ~= length(lambdaIn))
    
    error('spectrumRGB:unsupportedDimension', ...
          'Input must be a scalar or vector of wavelengths.')
    
end

% Get the color matching functions.
if (nargin == 1)
    
    matchingFcn = '1964_full';
    
elseif (nargin == 2)
    
    matchingFcn = varargin{1};
    
else
    
    error(nargchk(1, 2, nargin, 'struct'))
    
end

[lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn(matchingFcn);

% Interpolate the input wavelength in the color matching functions.
XYZ = interp1(lambdaMatch', [xFcn; yFcn; zFcn]', lambdaIn, 'pchip', 0);

% Reshape interpolated values to match standard image representation.
if (numel(lambdaIn) > 1)
    
    XYZ = permute(XYZ', [3 2 1]);
    
end

% Convert the XYZ values to sRGB.
XYZ2sRGB = makecform('xyz2srgb');
sRGB = applycform(XYZ, XYZ2sRGB);
