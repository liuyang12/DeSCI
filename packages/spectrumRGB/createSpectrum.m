function [lambda, RGB] = createSpectrum(matchingFcn)
%createSpectrum   Get the wavelengths and RGB values for a spectrum.
%
%    [LAMBDA, RGB] = CREATESPECTRUM(MATCHINGFCN) returns the wavelengths
%    LAMBDA (in nanometers) and RGB values for the visible spectrum.  The
%    string MATCHINGFCN determines which color matching functions are used
%    to generate the spectrum.
%
%    See also colorMatchFcn, spectrumRGB.

%    Copyright 1993-2005 The MathWorks, Inc.

% Get the color matching functions.
[lambda, xFcn, yFcn, zFcn] = colorMatchFcn(matchingFcn);

% For the equal-energy illuminant, the XYZ values are the same as the
% concatenated matching functions.
XYZ = permute([xFcn; yFcn; zFcn], [2 3 1]);

% Convert the XYZ values to RGB.
XYZ2sRGB = makecform('xyz2srgb');
RGB = applycform(XYZ, XYZ2sRGB);
RGB = repmat(RGB, [1, 30, 1]);
RGB = permute(RGB, [2 1 3]);
