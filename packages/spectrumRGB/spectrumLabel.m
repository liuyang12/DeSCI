function varargout = spectrumLabel(hAxesTarget)
%spectrumLabel   Add a spectrum to the X-Axis of a plot.
%
%    H = SPECTRUMLABEL(HTARGET) adds a spectrum colorbar and wavelength
%    labels to the X axis of the axes object HTARGET.  For plots whose X
%    values are in the range 300 to 830, the colorbar will show the ROYGBIV
%    spectrum based on the 1964 CIE 10-degree observer.  Values outside the
%    range will be shown as black.  H contains a handle to the axes object
%    containing the spectrum.
%
%    Examples:
%    ---------
%    % Show the D65 illuminant.
%    [lambda, D65] = illuminant('d65');
%    figure;
%    plot(lambda, D65)
%    title('D_{65} illuminant')
%    spectrumLabel(gca)
%
%    % Show the CIE 1964 10-degree observer color matching functions.
%    [lambda, x, y, z] = colorMatchFcn('1964_full');
%    figure;
%    plot(lambda, x, 'r', lambda, y, 'g', lambda, z, 'b')
%    title('CIE 1964 10-degree observer')
%    spectrumLabel(gca)
%
%    See also colorMatchFcn, illuminant.

%    Copyright 1993-2005 The MathWorks, Inc.


hFig = ancestor(hAxesTarget, 'figure');

% Add a new axes to the same figure as the target axes.
figure(hFig);
hAxesSpectrum = axes;
set(hAxesSpectrum, 'visible', 'off')

% Areas outside the visible spectrum are black.
set(hAxesSpectrum, 'color', [0 0 0])

% Position the axes as appropriate.
targetPosition = get(hAxesTarget, 'position');
spectrumPosition = [targetPosition(1), ...
                    targetPosition(2), targetPosition(3), 1];
set(hAxesSpectrum, 'position', spectrumPosition)
set(hAxesSpectrum, 'units', 'pixels')

spectrumPosition = get(hAxesSpectrum, 'position');
set(hAxesSpectrum, 'position', [spectrumPosition(1), ...
                                spectrumPosition(2) - 20, ...
                                spectrumPosition(3), ...
                                20])

% Line the X limits of the two axes up and display the spectrum.
set(hAxesSpectrum, 'xtick', get(hAxesTarget, 'xtick'))

[lambda, RGB] = createSpectrum('1964_full');
axes(hAxesSpectrum);
image(lambda, 1:size(RGB,1), RGB);

% Turn off the unneeded axes labels.
set(hAxesTarget, 'xtick', []);
set(hAxesSpectrum, 'ytick', []);

% Make the figure visible.
set(hAxesSpectrum, 'units', 'normalized')
set(hAxesSpectrum, 'visible', 'on')

% Return the handle if requested.
if (nargout == 1)
    
    varargout{1} = hAxesSpectrum;
    
end

% Add a callback to update the label if the X limits of the target changes.
% TO DO.