function [Sxx, faxis] = getSpectrum(d, fs)
% Grabs the frequency spectrum of d, a given time series

% Adapted from Kramer 2013

dt = 1/fs;
T = size(d,2) / 500;

xf = fft(d);                        % Fourier transform of x
Sxx = (2*dt^2)/T * xf .* conj(xf);  % Compute the power spectrum
negVals = floor(length(d)/2+1);
Sxx = Sxx(1:negVals);               % Ignore negative frequencies


% Build the frequency axis
df = 1/max(T);                      % frequency resolution 
fNQ = 1/dt/2;                       % Nyquist frequency 
faxis = (0:df:fNQ);                

end


