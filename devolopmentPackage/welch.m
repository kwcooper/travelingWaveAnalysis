

function [Px, f] = welch(x, fs, win, overlap, Nfft, oneside)
% x=vector of data samples
% fs=sampling frequency (in Hz)
% win=vector containing the window
% overlap=number of points to overlap
% Nfft=length of FFT to use
% oneside=a logical variable indicating whether 
%    scaling is for 1-sided or 2-sided PSD (true if 1-sided)
% Px=estimate of the power spectral density
% f=vector of frequencies (in Hz) for plotting
if nargin < 6
  error('Insufficient input arguments')
end

x = x(:);
win = win(:);

Nw = length(win);
Ns = length(x);

% Compute the frequency vector and scale factor
if oneside == true
  sfactor = 2 ./ (fs * sum(win.^2));
  f = fs * (0:Nfft/2).'./Nfft;
else
  sfactor = 1 ./ (fs * sum(win.^2));
  f = fs * (-Nfft/2:Nfft/2-1).'./Nfft;
end

% Clip excess data
x = x(1:floor(Ns/(Nw-overlap)) * (Nw-overlap));

% Perform overlap by appending a shifted version
% of the data matrix
y = reshape(x, Nw-overlap, []);
y = [y; circshift(y(1:overlap, :), [0 -1])];
y = bsxfun(@times, y, win);

% get rid of last column
if size(y, 2) > 1
  y = y(:,1:end-1);
end

% Compute the PSD
%Px = sfactor * mean( abs(fft(y, Nfft)).^2, 2);
Px = sfactor * abs(fft(y, Nfft)).^2;

% Return full or onesided result
if oneside == true
  Px = Px(1:end/2+1);
else
  Px = fftshift(Px);
end